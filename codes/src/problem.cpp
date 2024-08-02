#include "src/problem.h"
#include "src/general.h"

Problem::Problem(std::string file_path, bool multithreading) : problem_(file_path), multithreading_(multithreading)
{
    Reset();
}

Problem::~Problem()
{
    ResetCplex();
}

void Problem::InitCplex()
{
    ResetCplex();
    env_ = new IloEnv();
    model_ = new IloModel(*env_);
    cplex_ = new IloCplex(*env_);
    obj_ = new IloObjective(*env_);
    vars_ = new IloNumVarArray(*env_);
    rng_ = new IloRangeArray(*env_);
    cplex_->setOut(env_->getNullStream());
    cplex_->setWarning(env_->getNullStream());
    cplex_->setError(env_->getNullStream());
    cplex_->extract(*model_);
    curr_status_ = IloCplex::Unknown;
    curr_obj_value_ = std::nullopt;
    set_multithreading(multithreading_);
}

void Problem::ResetCplex()
{
    if (cplex_)
    {
        // cplex_->end();
        delete cplex_;
        cplex_ = nullptr;
    }

    if (model_)
    {
        // model_->end();
        delete model_;
        model_ = nullptr;
    }

    if (relaxed_model_)
    {
        // relaxed_model_->end();
        delete relaxed_model_;
        relaxed_model_ = nullptr;
    }

    if (obj_)
    {
        // obj_->end();
        delete obj_;
        obj_ = nullptr;
    }
    if (vars_)
    {
        // vars_->endElements();
        // vars_->end();
        delete vars_;
        vars_ = nullptr;
    }
    if (rng_)
    {
        // rng_->endElements();
        // rng_->end();
        delete rng_;
        rng_ = nullptr;
    }
    if (env_)
    {
        env_->end();
        delete env_;
        env_ = nullptr;
    }
    num_vars_ = 0;
    num_binary_vars_ = 0;
}

void Problem::BuildModel(bool linearly_relaxed, bool export_model)
{
    cplex_->importModel(*model_, problem_.c_str(), *obj_, *vars_, *rng_);
    num_vars_ = vars_->getSize();

    binary_vars_ = boost::dynamic_bitset<>(num_vars_, 0);

    // set binary variables.
    for (int i = 0; i < num_vars_; ++i)
    {
        // std::cout << i << " " << (*vars_)[i] << " " << (*vars_)[i].getType() << std::endl;
        if ((*vars_)[i].getType() == IloNumVar::Type::Bool)
            binary_vars_[i] = 1;
    }

    // all binary vars are active at first!
    curr_active_binary_vars_ = binary_vars_;

    num_binary_vars_ = binary_vars_.count();
    std::cout << num_binary_vars_ << " / " << num_vars_ << std::endl;

    if (linearly_relaxed)
    {
        relaxed_model_ = new IloModel(*env_);
        relaxed_model_->add(*model_);
        relaxed_model_->add(IloConversion(*env_, *vars_, ILOFLOAT));
    }

    if (export_model)
        cplex_->exportModel("problem.lp");
}

bool Problem::set_multithreading(bool multithreading)
{
    multithreading_ = multithreading;
    if (cplex_)
    {
        multithreading ? cplex_->setParam(IloCplex::Param::Threads, cplex_->getNumCores()) : cplex_->setParam(IloCplex::Param::Threads, 1);
        return true;
    }
    return false;
}

void Problem::Reset()
{
    InitCplex();
    BuildModel(true, false);
    original_sense_ = obj_->getSense();
    original_obj_expr_ = obj_->getExpr();

    for (IloExpr::LinearIterator it = original_obj_expr_.getLinearIterator(); it.ok(); ++it)
    {
        original_obj_norm_ += pow(1.0 * it.getCoef(), 2.0);
        // std::cout << "variable: " << it.getVar() << " coef: " << it.getCoef() << std::endl;
    }
    // std::cout << original_obj_norm_ << std::endl;
    // std::cout << original_obj_expr_ << std::endl;
    original_obj_norm_ = sqrt(original_obj_norm_);
    // std::cout << "original_obj_norm: " << original_obj_norm_ << std::endl;
}

bool Problem::Solve(bool solve_relaxed, bool stop_at_first_feas_integer_solution)
{
    bool result = false;
    if (cplex_)
    {
        if (solve_relaxed)
        {
            if (relaxed_model_)
            {
                cplex_->extract(*relaxed_model_);
                result = cplex_->solve();
                curr_status_ = cplex_->getCplexStatus();
                if (curr_status_ == IloCplex::Status::Optimal)
                    curr_obj_value_ = std::optional<double>(cplex_->getObjValue());
                else
                    curr_obj_value_ = std::nullopt;
                return result;
            }
        }
        else if (model_)
        {
            if (stop_at_first_feas_integer_solution)
            {
                cplex_->setParam(IloCplex::Param::MIP::Limits::Solutions, 1);
            }
            cplex_->extract(*model_);
            result = cplex_->solve();
            curr_status_ = cplex_->getCplexStatus();
            curr_obj_value_ = std::nullopt;
            if (curr_status_ != IloCplex::Status::Infeasible && curr_status_ != IloCplex::Status::InfOrUnbd)
            {
                curr_obj_value_ = std::optional<double>(cplex_->getObjValue());
                // IloNumArray values(*env_);
                // cplex_->getValues(values, *vars_);
                // std::cout << "int sol:" << std::endl;
                // auto sol = boost::dynamic_bitset<>(num_vars_, 0);
                // for (int var_index = binary_vars_.find_first(); var_index != boost::dynamic_bitset<>::npos; var_index = binary_vars_.find_next(var_index))
                // {
                //     if (double_equals(values[var_index], 1))
                //     {
                //         sol[var_index] = 1;
                //         std::cout << var_index << " ";
                //     }
                // }
                // std::cout << std::endl;
                // values.end();

                // std::cout << "cost: " << *curr_obj_value_ << std::endl;
            }
            else
            {
                curr_obj_value_ = std::nullopt;
            }
            return result;
        }
    }
    return result;
}

std::optional<double> Problem::ComputeSolutionValue(boost::dynamic_bitset<> solution)
{
    for (int var_index = binary_vars_.find_first(); var_index != boost::dynamic_bitset<>::npos; var_index = binary_vars_.find_next(var_index))
    {
        SetVarLB(var_index, solution[var_index]);
        SetVarUB(var_index, solution[var_index]);
    }

    IloExpr current_obj_expr = obj_->getExpr();
    SetObjectiveExpression(original_obj_expr_);

    Solve(true, false);

    // after computing, restore to original state of active variables.
    SetObjectiveExpression(current_obj_expr);
    for (int var_index = binary_vars_.find_first(); var_index != boost::dynamic_bitset<>::npos; var_index = binary_vars_.find_next(var_index))
    {
        SetVarLB(var_index, 0);
        SetVarUB(var_index, curr_active_binary_vars_[var_index]);
    }

    return curr_obj_value_;
}

void Problem::GetValues(IloNumArray &values)
{
    cplex_->getValues(values, *vars_);
}

void Problem::GetReducedCosts(IloNumArray &red_costs)
{
    cplex_->getReducedCosts(red_costs, *vars_);
}

void Problem::UpdateModelVarBounds(std::optional<boost::dynamic_bitset<>> vars_entering_problem, std::optional<boost::dynamic_bitset<>> vars_leaving_problem)
{
    // activate new variables to kernel.
    if (vars_entering_problem)
    {
        // std::cout << "entering: ";
        for (int var_index = vars_entering_problem->find_first(); var_index != boost::dynamic_bitset<>::npos; var_index = vars_entering_problem->find_next(var_index))
        {
            SetVarUB(var_index, 1.0);
            curr_active_binary_vars_[var_index] = 1;
            // std::cout << var_index << " ";
        }
        // std::cout << std::endl;
    }

    // deactivate variables leaving kernel.
    if (vars_leaving_problem)
    {
        // std::cout << "leaving: ";
        for (int var_index = vars_leaving_problem->find_first(); var_index != boost::dynamic_bitset<>::npos; var_index = vars_leaving_problem->find_next(var_index))
        {
            SetVarUB(var_index, 0.0);
            curr_active_binary_vars_[var_index] = 0;
            // std::cout << var_index << " ";
        }
        // std::cout << std::endl;
    }
}

void Problem::DeactivateAllBinaryVariables()
{
    UpdateModelVarBounds(std::nullopt, binary_vars_);
}