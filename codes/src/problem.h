#pragma once

#include <iostream>
#include <string>
#include <ilcplex/ilocplex.h>
#include <boost/dynamic_bitset.hpp>
#include <optional>

class Problem final
{
public:
    explicit Problem(std::string file_path, bool multithreading);
    virtual ~Problem();

    bool set_multithreading(bool multithreading);
    int num_vars() const { return num_vars_; }
    int num_binary_vars() const { return num_binary_vars_; }
    const boost::dynamic_bitset<> &binary_vars() const { return binary_vars_; }
    const boost::dynamic_bitset<> &curr_active_binary_vars() const { return curr_active_binary_vars_; }
    IloEnv *env() { return env_; }
    IloNumVarArray &vars() { return *vars_; }
    IloCplex::CplexStatus curr_status() { return curr_status_; }
    std::optional<double> curr_obj_value() { return curr_obj_value_; }
    bool IsMinimization() { return original_sense_ == IloObjective::Sense::Minimize; }
    bool Solve(bool solve_relaxed, bool stop_at_first_feas_integer_solution);
    void GetValues(IloNumArray &);
    void GetReducedCosts(IloNumArray &);
    void DeactivateAllBinaryVariables();
    void UpdateModelVarBounds(std::optional<boost::dynamic_bitset<>> vars_entering_problem, std::optional<boost::dynamic_bitset<>> vars_leaving_problem);
    void SetObjectiveExpression(IloExpr new_exp)
    {
        obj_->setExpr(new_exp);
    }

    void SetObjectiveSense(IloObjective::Sense sense)
    {
        obj_->setSense(sense);
    }

    IloExpr &original_obj_expr() { return original_obj_expr_; }
    double original_obj_norm() { return original_obj_norm_; }
    std::optional<double> ComputeSolutionValue(boost::dynamic_bitset<> solution);
    void Reset();

private:
    void InitCplex();
    void ResetCplex();
    void BuildModel(bool linearly_relaxed, bool export_model);
    void SetVarUB(int var_index, double ub) { (*vars_)[var_index].setUB(ub); }
    void SetVarLB(int var_index, double lb) { (*vars_)[var_index].setLB(lb); }

    IloEnv *env_ = nullptr;             // Cplex environment
    IloCplex *cplex_ = nullptr;         // Cplex solver
    IloModel *model_ = nullptr;         // Cplex model
    IloModel *relaxed_model_ = nullptr; // Cplex model
    IloObjective *obj_ = nullptr;
    IloNumVarArray *vars_ = nullptr;
    IloRangeArray *rng_ = nullptr;
    IloExpr original_obj_expr_;
    IloObjective::Sense original_sense_ = IloObjective::Sense::Minimize;
    double original_obj_norm_ = 0.0;
    IloCplex::CplexStatus curr_status_ = IloCplex::Unknown;
    std::optional<double> curr_obj_value_ = std::nullopt;

    int num_vars_ = 0;
    int num_binary_vars_ = 0;
    boost::dynamic_bitset<> binary_vars_;
    boost::dynamic_bitset<> curr_active_binary_vars_;
    std::string problem_;
    bool multithreading_ = false;
};