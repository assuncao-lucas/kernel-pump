#ifndef FOMRULATIONS_H_
#define FORMULATIONS_H_

#include <vector>
#include <unordered_map>
#include <lemon/list_graph.h>
#include <lemon/preflow.h>
#include <ilcplex/ilocplex.h>
#include "src/instance.h"
#include "src/timer.h"
#include "src/solution.hpp"
#include "src/benders_generic_callback.h"
#include "src/matrix.hpp"
#include "src/general.h"
#include "src/heuristic_solution.h"

typedef lemon::ListDigraph LemonGraph;
typedef int LimitValueType;

int a_var_to_index(int vertex, int budget, int num_vertices);
int f_var_to_index(int arc_pos, int budget, int num_arcs);

std::pair<int, int> index_to_a_var(int index, int num_vertices);

template <class T>
void SetSolutionStatus(IloCplex &cplex, Solution<T> &solution, bool solve_relax)
{
  // std::cout << cplex.getCplexStatus() << std::endl;
  if (solve_relax)
  {
    if ((cplex.getCplexStatus() == IloCplex::Infeasible) || (cplex.getCplexStatus() == IloCplex::InfOrUnbd))
      solution.lp_ = -1;
    else
      solution.lp_ = cplex.getObjValue();
    // if((cplex.getCplexStatus() == IloCplex::Optimal)||(cplex.getCplexStatus() == IloCplex::OptimalTol))
    // {
    //   solution.is_optimal_ =  true;
    // }
  }
  else
  {
    if ((cplex.getCplexStatus() == IloCplex::Infeasible) || (cplex.getCplexStatus() == IloCplex::InfOrUnbd))
      solution.is_feasible_ = false;
    else
    {
      double cost = cplex.getObjValue();
      if ((cplex.getCplexStatus() == IloCplex::Optimal) || (cplex.getCplexStatus() == IloCplex::OptimalTol))
      {
        solution.lb_ = cost;
        solution.ub_ = cost;
        solution.is_optimal_ = true;
      }
      else
      {
        solution.lb_ = cost;
        solution.ub_ = cplex.getBestObjValue();
        if (!(double_less(solution.lb_, solution.ub_)))
          solution.is_optimal_ = true;
      }
    }
  }
  solution.num_nodes_ = cplex.getNnodes();
}
void SetPriorityOrder(IloCplex &cplex, IloEnv &env, IloNumVarArray &x, Instance &instance, double *R0);
static void addArcVertexInferenceCuts(IloCplex &cplex, IloModel &model, IloEnv &env, IloNumVarArray &x, IloNumVarArray &y, Instance &instance, bool solve_relax, Solution<double> &solution);
static bool FindAndAddValidInqualities(IloCplex &cplex, IloEnv &env, IloModel &model, IloNumArray &y_values, IloNumArray &x_values, Instance &instance, Solution<double> &sol, std::list<UserCutGeneral *> *root_cuts);
static UserCut *GenerateCliqueConflictCuts(Instance &instance, std::vector<bool> &visited_nodes, std::vector<double> &nodes_sum,
                                           IloNumArray &lps,
                                           std::list<int> &visited_nodes_list, std::unordered_map<int, int> &subgraph_to_graph_map,
                                           std::unordered_map<int, int> &graph_to_subgraph_map, LemonGraph &g,
                                           std::vector<LemonGraph::Node> &l_nodes, lemon::ListDigraph::ArcMap<LimitValueType> &capacity,
                                           LemonGraph &g_inv, std::vector<LemonGraph::Node> &l_nodes_inv,
                                           lemon::ListDigraph::ArcMap<LimitValueType> &capacity_inv,
                                           Solution<double> &sol, std::list<UserCutGeneral *> &cuts, bool root,
                                           boost::dynamic_bitset<> clique_is_active);
static bool ConflictIsActive(std::list<int> &curr_conflict, std::vector<double> &nodes_sum, double &curr_conflict_nodes_sum);

static void PopulateByRowCommon(IloEnv &env, IloModel &model, MasterVariables &master_vars, const Instance &instance, bool force_use_all_vehicles);
static void PopulateByRowCompactBaselineContinuousSpace(IloEnv &env, IloModel &model, IloNumVarArray &x, IloNumVarArray &y, IloNumVarArray &a, std::optional<std::reference_wrapper<IloNumArray>> x_values, std::optional<std::reference_wrapper<IloNumArray>> y_values, const Instance &instance);
void PopulateByRowCompactBaseline(IloCplex &cplex, IloEnv &env, IloModel &model, MasterVariables &master_vars, IloNumVarArray &a, const Instance &instance, double *R0, double *Rn, bool force_use_all_vehicles, bool export_model);

static void PopulateByRowCompactSingleCommodityContinuousSpace(IloEnv &env, IloModel &model, IloNumVarArray &x, IloNumVarArray &y, IloNumVarArray &f, std::optional<std::reference_wrapper<IloNumArray>> x_values, std::optional<std::reference_wrapper<IloNumArray>> y_values, double *R0, double *Rn, const Instance &instance);
void PopulateByRowCompactSingleCommodity(IloCplex &cplex, IloEnv &env, IloModel &model, MasterVariables &master_vars, IloNumVarArray &f, const Instance &instance, double *R0, double *Rn, bool force_use_all_vehicles, bool export_model);

void AllocateMasterVariablesBaseline(IloEnv &env, MasterVariables &master_vars, const Instance &instance, bool force_use_all_vehicles, bool solve_relax, bool disable_all_binary_vars = false);
void AllocateMasterVariablesSingleCommodity(IloEnv &env, MasterVariables &master_vars, const Instance &instance, bool force_use_all_vehicles, bool solve_relax, bool disable_all_binary_vars = false);

void Benders(Instance &inst, Formulation formumlation, double *R0, double *Rn, double time_limit, bool apply_benders_generic_callback, bool combine_feas_op_cuts, bool separate_benders_cuts_relaxation, bool solve_relax, bool use_valid_inequalities, bool find_root_cuts, std::list<UserCutGeneral *> *initial_cuts, HeuristicSolution *initial_sol, bool force_use_all_vehicles, bool export_model, std::list<UserCutGeneral *> *root_cuts, Solution<double> &);

void CompactBaseline(Instance &inst, double *R0, double *Rn, double time_limit, bool solve_relax, bool use_valid_inequalities, bool find_root_cuts, std::list<UserCutGeneral *> *initial_cuts, HeuristicSolution *initial_sol, bool force_use_all_vehicles, bool export_model, std::list<UserCutGeneral *> *root_cuts, Solution<double> &);
void PrimalSubproblemCompactBaseline(Instance &inst, IloNumArray &x_values, IloNumArray &y_values, double *R0, double *Rn, double time_limit, bool export_model, Solution<double> &solution);
void DualSubproblemCompactBaseline(Instance &inst, IloNumArray &x_values, IloNumArray &y_values, double *R0, double *Rn, double time_limit, bool export_model, Solution<double> &solution);

void CompactSingleCommodity(Instance &inst, double *R0, double *Rn, double time_limit, bool solve_relax, bool use_valid_inequalities, bool find_root_cuts, std::list<UserCutGeneral *> *initial_cuts, HeuristicSolution *initial_sol, bool force_use_all_vehicles, bool export_model, std::list<UserCutGeneral *> *root_cuts, Solution<double> &);
void PrimalSubproblemCompactSingleCommodity(Instance &inst, IloNumArray &x_values, IloNumArray &y_values, double *R0, double *Rn, double time_limit, bool export_model, Solution<double> &);
void DualSubproblemCompactSingleCommodity(Instance &inst, IloNumArray &x_values, IloNumArray &y_values, double *R0, double *Rn, double time_limit, bool export_model, Solution<double> &);

void optimize(IloCplex &cplex, IloEnv &env, IloModel &model, std::optional<Formulation> formulation, MasterVariables &master_vars, Instance &instance, double total_time_limit, bool solve_relax, bool apply_benders, bool apply_benders_generic_callback, bool combine_feas_op_cuts, bool separate_benders_cuts_relaxation, bool use_valid_inequalities, bool find_root_cuts, double *R0, double *Rn, std::list<UserCutGeneral *> *initial_cuts, HeuristicSolution *initial_sol, bool export_model, std::list<UserCutGeneral *> *root_cuts, Solution<double> &solution);
void optimizeLP(IloCplex &cplex, IloEnv &env, IloModel &model, Instance &instance, double total_time_limit, double *R0, double *Rn, Solution<double> &);

#endif