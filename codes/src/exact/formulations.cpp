#include <ilcplex/ilocplex.h>
#include "src/instance.h"
#include "src/general.h"
#include "src/exact/formulations.h"
#include "src/user_cut.h"
#include "src/initial_solution/initial_solution.h"

ILOSTLBEGIN

ILOUSERCUTCALLBACK7(BendersUserCutCallback, IloCplex &, worker_cplex, IloObjective &, worker_obj,
                    DualVariablesBaseline *, dual_vars, MasterVariables &, master_vars, Instance &, instance,
                    bool, combine_feas_op_cuts, Solution<double> &, solution)
{
  // only separate user cuts ad root node
  if (getNnodes() > 0)
    return;
  // Skip the separation if not at the end of the cut loop
  if (!isAfterCutLoop())
    return;

  IloEnv master_env = getEnv();
  // Get the current solution.
  IloNumArray x_values(master_env);
  IloNumArray y_values(master_env);
  IloNum dual_bound_value = getValue(master_vars.dual_bound);
  getValues(y_values, master_vars.y);
  getValues(x_values, master_vars.x);

  // Benders' cut separation.
  IloExpr cut_expr(master_env);
  IloBool sep_status = SeparateBendersCutBaseline(master_env, master_vars.x, master_vars.y, master_vars.dual_bound, x_values, y_values, dual_bound_value, worker_cplex, dual_vars, instance, worker_obj, cut_expr, combine_feas_op_cuts, solution);
  if (sep_status)
  {
    // std::cout << "found new Cut" << std::endl;
    add(cut_expr >= 0).end();
  }

  // Free memory.
  cut_expr.end();
  x_values.end();
  y_values.end();

  return;

} // END BendersUserCutCallback

// Implementation class for the user-defined lazy constraint callback.
// The function BendersLazyCallback allows to add Benders' cuts as lazy constraints.
//
ILOLAZYCONSTRAINTCALLBACK7(BendersLazyCallbackBaseline, IloCplex &, worker_cplex, IloObjective &, worker_obj,
                           DualVariablesBaseline *, dual_vars, MasterVariables &, master_vars, Instance &, instance,
                           bool, combine_feas_op_cuts, Solution<double> &, solution)
{
  IloEnv master_env = getEnv();
  // Get the current solution.
  IloNumArray x_values(master_env);
  IloNumArray y_values(master_env);
  IloNum dual_bound_value = getValue(master_vars.dual_bound);
  getValues(y_values, master_vars.y);
  getValues(x_values, master_vars.x);

  // std::cout << "current master obj " << getObjValue() << std::endl;
  // std::cout << "current dual bound " << dual_bound_value << std::endl;

  // Benders' cut separation.
  IloExpr cut_expr(master_env);
  IloBool sep_status = SeparateBendersCutBaseline(master_env, master_vars.x, master_vars.y, master_vars.dual_bound, x_values, y_values, dual_bound_value, worker_cplex, dual_vars, instance, worker_obj, cut_expr, combine_feas_op_cuts, solution);
  if (sep_status)
  {
    // std::cout << "found new Cut" << std::endl;
    add(cut_expr >= 0).end();
  }

  // Free memory.
  cut_expr.end();
  x_values.end();
  y_values.end();

  return;

} // END BendersLazyCallback

struct CallbackArguments
{
  const double *R0;
  const double *Rn;
  const Instance *instance;
};

ILOLAZYCONSTRAINTCALLBACK7(BendersLazyCallbackSingleCommodity, IloCplex &, worker_cplex, IloObjective &, worker_obj,
                           DualVariablesSingleCommodity *, dual_vars, MasterVariables &, master_vars, CallbackArguments &, arguments,
                           bool, combine_feas_op_cuts, Solution<double> &, solution)
{
  IloEnv master_env = getEnv();
  // Get the current solution.
  IloNumArray x_values(master_env);
  IloNumArray y_values(master_env);
  IloNum dual_bound_value = getValue(master_vars.dual_bound);

  // std::cout << "current master obj " << getObjValue() << std::endl;
  // std::cout << "current dual bound " << dual_bound_value << std::endl;
  getValues(y_values, master_vars.y);
  getValues(x_values, master_vars.x);

  // for (int i =0 ; i < y_values.getSize(); ++i)
  //   std::cout << master_vars.y[i].getName() << " " << y_values[i] << std::endl;

  // for (int i =0 ; i < x_values.getSize(); ++i)
  //   std::cout << master_vars.x[i].getName() << " " << x_values[i] << std::endl;

  // Benders' cut separation.
  IloExpr cut_expr(master_env);
  IloBool sep_status = SeparateBendersCutSingleCommodity(master_env, master_vars.x, master_vars.y, master_vars.dual_bound, x_values, y_values, dual_bound_value, worker_cplex, dual_vars, arguments.R0, arguments.Rn, *(arguments.instance), worker_obj, cut_expr, combine_feas_op_cuts, solution);
  if (sep_status)
  {
    // std::cout << "found new Cut" << std::endl;
    add(cut_expr >= 0).end();
  }

  // Free memory.
  cut_expr.end();
  x_values.end();
  y_values.end();

  return;

} // END BendersLazyCallback

static bool ConflictIsActive(std::list<int> &curr_conflict, std::vector<double> &nodes_sum, double &curr_conflict_nodes_sum)
{
  int cont = 0;
  curr_conflict_nodes_sum = 0.0;
  for (std::list<int>::iterator it = curr_conflict.begin(); it != curr_conflict.end(); ++it)
  {
    if (double_greater(nodes_sum[*it], 0.0))
    {
      ++cont;
      curr_conflict_nodes_sum += (nodes_sum[*it]);
    }
  }
  return (cont >= 1);
}

static UserCut *GenerateCliqueConflictCuts(Instance &instance, std::vector<bool> &visited_nodes, std::vector<double> &nodes_sum,
                                           IloNumArray &lps,
                                           std::list<int> &visited_nodes_list, std::unordered_map<int, int> &subgraph_to_graph_map,
                                           std::unordered_map<int, int> &graph_to_subgraph_map, LemonGraph &g,
                                           std::vector<LemonGraph::Node> &l_nodes, lemon::ListDigraph::ArcMap<LimitValueType> &capacity,
                                           LemonGraph &g_inv, std::vector<LemonGraph::Node> &l_nodes_inv,
                                           lemon::ListDigraph::ArcMap<LimitValueType> &capacity_inv,
                                           Solution<double> &sol, std::list<UserCutGeneral *> &cuts, bool root,
                                           boost::dynamic_bitset<> clique_is_active)
{
  const Graph *graph = instance.graph();
  UserCut *best_cut = nullptr, *curr_cut = nullptr;
  double curr_max_flow = 0.0;
  bool solve_from_source = true, solve_from_destination = true;
  const auto &map_vertices_to_cliques = instance.map_vertices_to_cliques();

  /*std::cout << "G" << std::endl;
  std::cout << "Nodes:";
  for (lemon::ListDigraph::NodeIt i(g); i!=lemon::INVALID; ++i)
  std::cout << " " << g.id(i);
  std::cout << std::endl;

  std::cout << "Arcs:";
  for (lemon::ListDigraph::ArcIt i(g); i!=lemon::INVALID; ++i)
  std::cout << " (" << g.id(g.source(i)) << "," << g.id(g.target(i)) << ")[" << capacity[i] << "]";
  std::cout << std::endl;

  std::cout << "G inv" << std::endl;
  std::cout << "Nodes:";
  for (lemon::ListDigraph::NodeIt i(g_inv); i!=lemon::INVALID; ++i)
  std::cout << " " << g_inv.id(i);
  std::cout << std::endl;

  std::cout << "Arcs:";
  for (lemon::ListDigraph::ArcIt i(g_inv); i!=lemon::INVALID; ++i)
  std::cout << " (" << g_inv.id(g_inv.source(i)) << "," << g_inv.id(g_inv.target(i)) << ")[" << capacity_inv[i] << "]";
  std::cout << std::endl;

  getchar();getchar();*/

  int v3 = 0, v4 = 0;

  int num_vertices = graph->num_vertices();
  int num_arcs = graph->num_arcs();
  int num_vehicles = instance.num_vehicles();
  // std::list<std::pair<std::pair<int,int>,double>> visited_arcs;

  LemonGraph::Node artificial_node = g.addNode(); // artificial node of conflict tuple (of nodes)
  // LemonGraph::Arc arc_artifical_1, arc_artifical_2; // one from/to each node in conflict pair

  LemonGraph::Node artificial_node_inv = g_inv.addNode(); // artificial node of conflict tuple (of nodes)
  // LemonGraph::Arc arc_artifical_1_inv, arc_artifical_2_inv; // one from/to each node in conflict pair

  double curr_conflict_nodes_sum = 0.0;
  int clique_counter = -1;
  // boost::dynamic_bitset<> clique_is_active((instance.conflicts_list_).size());
  // clique_is_active.set();
  // std::vector<bool> clique_is_active((instance.conflicts_list_).size(),true);
  // std::cout << "init search" << std::endl;
  // std::cout << "num active cliques: " << clique_is_active.count() << "/" << clique_is_active.size() << std::endl;
  bool found_from_source = false;

  auto conflicts_list = instance.conflicts_list();
  for (std::vector<std::list<int>>::iterator it = conflicts_list.begin(); it != conflicts_list.end(); ++it)
  {
    ++clique_counter;
    if (((*it).size() == 1) && ((*it).front() == 0))
      clique_is_active[clique_counter] = false;
    if (clique_is_active[clique_counter])
    {
      std::list<int> curr_conflict = *it;
      if (ConflictIsActive(curr_conflict, nodes_sum, curr_conflict_nodes_sum)) // if at least one conflicting vertex is in the solution
      {
        if (solve_from_source)
        {
          found_from_source = false;
          // curr_conflict_nodes_sum = nodes_sum[v1] + nodes_sum[v2];
          std::list<LemonGraph::Arc> artificial_arcs;
          // int cont = 0;
          for (std::list<int>::iterator it2 = curr_conflict.begin(); it2 != curr_conflict.end(); ++it2)
          {
            if (visited_nodes[*it2])
            {
              artificial_arcs.push_back(g.addArc(l_nodes[graph_to_subgraph_map[*it2]], artificial_node));
              capacity[artificial_arcs.back()] = (int)(round(num_vehicles * K_PRECISION)); // it's sufficiently big by the definition of the problem (model)
            }
          }

          lemon::Preflow<LemonGraph> preflow(g, capacity, l_nodes[graph_to_subgraph_map[0]], artificial_node);
          preflow.init();
          preflow.startFirstPhase();

          curr_max_flow = preflow.flowValue() / K_PRECISION;

          if (double_less(curr_max_flow, 0.0))
          {
            // std::cout << "negativo no source!!" << std::endl;
            // getchar(); getchar();
            for (auto it3 = artificial_arcs.begin(); it3 != artificial_arcs.end(); ++it3)
              g.erase(*it3);
            continue;
          }

          // if found a violated cut and conflict nodes do not belong to this cut on the source side.
          if ((double_less(curr_max_flow, curr_conflict_nodes_sum, K_CLIQUE_CONFLICT_TOLERANCE))
              //&& (!(preflow.minCut(l_nodes[graph_to_subgraph_map[v1]])))
              //&& (!(preflow.minCut(l_nodes[graph_to_subgraph_map[v2]])))
              //&& (!(preflow.minCut(artificial_node)))
          )
          {
            found_from_source = true;
            // std::cout << "achou CCC: " << curr_max_flow << " " << curr_conflict_nodes_sum << std::endl;
            curr_cut = new UserCut(num_arcs, num_vertices, curr_conflict_nodes_sum - curr_max_flow, K_TYPE_CLIQUE_CONFLICT_CUT);
            // curr_cut->rhs_nonzero_coefficients_indexes_.push_back(v1);
            // curr_cut->rhs_nonzero_coefficients_indexes_.push_back(v2);
            for (std::list<int>::iterator it2 = curr_conflict.begin(); it2 != curr_conflict.end(); ++it2)
            {
              // std::cout << *it2 << std::endl;
              curr_cut->AddRhsElement(*it2);
              if (visited_nodes[*it2])
              {
                // curr_cut->AddRhsElement(*it2);
                // disable cliques which contain vertices of the current violated clique with non-zero relaxation values
                // std::cout << "before #active cliques: " << clique_is_active.count() << std::endl;
                for (auto it3 = (map_vertices_to_cliques[*it2]).begin(); it3 != (map_vertices_to_cliques[*it2]).end(); ++it3)
                  clique_is_active[*it3] = 0;
                // std::cout << "# cliques of vertex " << *it2 <<  ": " << ((instance.map_vertices_to_cliques_)[*it2]).size() << std::endl;
                // std::cout << "after #active cliques: " << clique_is_active.count() << std::endl;
                // getchar();getchar();
              }
            }
            // std::cout << "Achou violação 1!!" << std::endl;

            // std::cout << "** " << curr_max_flow << " < " << max_node_sum << "==" << nodes_sum[max_node_sum_index] << std::endl;
            //  Generates new constraint
            // IloExpr exp(env);
            // IloExpr exp2(env);
            // double sum1 = 0.0;
            // std::cout << "S(lemon): ";
            for (v3 = 0; v3 < num_vertices; v3++)
            {
              auto &v3_out_vertices = graph->AdjVerticesOut(v3);
              for (std::list<int>::iterator it3 = v3_out_vertices.begin(); it3 != v3_out_vertices.end(); ++it3)
              {
                v4 = *it3;
                // if the arc traverses the cut
                if (((visited_nodes[v3]) && (preflow.minCut(l_nodes[graph_to_subgraph_map[v3]]))) &&
                    (!(visited_nodes[v4]) || (!(preflow.minCut(l_nodes[graph_to_subgraph_map[v4]])))))
                {
                  //(curr_cut->lhs_nonzero_coefficients_indexes_).push_back(std::pair<int,int>(v3,v4));
                  //(curr_cut->lhs_coefficients_)[graph->pos(v3,v4)] = 1;
                  curr_cut->AddLhsElement(v3, v4, graph->pos(v3, v4));
                  // std::cout << "add " << v3 << ", " << v4 << std::endl;
                  // sum1 += lps[graph->pos(v3,v4)];
                  //(curr_cut->lhs_num_nonzero_coefs_)++;
                }
              }
            }
            // std::cout << "* " << max_node_sum_index << std::endl;
            // std::cout << sum1 << " == " << curr_max_flow << " < " << curr_conflict_nodes_sum << std::endl;
            // getchar();getchar();
            // add(exp >= exp2).end();
            //(sol.num_cuts_added_)++;

            cuts.push_back(curr_cut);
            sol.set_cut_found(curr_cut->type_, root);
            curr_cut->UpdateMeasures();

            if (curr_cut->isBetterThan(best_cut))
              best_cut = curr_cut;

            // std::cout << sol.num_cuts_added_ << std::endl;

            // exp.end();
            // exp2.end();
          }

          for (auto it3 = artificial_arcs.begin(); it3 != artificial_arcs.end(); ++it3)
            g.erase(*it3);
        }

        if (solve_from_destination && !found_from_source)
        {

          std::list<LemonGraph::Arc> artificial_arcs_inv;
          for (std::list<int>::iterator it2 = curr_conflict.begin(); it2 != curr_conflict.end(); ++it2)
          {
            if (visited_nodes[*it2])
            {
              artificial_arcs_inv.push_back(g_inv.addArc(l_nodes_inv[graph_to_subgraph_map[*it2]], artificial_node_inv));
              capacity_inv[artificial_arcs_inv.back()] = (int)(round(num_vehicles * K_PRECISION)); // it's sufficiently big by the definition of the problem (model)
            }
          }

          lemon::Preflow<LemonGraph> preflow(g_inv, capacity_inv, l_nodes_inv[graph_to_subgraph_map[0]], artificial_node_inv);
          preflow.init();
          preflow.startFirstPhase();

          curr_max_flow = preflow.flowValue() / K_PRECISION;

          if (double_less(curr_max_flow, 0.0))
          {
            // std::cout << "negativo no sink!!" << std::endl;
            for (auto it3 = artificial_arcs_inv.begin(); it3 != artificial_arcs_inv.end(); ++it3)
              g_inv.erase(*it3);
            continue;
          }
          // if found a violated cut and conflict nodes and source belong to this cut on the conflict artificial node side.

          if ((double_less(curr_max_flow, curr_conflict_nodes_sum, K_CLIQUE_CONFLICT_TOLERANCE))
              //  && (!(preflow.minCut(l_nodes_inv[graph_to_subgraph_map[v1]])))
              //  && (!(preflow.minCut(l_nodes_inv[graph_to_subgraph_map[v2]])))
              //  && (!(preflow.minCut(artificial_node_inv)))
          )
          {
            // std::cout << "achou!!!" << std::endl;
            // std::cout << "achou CCC 2: " << curr_max_flow << " " << curr_conflict_nodes_sum << std::endl;
            // std::cout << curr_max_flow << " < " << curr_conflict_nodes_sum << std::endl;
            curr_cut = new UserCut(num_arcs, num_vertices, curr_conflict_nodes_sum - curr_max_flow, K_TYPE_CLIQUE_CONFLICT_CUT);
            // curr_cut->rhs_nonzero_coefficients_indexes_.push_back(v1);
            // curr_cut->rhs_nonzero_coefficients_indexes_.push_back(v2);
            for (std::list<int>::iterator it2 = curr_conflict.begin(); it2 != curr_conflict.end(); ++it2)
            {
              // std::cout << *it2 << std::endl;
              curr_cut->AddRhsElement(*it2);
              if (visited_nodes[*it2])
              {
                // curr_cut->AddRhsElement(*it2);
                // disable cliques which contain vertices of the current violated clique with non-zero relaxation values
                for (auto it3 = (map_vertices_to_cliques[*it2]).begin(); it3 != ((map_vertices_to_cliques)[*it2]).end(); ++it3)
                  clique_is_active[*it3] = 0;
              }
            }
            // std::cout << "Achou violação 2!!" << std::endl;
            // std::cout << "max flow from sink: " << curr_max_flow << std::endl;
            // std::cout << "conflict nodes sum: " << curr_conflict_nodes_sum << std::endl;
            // std::cout << "v1: " << v1 << std::endl;
            // std::cout << "v2: " << v2 << std::endl;
            // std::cout << "desvio: " << curr_conflict_nodes_sum - curr_max_flow << std::endl;
            // getchar();getchar();
            // getchar();getchar();
            // IloExpr exp(env);
            // IloExpr exp2(env);
            // double sum1 = 0, sum2 = 0, sum3 =0;
            // std::cout << "S(lemon): ";
            for (v3 = 0; v3 < num_vertices; ++v3)
            {
              auto &out_vertices = graph->AdjVerticesOut(v3);
              for (std::list<int>::iterator it3 = out_vertices.begin(); it3 != out_vertices.end(); ++it3)
              {
                v4 = *it3;
                // if the arc traverses the cut
                if ((!(visited_nodes[v3]) || !(preflow.minCut(l_nodes_inv[graph_to_subgraph_map[v3]]))) &&
                    (((visited_nodes[v4]) && (preflow.minCut(l_nodes_inv[graph_to_subgraph_map[v4]])))))
                {
                  //(curr_cut->lhs_nonzero_coefficients_indexes_).push_back(std::pair<int,int>(v3,v4));
                  //(curr_cut->lhs_coefficients_)[graph->pos(v3,v4)] = 1;
                  curr_cut->AddLhsElement(v3, v4, graph->pos(v3, v4));
                  //(curr_cut->lhs_num_nonzero_coefs_)++;
                }
              }
            }

            // std::cout << "* " << max_node_sum_index << std::endl;
            // std::cout << sum1 << "==" << curr_max_flow << ">=" << sum2 << "==" << sum3 << "==" << nodes_sum[max_node_sum_index] << "==" << max_node_sum << std::endl;
            // getchar();getchar();
            // add(exp >= exp2).end();
            //(sol.num_cuts_added_)++;

            cuts.push_back(curr_cut);
            sol.set_cut_found(curr_cut->type_, root);

            curr_cut->UpdateMeasures();
            if (curr_cut->isBetterThan(best_cut))
              best_cut = curr_cut;

            // std::cout << sol.num_cuts_added_ << std::endl;

            // exp.end();
            // exp2.end();
          }

          for (auto it3 = artificial_arcs_inv.begin(); it3 != artificial_arcs_inv.end(); ++it3)
            g_inv.erase(*it3);
        }
      }
    }
  }

  g.erase(artificial_node);
  g_inv.erase(artificial_node_inv);

  return best_cut;
}

static bool FindAndAddValidInqualities(IloCplex &cplex, IloEnv &env, IloModel &model, IloNumVarArray &y, IloNumVarArray &x, IloNumArray &values_y, IloNumArray &values_x, Instance &instance, Solution<double> &sol, std::list<UserCutGeneral *> *root_cuts)
{
  // std::cout << " ***************** looking for cuts" << std::endl;
  std::vector<bool> *CALLBACKS_SELECTION = GetCallbackSelection();
  bool found_cut = false;
  (sol.num_calls_to_callback_lp_) += 1;

  std::list<UserCutGeneral *> cuts;

  // std::cout << "=>>> CALLBACK: ";
  // std::cout << (sol.num_calls_to_callback_) << std::endl;

  const Graph *graph = instance.graph();

  double curr_value = 0.0;
  int v1 = 0, v2 = 0;
  // int curr_vertex = -1, previous_vertex = -1;

  UserCut *best_cut = nullptr, *curr_cut = nullptr, *local_best_cut = nullptr;

  const int num_vertices = graph->num_vertices(), num_vehicles = instance.num_vehicles();
  std::vector<bool> visited_nodes(num_vertices, false);
  std::vector<double> nodes_sum(num_vertices, 0.0);
  std::unordered_map<int, int> subgraph_to_graph_map;
  std::unordered_map<int, int> graph_to_subgraph_map;
  nodes_sum[0] = 1.0;
  std::list<int> visited_nodes_list;
  // std::list<std::pair<std::pair<int,int>,double>> visited_arcs;

  // std::list<std::list<int>> paths_list;
  std::list<std::pair<int, int>> visited_arcs;
  // std::vector<int> parent;

  // if((*CALLBACKS_SELECTION)[K_TYPE_PATH_BOUND_CUT]) parent = std::vector<int>(num_vertices,-1);

  // marks all nodes and arcs that must be considered in the max flow algorithm
  std::list<int> q;
  q.push_front(0);
  visited_nodes[0] = true;
  visited_nodes_list.push_front(0);
  subgraph_to_graph_map[0] = 0;
  graph_to_subgraph_map[0] = 0;
  int cont = 1;
  LemonGraph g, g_inv;
  std::vector<LemonGraph::Node> l_nodes, l_nodes_inv;
  lemon::ListDigraph::ArcMap<LimitValueType> capacity(g), capacity_inv(g_inv);
  l_nodes.push_back(g.addNode());         // node 0
  l_nodes_inv.push_back(g_inv.addNode()); // node 0

  do
  {
    v1 = q.front();
    q.pop_front();

    auto &v1_out_vertices = graph->AdjVerticesOut(v1);
    for (std::list<int>::iterator it = v1_out_vertices.begin(); it != v1_out_vertices.end(); ++it)
    {
      v2 = *it;
      curr_value = values_x[graph->pos(v1, v2)];

      if (!double_equals(curr_value, 0.0))
      {
        if (!(visited_nodes[v2]))
        {
          l_nodes.push_back(g.addNode());
          l_nodes_inv.push_back(g_inv.addNode());
          subgraph_to_graph_map[cont] = v2;
          graph_to_subgraph_map[v2] = cont;
          q.push_front(v2);
          visited_nodes[v2] = true;
          visited_nodes_list.push_front(v2);
          cont++;
        }

        LemonGraph::Arc curr_arc = g.addArc(l_nodes[graph_to_subgraph_map[v1]], l_nodes[graph_to_subgraph_map[v2]]);
        LemonGraph::Arc curr_arc_inv = g_inv.addArc(l_nodes_inv[graph_to_subgraph_map[v2]], l_nodes_inv[graph_to_subgraph_map[v1]]);
        capacity[curr_arc] = (int)(round(curr_value * K_PRECISION));

        capacity_inv[curr_arc_inv] = (int)(round(curr_value * K_PRECISION));

        if (v2 != 0)
          nodes_sum[v2] += curr_value;
      }
    }
  } while (!q.empty());

  best_cut = nullptr;
  if ((*CALLBACKS_SELECTION)[K_TYPE_CLIQUE_CONFLICT_CUT])
  {
    // find clique conflict cuts.
    local_best_cut = GenerateCliqueConflictCuts(instance, visited_nodes, nodes_sum, values_x,
                                                visited_nodes_list, subgraph_to_graph_map, graph_to_subgraph_map, g,
                                                l_nodes, capacity, g_inv, l_nodes_inv, capacity_inv, sol, cuts, true, instance.active_conflict_cliques());

    if ((local_best_cut != nullptr) && (local_best_cut->isBetterThan(best_cut)))
      best_cut = local_best_cut;

    if (best_cut != nullptr)
    {
      found_cut = true;

      // adds best cut (most violated).
      IloExpr exp(env);
      for (std::list<std::pair<int, int>>::iterator it = best_cut->lhs_nonzero_coefficients_indexes_.begin(); it != best_cut->lhs_nonzero_coefficients_indexes_.end(); ++it)
      {
        v1 = (*it).first;
        v2 = (*it).second;
        exp += operator*(x[graph->pos(v1, v2)], (best_cut->lhs_)[graph->pos(v1, v2)]);
      }

      for (std::list<int>::iterator it = best_cut->rhs_nonzero_coefficients_indexes_.begin(); it != best_cut->rhs_nonzero_coefficients_indexes_.end(); ++it)
      {
        exp -= y[*it];
      }

      sol.set_cut_added(best_cut->type_, true);
      model.add(exp >= 0);
      // std::cout << "added cut to LP: " << std::endl;
      // for (auto i : best_cut->rhs_nonzero_coefficients_indexes_)
      //   std::cout << i << std::endl;
      // std::cout << std::endl;

      // for (auto [i, j] : best_cut->lhs_nonzero_coefficients_indexes_)
      //   std::cout << i << ", " << j << std::endl;
      // std::cout << std::endl;

      if (root_cuts != nullptr)
      {
        root_cuts->push_back(best_cut);
      }

      exp.end();

      // checks the angle between the each Cut and the most violated one
      for (std::list<UserCutGeneral *>::iterator it = cuts.begin(); it != cuts.end(); ++it)
      {
        curr_cut = static_cast<UserCut *>(*it);

        if (curr_cut == best_cut)
          continue;

        // std::cout << *curr_cut << std::endl;
        //  adds only the ones sufficiently orthogonal to the best_cut
        // std::cout << (*curr_cut)*(*best_cut) << std::endl;
        // getchar(); getchar();
        bool is_sufficiently_orthogonal = double_less((*curr_cut) * (*best_cut), K_CUTS_ANGLE_COSIN_LIMIT);
        if (is_sufficiently_orthogonal)
        {
          IloExpr exp2(env);
          for (std::list<std::pair<int, int>>::iterator it2 = curr_cut->lhs_nonzero_coefficients_indexes_.begin(); it2 != curr_cut->lhs_nonzero_coefficients_indexes_.end(); it2++)
          {
            v1 = (*it2).first;
            v2 = (*it2).second;
            exp2 += operator*(x[graph->pos(v1, v2)], (curr_cut->lhs_)[graph->pos(v1, v2)]);
          }

          for (std::list<int>::iterator it2 = curr_cut->rhs_nonzero_coefficients_indexes_.begin(); it2 != curr_cut->rhs_nonzero_coefficients_indexes_.end(); it2++)
          {
            exp2 -= y[*it2];
          }

          sol.set_cut_added(curr_cut->type_, true);
          model.add(exp2 >= 0);
          // std::cout << "added cut to LP 2: " << std::endl;
          // for (auto i : curr_cut->rhs_nonzero_coefficients_indexes_)
          //   std::cout << i << std::endl;
          // std::cout << std::endl;

          // for (auto [i, j] : curr_cut->lhs_nonzero_coefficients_indexes_)
          //   std::cout << i << ", " << j << std::endl;
          // std::cout << std::endl;

          if (root_cuts != nullptr)
          {
            root_cuts->push_back(curr_cut);
            // sol.set_cut_added(curr_cut->type_,false);
          }

          exp2.end();
        }

        if ((!is_sufficiently_orthogonal) || (root_cuts == nullptr))
        {
          delete curr_cut;
          *it = nullptr;
        }
      }

      if (root_cuts == nullptr)
      {
        delete best_cut;
        best_cut = nullptr;
      }
      cuts.clear();
    }
  }

  return found_cut;
}

int a_var_to_index(int vertex, int budget, int num_vertices)
{
  return budget * num_vertices + vertex;
}

int f_var_to_index(int arc_pos, int budget, int num_arcs)
{
  return budget * num_arcs + arc_pos;
}

std::pair<int, int> index_to_a_var(int index, int num_vertices)
{
  int vertex = index % num_vertices;
  int budget = index / num_vertices;
  return std::pair<int, int>(vertex, budget);
}

static void addArcVertexInferenceCuts(IloCplex &cplex, IloModel &model, IloEnv &env, IloNumVarArray &x, IloNumVarArray &y, Instance &instance, bool solve_relax, Solution<double> &solution)
{
  const Graph *graph = instance.graph();
  GArc *curr_arc = nullptr, *curr_inv_arc = nullptr;
  int num_vertices = graph->num_vertices();
  int cont = 0;

  if (graph->num_arcs() == 0)
    return;

  IloRangeArray cuts(env);
  for (int i = 1; i < num_vertices; ++i)
  {
    for (int j = i + 1; j < num_vertices; ++j)
    {
      curr_arc = (*graph)[i][j];
      curr_inv_arc = (*graph)[j][i];
      if ((curr_arc != nullptr) && (curr_inv_arc != nullptr))
      {
        ++cont;
        IloExpr exp1(env);
        IloExpr exp2(env);
        exp1 = y[i] - y[j] + x[graph->pos(i, j)] + x[graph->pos(j, i)];
        exp2 = y[j] - y[i] + x[graph->pos(i, j)] + x[graph->pos(j, i)];
        cuts.add(exp1 <= 1);
        solution.set_cut_added(K_TYPE_INITIAL_ARC_VERTEX_INFERENCE_CUT, solve_relax);
        if (solve_relax)
          solution.set_cut_found(K_TYPE_INITIAL_ARC_VERTEX_INFERENCE_CUT, solve_relax);
        cuts.add(exp2 <= 1);
        solution.set_cut_added(K_TYPE_INITIAL_ARC_VERTEX_INFERENCE_CUT, solve_relax);
        if (solve_relax)
          solution.set_cut_found(K_TYPE_INITIAL_ARC_VERTEX_INFERENCE_CUT, solve_relax);
        exp1.end();
        exp2.end();
      }
    }
  }

  // cplex.exportModel("Testando.lp");
  // getchar();
  // getchar();
  if (cont > 0)
  {
    // if only solve relaxed problem, then add already as cuts/restrictions
    if (solve_relax)
    {
      model.add(cuts);
    }
    else
    {
      cplex.addUserCuts(cuts);
      cuts.endElements(); // only delete like this if user cuts...if constraints, this command would delete them from model!
    }
    cuts.end();
    // std::cout << "added arc vertex inference cuts!" << std::endl;
  }
  // cplex.exportModel("Testando.lp");
}

void SetPriorityOrder(IloCplex &cplex, IloEnv &env, IloNumVarArray &x, Instance &instance, double *R0)
{
  const Graph *graph = instance.graph();
  int v1 = 0, v2 = 0;

  for (std::list<int>::iterator it = (graph->AdjVerticesOut(v1)).begin(); it != (graph->AdjVerticesOut(v1)).end(); ++it)
  {
    v2 = *it;
    cplex.setPriority(x[graph->pos(v1, v2)], 1.0);
  }

  /*IloNumArray pri_order(env, graph->num_arcs());

  for(v1 = 0; v1 < num_vertices; v1++)
  {
  // highest priorities for arcs leaving source
  if(v1 == 0) curr_pri = 1.0;
  else curr_pri = 0.0;

  for(std::list<int>::iterator it = ((*adj_lists)[v1]).begin(); it != ((*adj_lists)[v1]).end(); it++)
  {
  v2 = *it;
  pri_order[graph->pos(v1,v2)] = curr_pri + 1.0/(R0[v1] + ((*graph)[v1][v2])->dist());
}
}

cplex.setPriorities(x, pri_order);*/
}

void optimize(IloCplex &cplex, IloEnv &env, IloModel &model, std::optional<Formulation> formulation, MasterVariables &master_vars, Instance &instance, double total_time_limit, bool solve_relax, bool apply_benders, bool apply_benders_generic_callback, bool combine_feas_op_cuts, bool separate_benders_cuts_relaxation, bool use_valid_inequalities, bool find_root_cuts, double *R0, double *Rn, std::list<UserCutGeneral *> *initial_cuts, HeuristicSolution *initial_sol, bool export_model, std::list<UserCutGeneral *> *root_cuts, Solution<double> &solution)
{
  // if initial solution already proves the problem to be infeasible, return.
  if (initial_sol != nullptr && initial_sol->is_infeasible_)
  {
    solve_relax ? solution.root_time_ = 0 : solution.milp_time_ += initial_sol->total_time_spent_;
    solution.is_feasible_ = false;
    return;
  }

  const Graph *graph = instance.graph();
  int num_vertices = graph->num_vertices();
  int num_arcs = graph->num_arcs();
  Timestamp *ti = NewTimestamp(), *tf = NewTimestamp();
  Timer *timer = GetTimer();
  timer->Clock(ti);
  bool found_cuts = false;
  // std::optional<IloEnv> worker_env = std::nullopt;
  // std::optional<IloCplex> worker_cplex = std::nullopt;
  // std::optional<IloObjective> worker_obj = std::nullopt;
  // std::optional<DualVariables> worker_vars = std::nullopt;
  WorkerI *worker = nullptr;
  BendersGenericCallbackI *generic_callback = nullptr;
  CallbackArguments arguments{.R0 = R0, .Rn = Rn, .instance = &instance};

  cplex.setParam(IloCplex::Param::WorkMem, 100000);
  std::cout << "limit of memory 100000MB" << std::endl;
  cplex.setParam(IloCplex::IloCplex::Param::MIP::Strategy::File, 3);
  cplex.setOut(env.getNullStream());

  std::vector<bool> *CALLBACKS_SELECTION = GetCallbackSelection();

  if (use_valid_inequalities || find_root_cuts)
  {
    if ((*CALLBACKS_SELECTION)[K_TYPE_CLIQUE_CONFLICT_CUT])
    {
      if (!instance.FoundConflictGraph())
      {
        instance.ComputeConflictGraph();
        solution.num_conflicts_ = instance.num_conflicts();
      }
      if (!(instance.found_maximal_cliques()))
      {
        instance.FindAllMaximalConflictCliquesTomita();
        instance.BuildVerticesToCliquesMapping();
        instance.ResetConflictsCliques();
        solution.num_maximal_cliques_ = instance.conflicts_list().size();

        // int count = 0;
        // for (auto &conflict : instance.conflicts_list())
        // {
        //   ++count;
        //   std::cout << "clique :" << count << " :";
        //   for (auto vertex : conflict)
        //     std::cout << " " << vertex;
        //   std::cout << std::endl;
        // }
      }
    }

    if ((*CALLBACKS_SELECTION)[K_TYPE_INITIAL_ARC_VERTEX_INFERENCE_CUT])
      addArcVertexInferenceCuts(cplex, model, env, master_vars.x, master_vars.y, instance, solve_relax, solution);
  }

  // if solving relaxed problem, force dual cplex method (multithreading showd to be slower! Dual cplex forces single threading and is also useful when
  // the main LP is solved several times, either on Benders or while separating valid inequalities).
  if (solve_relax)
  {
    cplex.setParam(IloCplex::Param::Threads, 1);
    cplex.setParam(IloCplex::Param::RootAlgorithm, IloCplex::Dual);
  }

  if (apply_benders)
  {
    // cplex.setParam(IloCplex::Param::Emphasis::Numerical,1);
    // std::cout << "ENFASE!" << std::endl;
    // create subproblem!
    // The subproblem will be always the same, except for the objective function, which relies
    // on the values of x and y variables. Then, we can create a unique model beforehand.
    if (solve_relax || !apply_benders_generic_callback)
    {
      switch (*formulation)
      {
      case Formulation::baseline:
        worker = new WorkerBaseline(env, instance, master_vars, combine_feas_op_cuts, export_model);
        break;
      case Formulation::single_commodity:
        worker = new WorkerSingleCommodity(env, instance, master_vars, R0, Rn, combine_feas_op_cuts, export_model);
        break;
      default:
        throw "invalid formulation";
      }
    }

    // only add lazy constraint callback if solving the integer problem.
    if (!solve_relax)
    {
      // select if solve Benders using obsolete lazy constraints or the new generic callback.
      if (apply_benders_generic_callback)
      {
        int num_threads = K_MULTI_THREAD ? cplex.getNumCores() : 1;
        // std::cout << "chamou!" << std::endl;
        //  Set up the callback to be used for separating Benders' cuts
        switch (*formulation)
        {
        case Formulation::baseline:
          generic_callback = new BendersGenericCallbackBaseline(env, instance, master_vars, combine_feas_op_cuts, export_model, solution, num_threads);
          break;
        case Formulation::single_commodity:
          generic_callback = new BendersGenericCallbackSingleCommodity(env, instance, master_vars, combine_feas_op_cuts, export_model, solution, R0, Rn, num_threads);
          break;
        default:
          throw "invalid formulation";
        }

        CPXLONG contextmask = IloCplex::Callback::Context::Id::Candidate | IloCplex::Callback::Context::Id::ThreadUp | IloCplex::Callback::Context::Id::ThreadDown;
        if (separate_benders_cuts_relaxation)
          contextmask |= IloCplex::Callback::Context::Id::Relaxation;
        cplex.use(generic_callback, contextmask);
      }
      else
      {
        cplex.setParam(IloCplex::Param::Preprocessing::Presolve, IloFalse);
        // cplex.setParam(IloCplex::Param::Threads, 1);
        // Turn on traditional search for use with control callbacks.
        cplex.setParam(IloCplex::Param::MIP::Strategy::Search,
                       IloCplex::Traditional);

        switch (*formulation)
        {
        case Formulation::baseline:
          cplex.use(BendersLazyCallbackBaseline(env, worker->worker_cplex(), worker->worker_obj(), static_cast<DualVariablesBaseline *>(worker->worker_vars()), master_vars, instance, combine_feas_op_cuts, solution));
          // cplex.use(BendersUserCutCallback(env,*worker_cplex,*worker_obj,*worker_vars,y,x,dual_bound_opt->get(),instance));
          break;
        case Formulation::single_commodity:
          cplex.use(BendersLazyCallbackSingleCommodity(env, worker->worker_cplex(), worker->worker_obj(), static_cast<DualVariablesSingleCommodity *>(worker->worker_vars()), master_vars, arguments, combine_feas_op_cuts, solution));
          break;
        default:
          throw "invalid formulation";
        }
      }
    }
  }

  // if (initial_sol != nullptr)
  // {
  //   IloEnv worker_env = worker->worker_cplex().getEnv();
  //   IloModel worker_model = worker->worker_cplex().getModel();

  //   // Update the objective function in the worker LP.
  //   worker_model.remove(worker->worker_obj());
  //   IloExpr obj_expr = worker->worker_obj().getExpr();
  //   obj_expr.clear();

  //   IloNumArray y_values(env, num_vertices);
  //   for (int i = 0; i < num_vertices; ++i)
  //     y_values[i] = initial_sol->y_values_[i];

  //   IloNumArray x_values(env, num_arcs);
  //   for (int i = 0; i < num_arcs; ++i)
  //     x_values[i] = initial_sol->x_values_[i];

  //   FillObjectiveExpressionDualCompactSingleCommodityContinuousSpace(obj_expr, *static_cast<DualVariablesSingleCommodity *>(worker->worker_vars()), x_values, y_values, std::nullopt, R0, Rn, instance, combine_feas_op_cuts);
  //   worker->worker_obj().setExpr(obj_expr);
  //   worker_model.add(worker->worker_obj());
  //   obj_expr.end();

  //   optimizeLP(worker->worker_cplex(), worker_env, model, instance, -1, R0, Rn, solution);

  //   initial_sol->dual_bound_ = worker->worker_cplex().getObjValue();
  // }

  if (!solve_relax && (initial_sol != nullptr) && (initial_sol->is_feasible_))
  {
    // std::cout << "antes" << std::endl;
    IloNumVarArray start_vars(env, num_vertices + num_arcs);
    IloNumArray start_values(env, num_vertices + num_arcs);

    for (int i = 0; i < num_vertices; ++i)
    {
      start_vars[i] = master_vars.y[i];
      start_values[i] = initial_sol->bitset_vertices_[i];
    }

    for (int i = 0; i < num_arcs; ++i)
    {
      start_vars[i + num_vertices] = master_vars.x[i];
      start_values[i + num_vertices] = initial_sol->bitset_arcs_[i];
    }

    // std::cout << initial_sol->bitset_vertices_ << std::endl;
    // std::cout << initial_sol->bitset_arcs_ << std::endl;
    // std::cout << start_values << std::endl;

    // std::cout << " *** antes de colocar solucao" << std::endl;
    cplex.addMIPStart(start_vars, start_values, IloCplex::MIPStartSolveMIP);
    start_vars.end();
    start_values.end();
    // cplex.setParam(IloCplex::Param::MIP::Tolerances::LowerCutoff, 1.0*(initial_sol->profits_sum_));
    // std::cout << "colocou solucao" << std::endl;

    // cplex.setParam(IloCplex::Param::Advance,2);
    // std::cout << "colocou presolve" << std::endl;
    // std::cout << "depois" << std::endl;

    if (!double_equals(total_time_limit, -1))
      total_time_limit = std::max(0.0, total_time_limit - initial_sol->total_time_spent_);

    solution.milp_time_ += initial_sol->total_time_spent_;
  }

  // if applying lazy constraints, only possible to be single thread! (really? double check!)
  if ((!K_MULTI_THREAD) || (!solve_relax && apply_benders && !apply_benders_generic_callback))
    cplex.setParam(IloCplex::Param::Threads, 1);
  else
  {
    std::cout << "multi thread parallel " << cplex.getNumCores() << std::endl;
    // cplex.setParam(IloCplex::Param::Parallel, IloCplex::Opportunistic);
  }

  if ((initial_cuts != nullptr) && (!(initial_cuts->empty())))
  {
    IloRangeArray root_cuts(env);
    for (std::list<UserCutGeneral *>::iterator it = initial_cuts->begin(); it != initial_cuts->end(); it++)
    {
      switch ((*it)->type_)
      {
      case K_TYPE_CLIQUE_CONFLICT_CUT:
      {
        UserCut *curr_user_cut = static_cast<UserCut *>((*it));
        IloExpr exp(env);
        // std::cout << "cut x: ";
        for (std::list<std::pair<int, int>>::iterator it2 = curr_user_cut->lhs_nonzero_coefficients_indexes_.begin(); it2 != curr_user_cut->lhs_nonzero_coefficients_indexes_.end(); it2++)
        {
          int v1 = (*it2).first;
          int v2 = (*it2).second;
          // std::cout << "(" << v1 << "," << v2 << ") ";
          exp += operator*(master_vars.x[graph->pos(v1, v2)], (curr_user_cut->lhs_)[graph->pos(v1, v2)]);
        }

        // std::cout << "y ";
        for (std::list<int>::iterator it2 = curr_user_cut->rhs_nonzero_coefficients_indexes_.begin(); it2 != curr_user_cut->rhs_nonzero_coefficients_indexes_.end(); it2++)
        {
          // std::cout << *it2 << " ";
          exp -= master_vars.y[*it2];
        }
        // std::cout << std::endl;

        root_cuts.add(exp >= 0);
        solution.set_cut_added((*it)->type_, false);
        // std::cout << "added clique cut" << std::endl;

        exp.end();
        break;
      }
      default:
      {
        std::cout << (*it)->type_ << std::endl;
        throw 5;
      }
      }

      // delete (*it);
      //*it = NULL;
    }

    solve_relax ? model.add(root_cuts) : cplex.addUserCuts(root_cuts);
    root_cuts.endElements();
    root_cuts.end();
  }

  if (solve_relax)
    total_time_limit = -1;

  if (!double_equals(total_time_limit, -1))
  {
    // std::cout << total_time_limit - instance.time_spent_in_preprocessing() << std::endl;
    cplex.setParam(IloCplex::Param::ClockType, 2);
    double time_left = total_time_limit - instance.time_spent_in_preprocessing() - timer->CurrentElapsedTime(ti);
    // std::cout << time_left << std::endl;
    cplex.setParam(IloCplex::Param::TimeLimit, time_left);
  }

  double previous_bound = 0.0, curr_bound = std::numeric_limits<double>::infinity();
  do
  {
    previous_bound = curr_bound;
    // Optimize the problem and obtain solution.
    if (!cplex.solve())
    {
      timer->Clock(tf);
      if (worker != nullptr)
      {
        delete worker;
        worker = nullptr;
      }

      if (generic_callback != nullptr)
      {
        delete generic_callback;
        generic_callback = nullptr;
      }
      // if(apply_benders)
      // {
      //   // delete subproblem objects.
      //   (*worker_cplex).end();
      //   (*worker_obj).end();
      //   (*worker_env).end();
      // }

      std::cout << cplex.getCplexStatus() << std::endl;

      if (solve_relax)
        solution.root_time_ = timer->ElapsedTime(ti, tf);
      else
        solution.milp_time_ += timer->ElapsedTime(ti, tf);
      if ((cplex.getCplexStatus() == IloCplex::Infeasible) || (cplex.getCplexStatus() == IloCplex::InfOrUnbd))
        solution.is_feasible_ = false;
      // SetSolutionStatus(cplex,solution,solve_relax);

      delete (ti);
      ti = nullptr;
      delete (tf);
      tf = nullptr;
      return;
    }
    // getchar(); getchar();
    // cont++;
    // std::cout << cont << std::endl;
    found_cuts = false;

    //  have to add ALL benders cuts, even if stuck at the same bound!
    if (solve_relax) // && double_less(curr_bound,previous_bound,K_TAILING_OFF_TOLERANCE))
    {
      curr_bound = cplex.getObjValue();
      // std::cout << cplex.getCplexStatus() << ": " << curr_bound << std::endl;
      // std::cout << " ****" << previous_bound << " " << curr_bound << std::endl;

      IloNumArray x_values(env);
      IloNumArray y_values(env);
      cplex.getValues(y_values, master_vars.y);
      cplex.getValues(x_values, master_vars.x);

      if (apply_benders)
      {
        IloNum dual_bound_value = cplex.getValue(master_vars.dual_bound);
        IloBool sep_status = IloFalse;

        // std::cout << "current master obj " << cplex.getObjValue() << std::endl;
        // std::cout << "current dual bound " << dual_bound_value << std::endl;

        // Benders' cut separation.
        IloExpr cut_expr(env);
        switch (*formulation)
        {
        case Formulation::baseline:
          sep_status = SeparateBendersCutBaseline(env, master_vars.x, master_vars.y, master_vars.dual_bound, x_values, y_values, dual_bound_value, worker->worker_cplex(), static_cast<DualVariablesBaseline *>(worker->worker_vars()), instance, worker->worker_obj(), cut_expr, combine_feas_op_cuts, solution);
          break;
        case Formulation::single_commodity:
          sep_status = SeparateBendersCutSingleCommodity(env, master_vars.x, master_vars.y, master_vars.dual_bound, x_values, y_values, dual_bound_value, worker->worker_cplex(), static_cast<DualVariablesSingleCommodity *>(worker->worker_vars()), R0, Rn, instance, worker->worker_obj(), cut_expr, combine_feas_op_cuts, solution);
          break;
        default:
          throw "invalid formulation";
        }

        if (sep_status)
        {
          // std::cout << "found new Cut" << std::endl;
          found_cuts = true;
          model.add(cut_expr >= 0);
        }
        cut_expr.end();

        if (export_model)
          cplex.exportModel("model_Benders_compact.lp");
      }

      if ((use_valid_inequalities || find_root_cuts) && double_less(curr_bound, previous_bound, K_TAILING_OFF_TOLERANCE))
        found_cuts |= FindAndAddValidInqualities(cplex, env, model, master_vars.y, master_vars.x, y_values, x_values, instance, solution, root_cuts);

      // Free memory.
      x_values.end();
      y_values.end();
    }
  } while ((found_cuts));

  // if((cplex.getCplexStatus() == IloCplex::Optimal)||(cplex.getCplexStatus() == IloCplex::OptimalTol))
  // {
  //   solution.bitset_vertices_ = std::vector<double>(num_vertices);
  //   solution.bitset_arcs_ = std::vector<double>(graph->num_arcs());

  //   // print solution.
  //   IloNumArray y_solution(env);
  //   cplex.getValues(y_solution,master_vars.y);
  //   for (int i = 0; i < num_vertices; ++i)
  //   {
  //     auto vertex_info = graph->vertices_info()[i];
  //     solution.bitset_vertices_[i] = y_solution[i];
  //     if(initial_sol != nullptr && !double_equals(initial_sol->y_values_[i],y_solution[i]))
  //       std::cout << "y[" << i << "]: " << initial_sol->y_values_[i] << " x " <<  y_solution[i] << " decay " << vertex_info.decay_ratio_ << " profit " << vertex_info.profit_ << std::endl;
  //   }

  //   IloNumArray x_solution(env);
  //   cplex.getValues(x_solution,master_vars.x);
  //   for (int i = 0; i < num_vertices; ++i)
  //   {
  //     for(const int &j: graph->AdjVerticesOut(i))
  //     {
  //       auto arc_pos = graph->pos(i,j);
  //       solution.bitset_arcs_[arc_pos] = x_solution[arc_pos];
  //       // if(initial_sol != nullptr && !double_equals(initial_sol->x_values_[arc_pos],x_solution[arc_pos]))
  //       //   std::cout << "x[" << i << "," << j << "]: " << initial_sol->x_values_[arc_pos] << " x " << x_solution[arc_pos] << std::endl;
  //     }
  //   }

  //   // if(apply_benders)
  //   // {
  //   //   double obj_value = 0.0;
  //   //   auto vertices_info = graph->vertices_info();

  //   //   for(int i = instance.num_mandatory() + 1; i < num_vertices; ++i)
  //   //   {
  //   //     const auto& vertex_info = vertices_info[i];
  //   //     obj_value += vertex_info.profit_ * y_solution[i];
  //   //     // if (formulation == Formulation::single_commodity)
  //   //      obj_value -= vertex_info.decay_ratio_ * instance.limit() * y_solution[i];
  //   //   }

  //   //   std::cout << " solution value without dual: " << obj_value << std::endl;
  //   //   std::cout << "status " << cplex.getCplexStatus() << " " << cplex.getObjValue() - cplex.getValue(master_vars.dual_bound) << " + " << cplex.getValue(master_vars.dual_bound) << " = " << cplex.getObjValue() << std::endl;
  //   // }else
  //   // {
  //   //   double obj_value = 0.0;
  //   //   auto vertices_info = graph->vertices_info();

  //   //   for(int i = instance.num_mandatory() + 1; i < num_vertices; ++i)
  //   //   {
  //   //     const auto& vertex_info = vertices_info[i];
  //   //     obj_value += vertex_info.profit_ * y_solution[i];
  //   //     // if (formulation == Formulation::single_commodity)
  //   //      obj_value -= vertex_info.decay_ratio_ * instance.limit() * y_solution[i];
  //   //   }

  //   //   std::cout << " solution value without dual: " << obj_value << std::endl;
  //   // }

  //   x_solution.end();
  //   y_solution.end();
  // }

  // std::cout << curr_bound << std::endl;

  timer->Clock(tf);
  if (solve_relax)
    solution.root_time_ = timer->ElapsedTime(ti, tf);
  else
    solution.milp_time_ += (timer->ElapsedTime(ti, tf) + instance.time_spent_in_preprocessing());

  // std::cout << "time spent: " << solution.milp_time_ << std::endl;
  SetSolutionStatus(cplex, solution, solve_relax);

  delete (ti);
  ti = nullptr;
  delete (tf);
  tf = nullptr;

  if (worker != nullptr)
  {
    delete worker;
    worker = nullptr;
  }

  if (generic_callback != nullptr)
  {
    delete generic_callback;
    generic_callback = nullptr;
  }

  // // delete subproblem objects.
  // if(apply_benders)
  // {
  //   (*worker_cplex).end();
  //   (*worker_obj).end();
  //   (*worker_env).end();
  // }
}

void optimizeLP(IloCplex &cplex, IloEnv &env, IloModel &model, Instance &instance, double total_time_limit, double *R0, double *Rn, Solution<double> &solution)
{
  const Graph *graph = instance.graph();
  int num_vertices = graph->num_vertices();
  Timestamp *ti = NewTimestamp(), *tf = NewTimestamp();
  Timer *timer = GetTimer();

  cplex.setParam(IloCplex::Param::WorkMem, 15000);
  cplex.setParam(IloCplex::IloCplex::Param::MIP::Strategy::File, 3);
  cplex.setOut(env.getNullStream());

  if (!K_MULTI_THREAD)
    cplex.setParam(IloCplex::Param::Threads, 1);

  if (!double_equals(total_time_limit, -1))
  {
    // std::cout << total_time_limit - instance.time_spent_in_preprocessing_ << std::endl;
    cplex.setParam(IloCplex::Param::ClockType, 2);
    cplex.setParam(IloCplex::Param::TimeLimit, total_time_limit - instance.time_spent_in_preprocessing());
    // std::cout << total_time_limit << " - " << instance.time_spent_in_preprocessing_ << std::endl;
  }

  timer->Clock(ti);
  double curr_bound = std::numeric_limits<double>::infinity();

  if (!cplex.solve())
  {
    timer->Clock(tf);
    timer->ElapsedTime(ti, tf);
    if ((cplex.getCplexStatus() == IloCplex::Infeasible) || (cplex.getCplexStatus() == IloCplex::InfOrUnbd))
      solution.is_feasible_ = false;
    return;
  }

  curr_bound = cplex.getObjValue();

  timer->Clock(tf);
  solution.root_time_ = timer->ElapsedTime(ti, tf);

  SetSolutionStatus(cplex, solution, true);

  delete (ti);
  ti = nullptr;
  delete (tf);
  tf = nullptr;
}

void Benders(Instance &inst, Formulation formulation, double *R0, double *Rn, double time_limit, bool apply_benders_generic_callback, bool combine_feas_op_cuts, bool separate_benders_cuts_relaxation, bool solve_relax, bool use_valid_inequalities, bool find_root_cuts, std::list<UserCutGeneral *> *initial_cuts, HeuristicSolution *initial_sol, bool force_use_all_vehicles, bool export_model, std::list<UserCutGeneral *> *root_cuts, Solution<double> &solution)
{
  const Graph *graph = inst.graph();
  const int num_vertices = graph->num_vertices();
  // const int num_arcs = graph->num_arcs();
  // const int num_routes = inst.num_vehicles();
  const int num_arcs_from_origin = size(graph->AdjVerticesOut(0));

  IloEnv env;
  IloCplex cplex(env);
  IloModel model(env);
  cplex.extract(model);
  MasterVariables master_vars{};

  if (formulation == Formulation::baseline)
    AllocateMasterVariablesBaseline(env, master_vars, inst, force_use_all_vehicles, solve_relax);
  if (formulation == Formulation::single_commodity)
    AllocateMasterVariablesSingleCommodity(env, master_vars, inst, force_use_all_vehicles, solve_relax);
  PopulateByRowCommon(env, model, master_vars, inst, force_use_all_vehicles);
  const int num_mandatory = inst.num_mandatory();
  const auto vertices_info = graph->vertices_info();
  const double route_limit = inst.limit();

  // add objective function.
  IloExpr obj(env);

  for (int i = num_mandatory + 1; i < num_vertices; ++i)
  {
    const auto &vertex_info = vertices_info[i];
    obj += operator*(vertex_info.profit_, master_vars.y[i]);
    if (formulation == Formulation::single_commodity)
      obj -= operator*(vertex_info.decay_ratio_ * route_limit, master_vars.y[i]);
  }

  obj += master_vars.dual_bound;
  model.add(IloMaximize(env, obj));
  obj.end();

  if (export_model)
  {
    // add name to variables.
    for (int i = 0; i < num_vertices; ++i)
    {
      char strnum3[15];
      sprintf(strnum3, "y(%d)", i);
      master_vars.y[i].setName(strnum3);

      for (const auto j : graph->AdjVerticesOut(i))
      {
        char strnum[26];
        sprintf(strnum, "x(%d)(%d)", i, j);
        master_vars.x[graph->pos(i, j)].setName(strnum);
      }
    }

    master_vars.slack.setName("slack");
    master_vars.dual_bound.setName("dual_bound");

    cplex.exportModel("model_Benders_compact.lp");
  }

  optimize(cplex, env, model, formulation, master_vars, inst, time_limit, solve_relax, true, apply_benders_generic_callback, combine_feas_op_cuts, separate_benders_cuts_relaxation, use_valid_inequalities, find_root_cuts, R0, Rn, initial_cuts, initial_sol, export_model, root_cuts, solution);

  // // print solution.
  // IloNumArray y_solution(env);
  // cplex.getValues(y_solution,y);
  // for (int i = 0; i < num_vertices; ++i)
  //   std::cout << "y[" << i << "]: " << y_solution[i] << std::endl;

  // IloNumArray x_solution(env);
  // cplex.getValues(x_solution,x);
  // for (int i = 0; i < num_vertices; ++i)
  //   for(const int &j: graph->AdjVerticesOut(i))
  //     std::cout << "x[" << i << "," << j << "]: " << x_solution[graph->pos(i,j)] << std::endl;

  // std::cout << "slack: " << cplex.getValue(slack) << std::endl;
  // std::cout << "dual bound: " << cplex.getValue(dual_bound) << std::endl;
  // x_solution.end();
  // y_solution.end();

  cplex.end();
  env.end();
}

void CompactBaseline(Instance &inst, double *R0, double *Rn, double time_limit, bool solve_relax, bool use_valid_inequalities, bool find_root_cuts, std::list<UserCutGeneral *> *initial_cuts, HeuristicSolution *initial_sol, bool force_use_all_vehicles, bool export_model, std::list<UserCutGeneral *> *root_cuts, Solution<double> &solution)
{
  const Graph *graph = inst.graph();
  const int num_vertices = graph->num_vertices();
  const int num_arcs = graph->num_arcs();
  const int num_routes = inst.num_vehicles();
  const int budget = inst.uncertainty_budget();

  IloEnv env;
  IloCplex cplex(env);
  // cplex.setParam(IloCplex::Param::Benders::Strategy, IloCplex::BendersStrategyType::BendersFull);
  IloModel model(env);
  cplex.extract(model);

  MasterVariables master_vars{};

  AllocateMasterVariablesBaseline(env, master_vars, inst, force_use_all_vehicles, solve_relax);
  // budget + 1 to consider level 0 of budget! 0,..., budget
  IloNumVarArray a(env, (budget + 1) * num_vertices, 0, IloInfinity, ILOFLOAT);

  PopulateByRowCompactBaseline(cplex, env, model, master_vars, a, inst, R0, Rn, force_use_all_vehicles, export_model);

  optimize(cplex, env, model, std::nullopt, master_vars, inst, time_limit, solve_relax, false, false, false, false, use_valid_inequalities, find_root_cuts, R0, Rn, initial_cuts, initial_sol, export_model, root_cuts, solution);

  // if((cplex.getCplexStatus() == IloCplex::Optimal)||(cplex.getCplexStatus() == IloCplex::OptimalTol))
  // {
  //   // print solution.
  //   IloNumArray y_solution(env);
  //   cplex.getValues(y_solution,y);
  //   for (int i = 0; i < num_vertices; ++i)
  //     std::cout << "y[" << i << "]: " << y_solution[i] << std::endl;

  //   IloNumArray x_solution(env);
  //   cplex.getValues(x_solution,x);
  //   for (int i = 0; i < num_vertices; ++i)
  //     for(const int &j: graph->AdjVerticesOut(i))
  //       std::cout << "x[" << i << "," << j << "]: " << x_solution[graph->pos(i,j)] << std::endl;

  //   IloNumArray a_solution(env);
  //   cplex.getValues(a_solution,a);
  //   for (int budget_iter = 0; budget_iter <= budget; ++budget_iter)
  //     for (int i = 0; i < num_vertices; ++i)
  //       std::cout << "a[" << i << "," << budget_iter << "]: " << a_solution[a_var_to_index(i,budget_iter,num_vertices)] << std::endl;

  //   x_solution.end();
  //   y_solution.end();
  //   a_solution.end();
  // }

  cplex.end();
  env.end();
}

void PrimalSubproblemCompactBaseline(Instance &inst, IloNumArray &x_values, IloNumArray &y_values, double *R0, double *Rn, double time_limit, bool export_model, Solution<double> &solution)
{
  const Graph *graph = inst.graph();
  const int num_vertices = graph->num_vertices();
  const int num_arcs = graph->num_arcs();
  const int budget = inst.uncertainty_budget();

  IloEnv env;
  IloCplex cplex(env);
  IloModel model(env);
  cplex.extract(model);

  // budget + 1 to consider level 0 of budget! 0,..., budget
  IloNumVarArray a(env, (budget + 1) * num_vertices, 0, IloInfinity, ILOFLOAT);

  // empty.
  IloNumVarArray y(env);
  IloNumVarArray x(env);

  PopulateByRowCompactBaselineContinuousSpace(env, model, x, y, a, x_values, y_values, inst);

  // add objective function.
  IloExpr obj(env);

  const int num_mandatory = inst.num_mandatory();
  const auto *vertices_info = inst.graph()->vertices_info();

  for (int i = num_mandatory + 1; i < num_vertices; ++i)
  {
    const auto &vertex_info = vertices_info[i];
    obj -= operator*(vertex_info.decay_ratio_, a[a_var_to_index(i, budget, num_vertices)]);
  }

  model.add(IloMaximize(env, obj));
  obj.end();

  if (export_model)
  {
    // add name to variables.
    for (int i = 0; i < num_vertices; ++i)
    {
      for (int budget_iter = 0; budget_iter <= budget; ++budget_iter)
      {
        char strnum[26];
        sprintf(strnum, "a(%d)(%d)", i, budget_iter);
        a[a_var_to_index(i, budget_iter, num_vertices)].setName(strnum);
      }
    }

    cplex.exportModel("primal_continuous_model_compact_baseline.lp");
  }

  optimizeLP(cplex, env, model, inst, time_limit, R0, Rn, solution);

  // // print solution.
  // IloNumArray a_solution(env);
  // cplex.getValues(a_solution,a);
  // for (int budget_iter = 0; budget_iter <= budget; ++budget_iter)
  //   for (int i = 0; i < num_vertices; ++i)
  //     std::cout << "a[" << i << "," << budget_iter << "]: " << a_solution[a_var_to_index(i,budget_iter,num_vertices)] << std::endl;
  // a_solution.end();
  cplex.end();
  env.end();
}

void PrimalSubproblemCompactSingleCommodity(Instance &inst, IloNumArray &x_values, IloNumArray &y_values, double *R0, double *Rn, double time_limit, bool export_model, Solution<double> &solution)
{
  const Graph *graph = inst.graph();
  const int num_vertices = graph->num_vertices();
  const int num_arcs = graph->num_arcs();
  const int budget = inst.uncertainty_budget();
  const int num_mandatory = inst.num_mandatory();
  const auto *vertices_info = inst.graph()->vertices_info();

  IloEnv env;
  IloCplex cplex(env);
  IloModel model(env);
  cplex.extract(model);

  // budget + 1 to consider level 0 of budget! 0,..., budget
  IloNumVarArray f(env, num_arcs * (budget + 1), 0, IloInfinity, ILOFLOAT);

  // empty.
  IloNumVarArray y(env);
  IloNumVarArray x(env);

  PopulateByRowCompactSingleCommodityContinuousSpace(env, model, x, y, f, x_values, y_values, R0, Rn, inst);

  // add objective function.
  IloExpr obj(env);

  for (int j = num_mandatory + 1; j < num_vertices; ++j)
  {
    const auto &vertex_info = vertices_info[j];
    for (const int &i : graph->AdjVerticesIn(j))
      obj += operator*(vertex_info.decay_ratio_, f[f_var_to_index(graph->pos(i, j), budget, num_arcs)]);
  }

  model.add(IloMaximize(env, obj));
  obj.end();

  if (export_model)
  {
    // add name to variables.
    for (int i = 0; i < num_vertices; ++i)
    {
      for (const auto j : graph->AdjVerticesOut(i))
      {
        for (int budget_iter = 0; budget_iter <= budget; ++budget_iter)
        {
          char strnum[29];
          sprintf(strnum, "f(%d)(%d)(%d)", i, j, budget_iter);
          f[f_var_to_index(graph->pos(i, j), budget_iter, num_arcs)].setName(strnum);
        }
      }
    }
    cplex.exportModel("primal_continuous_model_compact_single_commodity.lp");
  }

  optimizeLP(cplex, env, model, inst, time_limit, R0, Rn, solution);

  // print solution.
  // IloNumArray f_solution(env);
  // cplex.getValues(f_solution,f);
  // for (int budget_iter = 0; budget_iter <= budget; ++budget_iter)
  // {
  //   for(int i = 0; i < num_vertices; ++i)
  //   {
  //     for(const auto j: graph->AdjVerticesOut(i))
  //     {
  //       std::cout << "f[" << i << "," << j << "," << budget_iter << "]: " << f_solution[f_var_to_index(graph->pos(i,j),budget_iter,num_arcs)] << std::endl;
  //     }
  //   }
  // }
  // f_solution.end();

  cplex.end();
  env.end();
}

void AllocateMasterVariablesBaseline(IloEnv &env, MasterVariables &master_vars, const Instance &instance, bool force_use_all_vehicles, bool solve_relax, bool disable_all_binary_vars)
{
  const Graph *graph = instance.graph();
  const int num_vertices = graph->num_vertices();
  const int num_arcs = graph->num_arcs();
  const int num_routes = instance.num_vehicles();

  const double binary_upper_bound = disable_all_binary_vars ? 0.0 : 1.0;

  master_vars.slack = IloNumVar(env, 0, force_use_all_vehicles ? 0 : num_routes, ILOFLOAT);
  master_vars.dual_bound = IloNumVar(env, -num_arcs * instance.limit(), 0, ILOFLOAT);

  if (solve_relax)
  {
    master_vars.x = IloNumVarArray(env, num_arcs, 0.0, binary_upper_bound, ILOFLOAT);
    master_vars.y = IloNumVarArray(env, num_vertices, 0.0, binary_upper_bound, ILOFLOAT);
  }
  else
  {
    master_vars.x = IloNumVarArray(env, num_arcs, 0.0, binary_upper_bound, ILOINT);
    master_vars.y = IloNumVarArray(env, num_vertices, 0.0, binary_upper_bound, ILOINT);
  }
}

void AllocateMasterVariablesSingleCommodity(IloEnv &env, MasterVariables &master_vars, const Instance &instance, bool force_use_all_vehicles, bool solve_relax, bool disable_all_binary_vars)
{
  const Graph *graph = instance.graph();
  const int num_vertices = graph->num_vertices();
  const int num_arcs = graph->num_arcs();
  const int num_routes = instance.num_vehicles();

  const double binary_upper_bound = disable_all_binary_vars ? 0.0 : 1.0;

  master_vars.slack = IloNumVar(env, 0, force_use_all_vehicles ? 0 : num_routes, ILOFLOAT);
  master_vars.dual_bound = IloNumVar(env, 0, num_arcs * instance.limit(), ILOFLOAT);

  if (solve_relax)
  {
    master_vars.x = IloNumVarArray(env, num_arcs, 0.0, binary_upper_bound, ILOFLOAT);
    master_vars.y = IloNumVarArray(env, num_vertices, 0.0, binary_upper_bound, ILOFLOAT);
  }
  else
  {
    master_vars.x = IloNumVarArray(env, num_arcs, 0.0, binary_upper_bound, ILOINT);
    master_vars.y = IloNumVarArray(env, num_vertices, 0.0, binary_upper_bound, ILOINT);
  }
}

void DualSubproblemCompactBaseline(Instance &inst, IloNumArray &x_values, IloNumArray &y_values, double *R0, double *Rn, double time_limit, bool export_model, Solution<double> &solution)
{
  const Graph *graph = inst.graph();
  const int num_vertices = graph->num_vertices();
  const int num_arcs = graph->num_arcs();

  IloEnv env;
  IloCplex cplex(env);
  IloModel model(env);
  cplex.extract(model);

  DualVariablesBaseline *dual_vars = new DualVariablesBaseline(env, inst, false);

  PopulateByRowDualCompactBaselineContinuousSpace(env, model, dual_vars, inst, false);

  // add objective function.
  IloExpr obj_expr(env);
  FillObjectiveExpressionDualCompactBaselineContinuousSpace(obj_expr, *dual_vars, x_values, y_values, std::nullopt, inst, false);
  model.add(IloMinimize(env, obj_expr));
  obj_expr.end();

  if (export_model)
  {
    dual_vars->AddNamesToDualVariables();
    cplex.exportModel("dual_continuous_model_compact_baseline.lp");
  }

  optimizeLP(cplex, env, model, inst, time_limit, R0, Rn, solution);

  // // print solution.
  // IloNumArray a_solution(env);
  // cplex.getValues(a_solution,a);
  // for (int budget_iter = 0; budget_iter <= budget; ++budget_iter)
  //   for (int i = 0; i < num_vertices; ++i)
  //     std::cout << "a[" << i << "," << budget_iter << "]: " << a_solution[a_var_to_index(i,budget_iter,num_vertices)] << std::endl;

  delete dual_vars;
  dual_vars = nullptr;
  cplex.end();
  env.end();
}

void DualSubproblemCompactSingleCommodity(Instance &inst, IloNumArray &x_values, IloNumArray &y_values, double *R0, double *Rn, double time_limit, bool export_model, Solution<double> &solution)
{
  const Graph *graph = inst.graph();
  const int num_vertices = graph->num_vertices();
  const int num_arcs = graph->num_arcs();

  IloEnv env;
  IloCplex cplex(env);
  IloModel model(env);
  cplex.extract(model);

  DualVariablesSingleCommodity *dual_vars = new DualVariablesSingleCommodity(env, inst, false);
  PopulateByRowDualCompactSingleCommodityContinuousSpace(env, model, dual_vars, inst, false);

  // add objective function.
  IloExpr obj_expr(env);
  FillObjectiveExpressionDualCompactSingleCommodityContinuousSpace(obj_expr, *dual_vars, x_values, y_values, std::nullopt, R0, Rn, inst, false);
  model.add(IloMinimize(env, obj_expr));
  obj_expr.end();

  if (export_model)
  {
    dual_vars->AddNamesToDualVariables();
    cplex.exportModel("dual_continuous_model_compact_single_commodity.lp");
  }

  optimizeLP(cplex, env, model, inst, time_limit, R0, Rn, solution);

  // // print solution.
  // IloNumArray a_solution(env);
  // cplex.getValues(a_solution,a);
  // for (int budget_iter = 0; budget_iter <= budget; ++budget_iter)
  //   for (int i = 0; i < num_vertices; ++i)
  //     std::cout << "a[" << i << "," << budget_iter << "]: " << a_solution[a_var_to_index(i,budget_iter,num_vertices)] << std::endl;

  delete dual_vars;
  dual_vars = nullptr;

  cplex.end();
  env.end();
}

void CompactSingleCommodity(Instance &inst, double *R0, double *Rn, double time_limit, bool solve_relax, bool use_valid_inequalities, bool find_root_cuts, std::list<UserCutGeneral *> *initial_cuts, HeuristicSolution *initial_sol, bool force_use_all_vehicles, bool export_model, std::list<UserCutGeneral *> *root_cuts, Solution<double> &solution)
{
  const Graph *graph = inst.graph();
  int num_vertices = graph->num_vertices();
  int num_arcs = graph->num_arcs();
  int num_routes = inst.num_vehicles();
  int budget = inst.uncertainty_budget();

  IloEnv env;
  IloCplex cplex(env);
  IloModel model(env);
  cplex.extract(model);

  MasterVariables master_vars{};

  AllocateMasterVariablesSingleCommodity(env, master_vars, inst, force_use_all_vehicles, solve_relax);
  // budget + 1 to consider level 0 of budget! 0,..., budget
  IloNumVarArray f(env, num_arcs * (budget + 1), 0, IloInfinity, ILOFLOAT);

  PopulateByRowCompactSingleCommodity(cplex, env, model, master_vars, f, inst, R0, Rn, force_use_all_vehicles, export_model);

  // cplex.setParam(IloCplex::Param::Benders::Strategy, IloCplex::BendersStrategyType::BendersFull);
  // IloCplex::LongAnnotation decomp = cplex.newLongAnnotation(IloCplex::BendersAnnotation,CPX_BENDERS_MASTERVALUE);
  // for (IloInt i = 0; i < f.getSize(); ++i)
  //     cplex.setAnnotation(decomp, f[i], CPX_BENDERS_MASTERVALUE+1);

  optimize(cplex, env, model, std::nullopt, master_vars, inst, time_limit, solve_relax, false, false, false, false, use_valid_inequalities, find_root_cuts, R0, Rn, initial_cuts, initial_sol, export_model, root_cuts, solution);

  cplex.end();
  env.end();
}

void PopulateByRowCompactSingleCommodity(IloCplex &cplex, IloEnv &env, IloModel &model, MasterVariables &master_vars, IloNumVarArray &f, const Instance &instance, double *R0, double *Rn, bool force_use_all_vehicles, bool export_model)
{
  const Graph *graph = instance.graph();
  const int num_vehicles = instance.num_vehicles();
  const int num_vertices = graph->num_vertices();
  const int num_arcs = graph->num_arcs();
  const int num_mandatory = instance.num_mandatory();
  GArc *curr_arc = nullptr, *curr_inv_arc = nullptr;
  const auto *vertices_info = graph->vertices_info();
  const double route_limit = instance.limit();

  PopulateByRowCommon(env, model, master_vars, instance, force_use_all_vehicles);

  PopulateByRowCompactSingleCommodityContinuousSpace(env, model, master_vars.x, master_vars.y, f, std::nullopt, std::nullopt, R0, Rn, instance);

  // add objective function.
  IloExpr obj(env);
  const int budget = instance.uncertainty_budget();

  for (int j = num_mandatory + 1; j < num_vertices; ++j)
  {
    const auto &vertex_info = vertices_info[j];
    obj += operator*(vertex_info.profit_, master_vars.y[j]);
    obj -= operator*(vertex_info.decay_ratio_ * route_limit, master_vars.y[j]);

    for (const int &i : graph->AdjVerticesIn(j))
      obj += operator*(vertex_info.decay_ratio_, f[f_var_to_index(graph->pos(i, j), budget, num_arcs)]);
  }

  model.add(IloMaximize(env, obj));
  obj.end();

  if (export_model)
  {
    // add name to variables.
    for (int i = 0; i < num_vertices; ++i)
    {
      char strnum3[15];
      sprintf(strnum3, "y(%d)", i);
      master_vars.y[i].setName(strnum3);

      for (const auto j : graph->AdjVerticesOut(i))
      {
        char strnum[26];
        sprintf(strnum, "x(%d)(%d)", i, j);
        master_vars.x[graph->pos(i, j)].setName(strnum);

        for (int budget_iter = 0; budget_iter <= budget; ++budget_iter)
        {
          char strnum2[29];
          sprintf(strnum2, "f(%d)(%d)(%d)", i, j, budget_iter);
          f[f_var_to_index(graph->pos(i, j), budget_iter, num_arcs)].setName(strnum2);
        }
      }
    }

    master_vars.slack.setName("slack");

    cplex.exportModel("model_compact_single_commodity.lp");
  }
}

static void PopulateByRowCommon(IloEnv &env, IloModel &model, MasterVariables &master_vars, const Instance &instance, bool force_use_all_vehicles)
{
  int num_vehicles = instance.num_vehicles();
  int num_vertices = (instance.graph())->num_vertices();
  int num_mandatory = instance.num_mandatory();
  GArc *curr_arc = nullptr, *curr_inv_arc = nullptr;
  const Graph *graph = instance.graph();

  for (int i = 1; i < num_vertices; ++i)
  {
    IloExpr exp(env);
    auto adj_vertices = graph->AdjVerticesOut(i);
    size_t num_adj_arc = adj_vertices.size();

    for (const auto j : adj_vertices)
      exp += master_vars.x[graph->pos(i, j)];

    if (i <= num_mandatory)
    {
      if (num_adj_arc == 0)
        model.add(master_vars.y[i] == 0); // creates infeasibility if a mandatory node is disconnected from graph!
      else
        model.add(exp == 1);
      model.add(master_vars.y[i] == 1);
    }
    else
    {
      if (num_adj_arc == 0)
        master_vars.y[i].setUB(0);
      else
        model.add(exp == master_vars.y[i]);
    }
    exp.end();
  }

  model.add(master_vars.y[0] == 1);

  IloExpr expi(env);
  for (const auto &j : graph->AdjVerticesOut(0))
    expi += master_vars.x[graph->pos(0, j)];

  expi += master_vars.slack;

  model.add(expi == num_vehicles);
  expi.end();

  if (force_use_all_vehicles)
    master_vars.slack.setUB(0);

  for (int i = 0; i < num_vertices; ++i)
  {
    if (size(graph->AdjVerticesIn(i)) > 0 || size(graph->AdjVerticesOut(i)) > 0)
    {
      IloExpr exp(env);
      for (int j = 0; j < num_vertices; ++j)
      {
        if ((*(instance.graph()))[i][j])
          exp += master_vars.x[graph->pos(i, j)];
        if ((*(instance.graph()))[j][i])
          exp -= master_vars.x[graph->pos(j, i)];
      }

      model.add(exp == 0);
      exp.end();
    }
  }
}

static void PopulateByRowCompactSingleCommodityContinuousSpace(IloEnv &env, IloModel &model, IloNumVarArray &x, IloNumVarArray &y, IloNumVarArray &f, std::optional<std::reference_wrapper<IloNumArray>> x_values, std::optional<std::reference_wrapper<IloNumArray>> y_values, double *R0, double *Rn, const Instance &instance)
{
  const Graph *graph = instance.graph();
  const int num_vertices = graph->num_vertices();
  const int num_arcs = graph->num_arcs();
  const int num_mandatory = instance.num_mandatory();
  GArc *curr_arc = nullptr;
  int arc_pos = 0;
  const auto *vertices_info = graph->vertices_info();
  const int budget = instance.uncertainty_budget();
  const double route_time_limit = instance.limit();
  double vertex_deadline = 0.0;

  for (int budget_iter = 0; budget_iter <= budget; ++budget_iter)
  {
    for (int i = 0; i < num_vertices; ++i)
    {
      const auto &from_vertex_info = vertices_info[i];
      for (const int &j : graph->AdjVerticesOut(i))
      {
        curr_arc = (*graph)[i][j];
        arc_pos = graph->pos(i, j);
        double lb_coef = route_time_limit - R0[i] - from_vertex_info.nominal_service_time_ - curr_arc->distance();
        x_values.has_value() ? model.add(f[f_var_to_index(arc_pos, budget_iter, num_arcs)] <= (lb_coef) * (x_values->get())[arc_pos])
                             : model.add(f[f_var_to_index(arc_pos, budget_iter, num_arcs)] <= operator*(lb_coef, x[arc_pos]));

        // also add upper bound if destination vertex is mandatory or profitable.
        if (j > 0)
        {
          const auto &to_vertex_info = vertices_info[j];
          // if mandatory, the deadline of the vertex D_i is the time limit T.
          vertex_deadline = (j <= num_mandatory) ? route_time_limit : round_decimals(to_vertex_info.profit_ / to_vertex_info.decay_ratio_, 2);

          double ub_coef = max(route_time_limit - vertex_deadline, to_vertex_info.nominal_service_time_ + Rn[j]);
          x_values.has_value() ? model.add(f[f_var_to_index(arc_pos, budget_iter, num_arcs)] >= (ub_coef) * (x_values->get())[arc_pos])
                               : model.add(f[f_var_to_index(arc_pos, budget_iter, num_arcs)] >= operator*(ub_coef, x[arc_pos]));
        }
      }
    }

    // robust constraints.
    for (int i = 1; i < num_vertices; ++i)
    {
      const auto &vertex_info = vertices_info[i];
      // if mandatory, the deadline of the vertex D_i is the time limit T.
      const double vertex_deadline = (i <= num_mandatory) ? route_time_limit : round_decimals(vertex_info.profit_ / vertex_info.decay_ratio_, 2);

      for (const int &j : graph->AdjVerticesOut(i))
      {
        curr_arc = (*graph)[i][j];
        assert(curr_arc);
        arc_pos = graph->pos(i, j);

        IloExpr exp(env);
        exp += f[f_var_to_index(arc_pos, budget_iter, num_arcs)];

        for (const int &k : graph->AdjVerticesIn(i))
          exp -= f[f_var_to_index(graph->pos(k, i), budget_iter, num_arcs)];

        double coef = -(vertex_info.nominal_service_time_ + curr_arc->distance());
        x_values.has_value() ? model.add(exp <= coef * ((x_values->get())[graph->pos(i, j)]))
                             : model.add(exp <= operator*(coef, x[graph->pos(i, j)]));

        exp.end();

        if (budget_iter > 0)
        {
          IloExpr exp2(env);
          exp2 += f[f_var_to_index(arc_pos, budget_iter, num_arcs)];

          for (const int &k : graph->AdjVerticesIn(i))
            exp2 -= f[f_var_to_index(graph->pos(k, i), budget_iter - 1, num_arcs)];

          double coef2 = -(vertex_info.nominal_service_time_ + vertex_info.dev_service_time_ + curr_arc->distance());
          x_values.has_value() ? model.add(exp2 <= coef2 * (x_values->get())[graph->pos(i, j)])
                               : model.add(exp2 <= operator*(coef2, x[graph->pos(i, j)]));

          exp2.end();
        }
      }
    }
  }
}

static void PopulateByRowCompactBaselineContinuousSpace(IloEnv &env, IloModel &model, IloNumVarArray &x, IloNumVarArray &y, IloNumVarArray &a, std::optional<std::reference_wrapper<IloNumArray>> x_values, std::optional<std::reference_wrapper<IloNumArray>> y_values, const Instance &instance)
{
  const int num_vertices = (instance.graph())->num_vertices();
  const int num_mandatory = instance.num_mandatory();
  GArc *curr_arc = nullptr;
  const Graph *graph = instance.graph();
  const auto *vertices_info = graph->vertices_info();
  const int budget = instance.uncertainty_budget();
  const double route_time_limit = instance.limit();

  for (int budget_iter = 0; budget_iter <= budget; ++budget_iter)
  {
    for (const auto &vertex : graph->AdjVerticesOut(0))
    {
      curr_arc = *graph[0][vertex];
      assert(curr_arc);
      y_values.has_value() ? model.add(a[a_var_to_index(vertex, budget_iter, num_vertices)] >= curr_arc->distance() * (y_values->get())[vertex])
                           : model.add(a[a_var_to_index(vertex, budget_iter, num_vertices)] >= operator*(curr_arc->distance(), y[vertex]));
    }

    // i is zero or mandatory.
    for (int i = 0; i <= num_mandatory; ++i)
      model.add(a[a_var_to_index(i, budget_iter, num_vertices)] <= route_time_limit);

    // i is profitable.
    for (int i = num_mandatory + 1; i < num_vertices; ++i)
    {
      const auto &vertex_info = vertices_info[i];
      double vertex_deadline = round_decimals(vertex_info.profit_ / vertex_info.decay_ratio_, 2);
      y_values.has_value() ? model.add(a[a_var_to_index(i, budget_iter, num_vertices)] <= vertex_deadline * (y_values->get())[i])
                           : model.add(a[a_var_to_index(i, budget_iter, num_vertices)] <= operator*(vertex_deadline, y[i]));
    }

    // robust constraints.
    for (int i = 1; i < num_vertices; ++i)
    {
      const auto &vertex_info = vertices_info[i];
      // if mandatory, the deadline of the vertex D_i is the time limit T.
      const double vertex_deadline = (i <= num_mandatory) ? route_time_limit : round_decimals(vertex_info.profit_ / vertex_info.decay_ratio_, 2);

      for (const int &j : graph->AdjVerticesOut(i))
      {
        IloExpr exp(env);
        curr_arc = (*graph)[i][j];
        assert(curr_arc);
        exp += (a[a_var_to_index(i, budget_iter, num_vertices)] + vertex_info.nominal_service_time_ + curr_arc->distance() - a[a_var_to_index(j, budget_iter, num_vertices)]);

        x_values.has_value() ? exp -= (vertex_deadline + vertex_info.nominal_service_time_ + curr_arc->distance()) * (1 - ((x_values->get())[graph->pos(i, j)]))
                             : exp -= operator*(vertex_deadline + vertex_info.nominal_service_time_ + curr_arc->distance(), 1 - x[graph->pos(i, j)]);

        model.add(exp <= 0);
        exp.end();

        if (budget_iter > 0)
        {
          IloExpr exp(env);
          curr_arc = (*graph)[i][j];
          assert(curr_arc);
          exp += (a[a_var_to_index(i, budget_iter - 1, num_vertices)] + vertex_info.nominal_service_time_ + vertex_info.dev_service_time_ + curr_arc->distance() - a[a_var_to_index(j, budget_iter, num_vertices)]);

          x_values.has_value() ? exp -= (vertex_deadline + vertex_info.nominal_service_time_ + vertex_info.dev_service_time_ + curr_arc->distance()) * (1 - (x_values->get())[graph->pos(i, j)])
                               : exp -= operator*(vertex_deadline + vertex_info.nominal_service_time_ + vertex_info.dev_service_time_ + curr_arc->distance(), 1 - x[graph->pos(i, j)]);

          model.add(exp <= 0);
          exp.end();
        }
      }
    }
  }
}

void PopulateByRowCompactBaseline(IloCplex &cplex, IloEnv &env, IloModel &model, MasterVariables &master_vars, IloNumVarArray &a, const Instance &instance, double *R0, double *Rn, bool force_use_all_vehicles, bool export_model)
{
  int num_vertices = (instance.graph())->num_vertices();
  int num_mandatory = instance.num_mandatory();
  GArc *curr_arc = nullptr;
  const Graph *graph = instance.graph();
  const auto *vertices_info = graph->vertices_info();

  PopulateByRowCommon(env, model, master_vars, instance, force_use_all_vehicles);

  PopulateByRowCompactBaselineContinuousSpace(env, model, master_vars.x, master_vars.y, a, std::nullopt, std::nullopt, instance);

  // add objective function.
  IloExpr obj(env);
  const int budget = instance.uncertainty_budget();

  for (int i = num_mandatory + 1; i < num_vertices; ++i)
  {
    const auto &vertex_info = vertices_info[i];
    obj += operator*(vertex_info.profit_, master_vars.y[i]);
    obj -= operator*(vertex_info.decay_ratio_, a[a_var_to_index(i, budget, num_vertices)]);
  }

  model.add(IloMaximize(env, obj));
  obj.end();

  if (export_model)
  {
    // add name to variables.
    for (int i = 0; i < num_vertices; ++i)
    {
      char strnum3[15];
      sprintf(strnum3, "y(%d)", i);
      master_vars.y[i].setName(strnum3);

      for (const auto j : graph->AdjVerticesOut(i))
      {
        char strnum[26];
        sprintf(strnum, "x(%d)(%d)", i, j);
        master_vars.x[graph->pos(i, j)].setName(strnum);
      }

      for (int budget_iter = 0; budget_iter <= budget; ++budget_iter)
      {
        char strnum[26];
        sprintf(strnum, "a(%d)(%d)", i, budget_iter);
        a[a_var_to_index(i, budget_iter, num_vertices)].setName(strnum);
      }
    }

    master_vars.slack.setName("slack");

    cplex.exportModel("model_compact_baseline.lp");
  }
}