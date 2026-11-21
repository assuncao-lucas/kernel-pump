/**
 * @file submip_heuristic.cpp
 * @brief Sub-MIP feasibility heuristic + factory registration
 *
 * See docs/spec-submip-heuristic.md.
 */

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <stdexcept>
#include <iostream>
#include <string>
#include <vector>

#include <utils/consolelog.h>
#include <utils/fileconfig.h>
#include <utils/floats.h>

#include "kernelpump/submip_heuristic.h"

using namespace dominiqs;

namespace
{
	/* Default node cap per invocation. 500 follows Mexi, Berthold & Salvagnin (EJCO 2023), who use
	 * that limit for sub-MIP stages inside FP, matching SCIP's conventions for sub-MIP heuristics. */
	const int K_SUBMIP_DEFAULT_NODE_LIMIT = 500;
} // namespace

namespace dominiqs
{

	void SubMipHeuristic::readConfig()
	{
		node_limit_ = gConfig().get("submip.nodeLimit", K_SUBMIP_DEFAULT_NODE_LIMIT);
		feas_emphasis_ = gConfig().get("submip.feasEmphasis", true);
		multi_threading_ = gConfig().get("multiThreading", false);
		seed_ = gConfig().get<uint64_t>("seed", 0);

		/* The sub-MIP objective is NOT a key of its own: the sub-MIP solve is this heuristic's
		 * counterpart of FP's initial LP, so it honours exactly the keys
		 * FeasibilityPump::solveInitialLP honours, with the same precedence. That keeps KP-, KP+
		 * and KP+/- meaning the same thing whichever inner heuristic is driving -- the variant is a
		 * property of the configuration, not of the heuristic. See spec section 6.1. */
		force_null_obj_ = gConfig().get("fp.forceNullObjectiveInitialLP", false);
		force_sum_vars_obj_ = gConfig().get("fp.forceSumVarsObjectiveInitialLP", false);
		force_sum_vars_obj_max_sense_ = gConfig().get("fp.forceSumVarsObjectiveInitialLPMaxSense", false);
		reverse_obj_ = gConfig().get("fp.reverseObjectiveFunction", false);

		consoleInfo("[config submip]");
		consoleLog("submip.nodeLimit = {}", node_limit_);
		consoleLog("submip.feasEmphasis = {}", feas_emphasis_);
		consoleLog("objective = {}", force_null_obj_ ? "null"
								   : force_sum_vars_obj_
										 ? (force_sum_vars_obj_max_sense_ ? "sumBinaries/MAX (towards all-ones)"
																		  : "sumBinaries/MIN (towards all-zeros)")
										 : "original");
		if (reverse_obj_)
			consoleLog("objectiveSense = reversed");

	}

	void SubMipHeuristic::validateConfig(bool drivenByKernelPump) const
	{
		/* Standalone, this heuristic would solve the whole problem at SolutionLimit=1 -- which is
		 * exactly what method=solver already does, except capped at submip.nodeLimit and therefore
		 * strictly worse. Two ways to do one thing, one of them a trap: an instance that
		 * method=solver would solve can be reported as a failure here purely because the node cap
		 * was hit. Refuse and point at the existing method instead. */
		if (!drivenByKernelPump)
			throw std::runtime_error(
				"method=submip is not supported: standalone it would be method=solver capped at "
				"submip.nodeLimit, i.e. the same search but artificially limited. Use method=solver "
				"for a plain first-feasible solve, or kp.heuristic=submip to drive the Kernel Pump "
				"with it");

		/* Under KernelPump: submip reports no closest point of its own, so with selective
		 * absorption the kernel would never grow and every bucket would be tried against the
		 * initial kernel alone. Refuse rather than run a silently pointless configuration
		 * (spec section 7.1). */
		if (!gConfig().get("kp.forceBucketVarsIntoKernel", false))
			throw std::runtime_error(
				"kp.heuristic=submip requires kp.forceBucketVarsIntoKernel=1: submip reports no "
				"closest-to-feasible point of its own, so with selective absorption the kernel "
				"would never grow and every bucket would be retried against the initial kernel");
	}

	bool SubMipHeuristic::init(MIPModelPtr model, const std::vector<char> &)
	{
		DOMINIQS_ASSERT(model);

		/* Work on a clone, never the caller's object: the driver reuses its model across
		 * invocations and must not inherit our parameter changes or solve state. No presolve --
		 * Gurobi has none, and skipping it means every vector we return is already in the driver's
		 * space (spec section 4.2). */
		model_ = model->clone();
		DOMINIQS_ASSERT(model_);

		found_solution_ = false;
		solution_.clear();
		primal_bound_ = INFBOUND;
		return true;
	}

	void SubMipHeuristic::applyObjectiveOverride(MIPModelPtr &m, bool force_null_obj_,
												 bool force_sum_vars_obj_,
												 bool force_sum_vars_obj_max_sense_,
												 bool reverse_obj_)
	{
		const int num_vars = m->ncols();

		/* Mirrors FeasibilityPump::solveInitialLP exactly, INCLUDING the precedence: there it is
		 * `if (forceNull) ... else if (forceSumVars) ...`, so null wins when both are set. Written
		 * as the same if/else-if chain rather than a combined branch, because a combined branch
		 * silently inverted that precedence. */
		if (force_null_obj_)
		{
			std::vector<double> coefs(num_vars, 0.0);
			std::vector<int> indexes(num_vars);
			for (int i = 0; i < num_vars; ++i)
				indexes[i] = i;

			m->objcoefs(num_vars, &(indexes[0]), &(coefs[0]));
			m->objOffset(0.0); //< get rid of offset
		}
		else if (force_sum_vars_obj_)
		{
			std::vector<double> coefs(num_vars, 0.0);
			std::vector<int> indexes(num_vars);
			std::vector<char> xType(num_vars);
			m->ctypes(&(xType[0]));
			for (int i = 0; i < num_vars; ++i)
			{
				if (xType[i] == 'B')
					coefs[i] = 1.0;
				indexes[i] = i;
			}

			m->objcoefs(num_vars, &(indexes[0]), &(coefs[0]));
			m->objSense(force_sum_vars_obj_max_sense_ ? ObjSense::MAX : ObjSense::MIN);
			m->objOffset(0.0); //< get rid of offset
		}

		if (reverse_obj_) // applied last, on top of whatever sense is in force
		{
			const auto curr = m->objSense();
			m->objSense(curr == ObjSense::MAX ? ObjSense::MIN : ObjSense::MAX);
		}
	}

	HeuristicOutcome SubMipHeuristic::run(double timeLimit, const HeuristicStart &start,
										  bool /* stopOnNoImprovement: no stagnation notion here */)
	{
		DOMINIQS_ASSERT(model_);
		HeuristicOutcome outcome;

		found_solution_ = false;
		solution_.clear();
		primal_bound_ = INFBOUND;

		/* Pass the driver's own closest point straight back: a branch-and-bound call produces none
		 * of its own, and returning the unchanged value means the driver's improvement test is
		 * correctly false (spec section 7.1). */
		closest_frac_ = start.frac;
		closest_dist_ = start.dist;

		applyObjectiveOverride(model_, force_null_obj_, force_sum_vars_obj_,
							   force_sum_vars_obj_max_sense_, reverse_obj_);

		if (!multi_threading_)
			model_->intParam(IntParam::Threads, 1);
		model_->seed(static_cast<int>(seed_));
		model_->dblParam(DblParam::TimeLimit, std::max(timeLimit, 0.0));
		model_->intParam(IntParam::SolutionLimit, 1); // stop at the first feasible solution
		if (node_limit_ > 0)
			model_->intParam(IntParam::NodeLimit, node_limit_);
		if (feas_emphasis_)
			model_->intParam(IntParam::Emphasis, 1); // CPLEX feasibility emphasis / Gurobi MIPFocus

		model_->handleCtrlC(true);
		model_->mipopt(); // return value is meaningless on every backend; see the header
		model_->handleCtrlC(false);

		nodes_ += model_->intAttr(IntAttr::Nodes);

		if (model_->aborted())
			return outcome; // {false, false}

		/* Classify using only the portable predicates. isPrimalFeas() says a point is available;
		 * isSolutionFeasible() re-checks it against the rows, the way main.cpp double-checks every
		 * solution before reporting it. A phantom solution handed to KernelPump would corrupt the
		 * whole run, so both must hold. */
		if (model_->isPrimalFeas())
		{
			std::vector<double> x(model_->ncols());
			model_->sol(&(x[0]));
			if (model_->isSolutionFeasible(x))
			{
				solution_ = x;
				found_solution_ = true;
				primal_bound_ = model_->objval();
				outcome.foundSolution = true;
				outcome.progressed = true;
				return outcome;
			}
			consoleWarn("submip: solver reported a primal-feasible point that fails the row check; ignoring");
		}

		/* No usable solution. Proven infeasible means this sub-problem can never succeed, so the
		 * driver should absorb the whole reference kernel; a limit merely means it ran out of
		 * budget and the driver should carry on. isProvenInfeasible() keeps the solver codes behind
		 * the abstraction. */
		outcome.progressed = !model_->isProvenInfeasible();
		return outcome;
	}

	void SubMipHeuristic::getSolution(std::vector<double> &x) const
	{
		DOMINIQS_ASSERT(found_solution_);
		x = solution_; // already in the driver's space: this heuristic never presolves
	}

} // namespace dominiqs

// auto registration

class SUBMIP_FACTORY_RECORDER
{
public:
	SUBMIP_FACTORY_RECORDER()
	{
#ifndef SILENT_EXEC
		std::cout << "Registering SubMipHeuristic...";
#endif //< SILENT_EXEC
		HeuristicsFactory::getInstance().registerClass<SubMipHeuristic>("submip");
#ifndef SILENT_EXEC
		std::cout << "done" << std::endl;
#endif //< SILENT_EXEC
	}
};

SUBMIP_FACTORY_RECORDER my_submip_factory_recorder;
