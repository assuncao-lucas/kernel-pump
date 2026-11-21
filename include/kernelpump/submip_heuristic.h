/**
 * @file submip_heuristic.h
 * @brief Sub-MIP feasibility heuristic: hand the restricted problem to the backend's MIP solver
 *
 * See docs/spec-submip-heuristic.md.
 *
 * Registered as "submip", so it is available both as `kp.heuristic=submip` (driven per bucket by
 * the Kernel Pump, i.e. kernel growth guided by an exact sub-solver instead of a pump) and as
 * `method=submip` (standalone over the whole problem).
 *
 * PORTABILITY (the reason this class looks the way it does -- see spec section 4):
 *  - mipopt() returns true unconditionally on BOTH CPXModel and GRBModel, and status() returns raw
 *    solver codes that share neither values nor meanings. Outcomes are therefore classified only
 *    with isPrimalFeas(), isProvenInfeasible(), aborted() and the backend-independent
 *    isSolutionFeasible() helper. No solver status code appears in this class.
 *  - Gurobi has no MIP presolve (no crush/uncrush in the C API), so this heuristic never presolves.
 *    It works on a clone of the model it is given, which means every vector it returns is already
 *    in the driver's space and needs no mapping.
 *  - every parameter it sets (SolutionLimit, NodeLimit, TimeLimit, Emphasis, Threads, Seed) and
 *    every attribute it reads (IntAttr::Nodes) is supported by both backends.
 *
 * CROSS-INVOCATION STATE (required by the interface contract):
 * only the node total behind getIterations() accumulates. There is deliberately NO budget shared
 * across invocations -- in particular nothing analogous to FeasibilityPump's fp.iterLimit, which is
 * decremented across buckets. Each invocation gets the driver's time limit and a fresh node limit.
 * A reader who knows FP will assume otherwise.
 */

#ifndef SUBMIP_HEURISTIC_H
#define SUBMIP_HEURISTIC_H

#include <string>
#include <vector>

#include "kernelpump/feasibility_heuristic.h"

namespace dominiqs
{

	class SubMipHeuristic final : public FeasibilityHeuristic
	{
	public:
		void readConfig() override;
		bool init(MIPModelPtr model, const std::vector<char> &ctype = std::vector<char>()) override;
		HeuristicOutcome run(double timeLimit, const HeuristicStart &start,
							 bool stopOnNoImprovement) override;

		bool foundSolution() const override { return found_solution_; }
		void getSolution(std::vector<double> &x) const override;
		double getPrimalBound() const override { return primal_bound_; }

		/** Always the driver's own value: a branch-and-bound call produces no "closest to feasible"
		 *  fractional point of its own. See spec section 7.1 -- this is why submip only makes sense
		 *  with kp.forceBucketVarsIntoKernel=1. */
		double getClosestDist() const override { return closest_dist_; }
		void getClosestFrac(std::vector<double> &frac) const override { frac = closest_frac_; }

		/** Accumulated branch-and-bound nodes, NOT pumping iterations. Reported in the .sol file as
		 *  `# iterations`, so that field means nodes for a submip run. */
		int getIterations() const override { return nodes_; }

		std::string name() const override { return "submip"; }

		void validateConfig(bool drivenByKernelPump) const override;

		/**
		 * Apply the sub-MIP objective override in place: the same three-way switch, with the same
		 * precedence, that FeasibilityPump::solveInitialLP applies to FP's initial LP (null >
		 * sum-of-binaries > original, then the optional sense flip). Static and public so the
		 * precedence rules can be asserted directly rather than inferred from which solution a
		 * solver happens to return first.
		 */
		static void applyObjectiveOverride(MIPModelPtr &m, bool forceNull, bool forceSumVars,
										   bool sumVarsMaxSense, bool reverseSense);

	private:

		// the clone this heuristic works on; never the caller's object
		MIPModelPtr model_;

		// parameters
		int node_limit_ = 500;          // submip.nodeLimit, 0 = unlimited
		bool feas_emphasis_ = true;     // submip.feasEmphasis
		bool multi_threading_ = false;  // multiThreading (top level)
		uint64_t seed_ = 0;             // seed (top level)
		// objective override, read from the same fp.* keys FeasibilityPump::solveInitialLP uses
		bool force_null_obj_ = false;
		bool force_sum_vars_obj_ = false;
		bool force_sum_vars_obj_max_sense_ = false;
		bool reverse_obj_ = false;

		// last-invocation results
		bool found_solution_ = false;
		std::vector<double> solution_;
		double primal_bound_ = INFBOUND;
		double closest_dist_ = INFBOUND;
		std::vector<double> closest_frac_;
		int nodes_ = 0;
	};

} // namespace dominiqs

#endif /* SUBMIP_HEURISTIC_H */
