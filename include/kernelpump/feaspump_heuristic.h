/**
 * @file feaspump_heuristic.h
 * @brief FeasibilityHeuristic adapter for FeasibilityPump
 *
 * See docs/spec-heuristic-interface.md section 7.1.
 *
 * This is a thin forwarding wrapper, deliberately so: FeasibilityPump itself is left completely
 * untouched by the generalisation, which makes "the FP path behaves exactly as before" provable by
 * inspection (an empty diff on src/feaspump.cpp and include/kernelpump/feaspump.h) rather than only
 * by testing.
 *
 * CROSS-INVOCATION STATE (required by the interface contract, section 5.5):
 * FeasibilityPump is stateful across run() calls and KernelPump depends on it.
 *  - fp.iterLimit / fp.stage1IterLimit / fp.stage2IterLimit are decremented by the iterations
 *    consumed at the end of every pumpLoop (src/feaspump.cpp:1975-1985) and are NOT restored by
 *    init() -> resetPartial(). They are therefore GLOBAL budgets shared by every subproblem.
 *  - the iteration counter behind getIterations() accumulates across invocations; only
 *    resetTotal() clears it, and KernelPump never calls it.
 *  - getClosestDist()/getClosestFrac(), by contrast, are per-invocation: resetPartial() clears
 *    them, and pump() resets closestDist to INFBOUND when given a warm start.
 */

#ifndef FEASPUMP_HEURISTIC_H
#define FEASPUMP_HEURISTIC_H

#include <string>
#include <tuple>
#include <vector>

#include "kernelpump/feaspump.h"
#include "kernelpump/feasibility_heuristic.h"

namespace dominiqs
{

	/** The classical Feasibility Pump, behind the generic interface. Registered as "feaspump". */
	class FeasibilityPumpHeuristic final : public FeasibilityHeuristic
	{
	public:
		void readConfig() override { fp_.readConfig(); }

		bool init(MIPModelPtr model, const std::vector<char> &ctype = std::vector<char>()) override
		{
			return fp_.init(model, ctype);
		}

		HeuristicOutcome run(double timeLimit, const HeuristicStart &start,
							 bool stopOnNoImprovement) override
		{
			auto result = fp_.pump(timeLimit, stopOnNoImprovement, start.frac, start.dist,
								   start.primalFeasible);
			HeuristicOutcome outcome;
			outcome.foundSolution = std::get<0>(result);
			outcome.progressed = std::get<1>(result);
			return outcome;
		}

		bool foundSolution() const override { return fp_.foundSolution(); }
		void getSolution(std::vector<double> &x) const override { fp_.getSolution(x); }
		double getPrimalBound() const override { return fp_.getPrimalBound(); }
		double getClosestDist() const override { return fp_.getClosestDist(); }
		void getClosestFrac(std::vector<double> &frac) const override { fp_.getClosestFrac(frac); }
		int getIterations() const override { return fp_.getIterations(); }
		std::string name() const override { return "feaspump"; }

	private:
		FeasibilityPump fp_;
	};

} // namespace dominiqs

#endif /* FEASPUMP_HEURISTIC_H */
