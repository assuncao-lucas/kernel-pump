/**
 * @file feasibility_heuristic.h
 * @brief Generic feasibility-heuristic interface driven by the Kernel Pump
 *
 * See docs/spec-heuristic-interface.md for the full contract and rationale.
 *
 * The Kernel Pump does not depend on the Feasibility Pump as such: it needs some procedure that,
 * given the sub-MIP induced by the currently active kernel, either returns an integer-feasible
 * point or returns the fractional point it got "closest to feasible", so that KP can decide how to
 * grow the kernel. This header defines that contract; FeasibilityPump is one implementation of it
 * (see feaspump_heuristic.h).
 */

#ifndef FEASIBILITY_HEURISTIC_H
#define FEASIBILITY_HEURISTIC_H

#include <memory>
#include <string>
#include <vector>

#include <utils/factory.h>
#include <utils/floats.h>
#include <utils/singleton.h>

#include "kernelpump/mipmodel.h"

namespace dominiqs
{

	/**
	 * Warm start handed to a heuristic by its driver.
	 */
	struct HeuristicStart
	{
		std::vector<double> frac; /**< starting fractional point; empty => solve the relaxation */
		/**
		 * The driver's own distance-to-integrality for `frac`.
		 *
		 * NOTE: informational only, as far as the current FP implementation is concerned.
		 * FeasibilityPump::pump() logs it but then sets its own closestDist to INFBOUND
		 * (src/feaspump.cpp:691), so a warm start does NOT seed the heuristic's notion of the best
		 * distance seen. KernelPump keeps the running best itself and compares against it. An
		 * implementation may use this value, but no driver may assume that it does.
		 */
		double dist = INFBOUND;
		bool primalFeasible = false; /**< is `frac` primal feasible? */

		bool empty() const { return frac.empty(); }
	};

	/**
	 * What one invocation of a heuristic achieved.
	 */
	struct HeuristicOutcome
	{
		bool foundSolution = false; /**< an integer-feasible point is available via getSolution() */
		/**
		 * The heuristic produced at least one reference point, i.e. it was able to iterate at all.
		 * false => the driver should treat this subproblem as (locally) infeasible; KernelPump
		 * reacts by absorbing the whole reference kernel permanently.
		 *
		 * NOTE: for FeasibilityPump this is `!closestFrac.empty()` (src/feaspump.cpp:822), and a
		 * warm start assigns closestFrac immediately (src/feaspump.cpp:686). So once KP has a warm
		 * start to give, this is effectively always true; the false case shows up mainly on the
		 * first bucket or with kp.resetFPBasisAtNewPump=1.
		 */
		bool progressed = false;
	};

	/**
	 * A procedure that looks for an integer-feasible point of the model it was init'd with, and
	 * that can report a "closest to feasible" fractional point when it fails.
	 *
	 * LIFECYCLE CONTRACT (docs/spec-heuristic-interface.md section 5):
	 *  1. ONE instance per driver run. The driver constructs the heuristic once and reuses it for
	 *     every subproblem. Implementations may therefore carry state across invocations -- and
	 *     FeasibilityPump deliberately does: its iteration limits are decremented by the
	 *     iterations consumed and are NOT restored by init(), making fp.iterLimit a global budget
	 *     shared by all buckets. Constructing one instance per subproblem is a behaviour change,
	 *     not a refactor.
	 *  2. readConfig() exactly once, before the first init(). For FP this is where the RNG is
	 *     seeded and the rounder created; calling it per subproblem would re-seed the RNG.
	 *  3. init() then run(), once per subproblem, in that order. If init() returns false, run() is
	 *     not called.
	 *  4. getClosestDist()/getClosestFrac() describe the LAST invocation only, not the whole run.
	 *     The driver is responsible for keeping the running best across invocations.
	 *  5. An implementation must document its own cross-invocation state in its header.
	 */
	class FeasibilityHeuristic
	{
	public:
		virtual ~FeasibilityHeuristic() = default;

		/** Read own parameters from gConfig(). Called once, before the first init(). */
		virtual void readConfig() {}

		/**
		 * Prepare for a subproblem.
		 * @param model the (sub-)problem to work on; may be modified by the heuristic
		 * @param ctype optional variable-type info, if `model` is an LP
		 * @return false if the subproblem is already infeasible
		 */
		virtual bool init(MIPModelPtr model, const std::vector<char> &ctype = std::vector<char>()) = 0;

		/**
		 * Search for an integer-feasible point.
		 * @param timeLimit wall-clock budget for this invocation
		 * @param start warm start, or a default-constructed HeuristicStart
		 * @param stopOnNoImprovement may the heuristic give up early on stagnation?
		 */
		virtual HeuristicOutcome run(double timeLimit, const HeuristicStart &start,
									 bool stopOnNoImprovement) = 0;

		virtual bool foundSolution() const = 0;
		virtual void getSolution(std::vector<double> &x) const = 0;
		virtual double getPrimalBound() const = 0;

		/** Distance of the closest point to feasibility seen in the LAST invocation. */
		virtual double getClosestDist() const = 0;
		/** Projection of that closest point; empty if there is none. */
		virtual void getClosestFrac(std::vector<double> &frac) const = 0;

		/** Total iterations over all invocations so far (reported as `# iterations`). */
		virtual int getIterations() const = 0;

		/** Registered name, for logging. */
		virtual std::string name() const = 0;
	};

	typedef std::unique_ptr<FeasibilityHeuristic> FeasibilityHeuristicPtr;

	typedef SingletonHolder<Factory<FeasibilityHeuristic, std::string>> HeuristicsFactory;

	/**
	 * Look up a registered heuristic name, ignoring case.
	 *
	 * @return the registered name exactly as it was registered, or "" if there is no such
	 *         heuristic. Matching is case-insensitive so that `method=FEASPUMP` keeps working
	 *         exactly as it did when main.cpp compared upper-cased method strings.
	 *
	 * Register names in lower case; the comparison is case-insensitive either way, but the
	 * lower-case spelling is what appears in config files and in error messages.
	 */
	std::string resolveHeuristicName(const std::string &name);

	/** Comma-separated list of registered heuristic names, for error messages and logging. */
	std::string registeredHeuristicNames();

} // namespace dominiqs

#endif /* FEASIBILITY_HEURISTIC_H */
