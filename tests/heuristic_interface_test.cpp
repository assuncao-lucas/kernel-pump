/**
 * @file heuristic_interface_test.cpp
 * @brief Contract tests for the FeasibilityHeuristic interface as KernelPump drives it
 *        (milestone M5 of docs/spec-heuristic-interface.md).
 *
 * Usage: heuristic_interface_test <scratch dir>
 *
 * The golden-file suite proves the *outcome* of the generalisation is unchanged. This proves the
 * *mechanism*: that KP calls the interface in the order and with the arguments the contract
 * promises, and that each of the kernel-absorption branches still behaves as documented --
 * including the ones no config in settings/all reaches, because all of them set
 * kp.forceBucketVarsIntoKernel=1.
 *
 * Two test-only heuristics are registered with HeuristicsFactory:
 *   - "test-recording": wraps FeasibilityPumpHeuristic, forwards everything, logs every call.
 *                       KP therefore behaves exactly as with kp.heuristic=feaspump.
 *   - "test-scripted":  returns canned outcomes with no solver work, so each absorption branch
 *                       can be driven deliberately.
 *
 * Requires CPLEX (KernelPump needs a real MIPModelI to build the kernel).
 */

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <utils/fileconfig.h>
#include <utils/floats.h>

#include "kernelpump/cpxmodel.h"
#include "kernelpump/feaspump_heuristic.h"
#include "kernelpump/submip_heuristic.h"
#include "kernelpump/kernelpump.h"

using namespace dominiqs;

namespace
{
	int failures = 0;
	int checks = 0;

	void check(bool ok, const std::string &what)
	{
		++checks;
		if (!ok)
			++failures;
		std::cout << (ok ? "  PASS  " : "  FAIL  ") << what << "\n";
	}

	void section(const std::string &title) { std::cout << "\n[" << title << "]\n"; }

	/* ------------------------------------------------------------------ *
	 * call log, shared because the heuristic instances are owned by KP
	 * ------------------------------------------------------------------ */

	struct RunCall
	{
		double timeLimit;
		bool startEmpty;
		double startDist;
		bool startPrimalFeasible;
		bool stopOnNoImprovement;
	};

	struct CallLog
	{
		int constructions = 0;
		int readConfigCalls = 0;
		int initCalls = 0;
		bool readConfigBeforeFirstInit = false;
		bool initAlwaysPrecededRun = true;
		std::vector<const void *> instances;
		std::vector<RunCall> runs;
		std::vector<int> iterationsAfterRun; // to show the counter accumulates across buckets

		void clear() { *this = CallLog(); }
	};

	CallLog gLog;

	/* ------------------------------------------------------------------ *
	 * "test-recording": a transparent spy over the real FP adapter
	 * ------------------------------------------------------------------ */

	class RecordingHeuristic final : public FeasibilityHeuristic
	{
	public:
		RecordingHeuristic()
		{
			++gLog.constructions;
			gLog.instances.push_back(this);
		}

		void readConfig() override
		{
			++gLog.readConfigCalls;
			if (gLog.readConfigCalls == 1)
				gLog.readConfigBeforeFirstInit = (gLog.initCalls == 0);
			fp_.readConfig();
		}

		bool init(MIPModelPtr model, const std::vector<char> &ctype = std::vector<char>()) override
		{
			++gLog.initCalls;
			pendingRun_ = true;
			return fp_.init(model, ctype);
		}

		HeuristicOutcome run(double timeLimit, const HeuristicStart &start,
							 bool stopOnNoImprovement) override
		{
			if (!pendingRun_)
				gLog.initAlwaysPrecededRun = false;
			pendingRun_ = false;
			gLog.runs.push_back({timeLimit, start.empty(), start.dist, start.primalFeasible,
								 stopOnNoImprovement});
			auto outcome = fp_.run(timeLimit, start, stopOnNoImprovement);
			gLog.iterationsAfterRun.push_back(fp_.getIterations());
			return outcome;
		}

		bool foundSolution() const override { return fp_.foundSolution(); }
		void getSolution(std::vector<double> &x) const override { fp_.getSolution(x); }
		double getPrimalBound() const override { return fp_.getPrimalBound(); }
		double getClosestDist() const override { return fp_.getClosestDist(); }
		void getClosestFrac(std::vector<double> &f) const override { fp_.getClosestFrac(f); }
		int getIterations() const override { return fp_.getIterations(); }
		std::string name() const override { return "test-recording"; }

	private:
		FeasibilityPumpHeuristic fp_;
		bool pendingRun_ = false;
	};

	/* ------------------------------------------------------------------ *
	 * "test-scripted": canned outcomes, no solver work
	 * ------------------------------------------------------------------ */

	struct Script
	{
		// one entry per subproblem; the last entry is reused if KP asks for more
		std::vector<bool> initResults{true};
		std::vector<HeuristicOutcome> outcomes{HeuristicOutcome{}};
		std::vector<double> closestDists{INFBOUND};
		// which binaries the reported closest fractional point puts above zero
		std::vector<int> fracSupport;

		void clear() { *this = Script(); }
	};

	Script gScript;

	// Returns BY VALUE on purpose: std::vector<bool>::operator[] yields a proxy, not a bool&, so a
	// `const T&` return here would bind to a dangling temporary.
	template <typename T>
	T at_clamped(const std::vector<T> &v, int i)
	{
		if (v.empty())
			throw std::runtime_error("script vector is empty");
		return v[std::min<int>(i, (int)v.size() - 1)];
	}

	class ScriptedHeuristic final : public FeasibilityHeuristic
	{
	public:
		bool init(MIPModelPtr model, const std::vector<char> & = std::vector<char>()) override
		{
			ncols_ = model->ncols();
			idx_ = step_++;
			return at_clamped(gScript.initResults, idx_);
		}

		HeuristicOutcome run(double, const HeuristicStart &, bool) override
		{
			return at_clamped(gScript.outcomes, idx_);
		}

		bool foundSolution() const override { return at_clamped(gScript.outcomes, idx_).foundSolution; }

		void getSolution(std::vector<double> &x) const override
		{
			x.assign(ncols_, 0.0); // a vector of the right length is all KP requires here
		}

		double getPrimalBound() const override { return 0.0; }
		double getClosestDist() const override { return at_clamped(gScript.closestDists, idx_); }

		void getClosestFrac(std::vector<double> &f) const override
		{
			f.assign(ncols_, 0.0);
			for (int j : gScript.fracSupport)
				if (j >= 0 && j < ncols_)
					f[j] = 0.5;
		}

		// a sentinel, so the test can prove KP forwards this value rather than inventing one
		int getIterations() const override { return 1234; }
		std::string name() const override { return "test-scripted"; }

	private:
		int ncols_ = 0;
		int idx_ = 0;
		int step_ = 0;
	};

	/* ------------------------------------------------------------------ *
	 * fixtures
	 * ------------------------------------------------------------------ */

	std::string writeFile(const std::string &dir, const std::string &name, const std::string &content)
	{
		std::string path = dir + "/" + name;
		{
			std::ofstream out(path);
			out << content;
			/* Fail here, naming the directory, rather than letting the solver report a baffling
			 * "could not open file" later: a missing or unwritable scratch directory is by far the
			 * most likely reason this test cannot run. */
			if (!out)
				throw std::runtime_error("cannot write " + path +
										 " -- does the scratch directory exist and is it writable?");
		}
		return path;
	}

	/* 6 binaries, equality forcing 5 of them to 1.
	 *
	 * With kp.maxBucketSize=2 and flat bucketing, the initial kernel holds 2 binaries and the rest
	 * form 2 buckets of 2. Only 2 (then 4) binaries can be nonzero, so the first two subproblems
	 * are infeasible and the third is not: KP is forced to walk every bucket instead of finishing
	 * on the initial kernel, which is what makes the lifecycle observable. */
	const char *kModelLp =
		"Maximize\n"
		" obj: 6 x1 + 5 x2 + 4 x3 + 3 x4 + 2 x5 + 1 x6\n"
		"Subject To\n"
		" c1: x1 + x2 + x3 + x4 + x5 + x6 = 5\n"
		"Bounds\n"
		"Binaries\n"
		" x1 x2 x3 x4 x5 x6\n"
		"End\n";

	/* Base config. kp.tryEnforceFeasibilityInitialKernel must stay 0: the repair loop would grow
	 * the initial kernel until its LP is feasible, which is exactly what we do not want here.
	 * mipPresolve=0 keeps the toy model intact so the bucket structure stays predictable. */
	std::string baseConfig(const std::string &heuristic, bool forceBucketVars)
	{
		return std::string("solver=cpx\n"
						   "runName=heuristic_interface_test\n"
						   "method=kernelpump\n"
						   "timeLimit=30\n"
						   "mipPresolve=0\n"
						   "multiThreading=0\n"
						   "modelLogging=0\n"
						   "printSol=0\n"
						   "seed=1\n"
						   "kp.tryEnforceFeasibilityInitialKernel=0\n"
						   "kp.buildBucketsByRelaxationLayers=0\n"
						   "kp.buildBucketsConsideringVariableDependency=0\n"
						   "kp.resetFPBasisAtNewPump=0\n"
						   "kp.sortByFractionalPart=0\n"
						   "kp.maxBucketSize=2\n") +
			   "kp.forceBucketVarsIntoKernel=" + (forceBucketVars ? "1" : "0") + "\n" +
			   "kp.heuristic=" + heuristic + "\n" +
			   "fp.frac2int=propround\n"
			   "fp.alpha=0.0\n"
			   "fp.iterLimit=2000\n"
			   "fp.stage1IterLimit=1500\n"
			   "fp.stage2IterLimit=500\n"
			   "fp.stage1NoImprIterLimit=40\n"
			   "fp.stage2NoImprIterLimit=40\n"
			   "fp.imprThreshold=0.9\n"
			   "fp.doStage3=0\n"
			   "fp.numFracsObj=1\n"
			   "fp.firstOptMethod=default\n"
			   "fp.reOptMethod=default\n";
	}

	struct RunResult
	{
		bool initOk = false;
		bool foundSolution = false;
		int numBuckets = 0;
		int lastBucketVisited = 0;
		int firstBucketToIterPump = 0;
		int varsInKernel = 0;
		int iterations = 0;
	};

	// Fresh config + fresh model + fresh KernelPump, so each case is independent.
	RunResult runKp(const std::string &scratch, const std::string &cfgName, const std::string &cfg,
					const std::string &modelPath)
	{
		std::string cfgPath = writeFile(scratch, cfgName, cfg);
		gConfig().load(cfgPath, /*merge=*/false);

		MIPModelPtr model(new CPXModel());
		model->intParam(IntParam::Threads, 1);
		model->readModel(modelPath);

		RunResult r;
		KernelPump kp;
		kp.readConfig();
		r.initOk = kp.Init(model);
		if (r.initOk)
			kp.Run(gConfig().get("timeLimit", 30.0));
		r.foundSolution = kp.foundSolution();
		r.numBuckets = kp.getNumBuckets();
		r.lastBucketVisited = kp.getLastBucketVisited();
		r.firstBucketToIterPump = kp.getFirstBucketToIterPump();
		r.varsInKernel = kp.getNumVarsInKernel();
		r.iterations = kp.getIterations();
		return r;
	}

	/* ------------------------------------------------------------------ *
	 * submip fixtures
	 * ------------------------------------------------------------------ */

	/* Infeasible: two binaries cannot sum to 5. Presolve proves it, so the outcome must be
	 * "proven infeasible", not "limit reached". */
	const char *kInfeasibleLp =
		"Maximize\n"
		" obj: x1 + x2\n"
		"Subject To\n"
		" c1: x1 + x2 = 5\n"
		"Binaries\n"
		" x1 x2\n"
		"End\n";

	/* Feasible with a free choice of how many binaries are set: the MIN-sum objective wants 2,
	 * the MAX-sum objective wants all 6. Used to show the objective override reaches the sub-MIP. */
	const char *kChoiceLp =
		"Minimize\n"
		" obj: 0 x1\n"
		"Subject To\n"
		" c1: x1 + x2 + x3 + x4 + x5 + x6 >= 2\n"
		"Binaries\n"
		" x1 x2 x3 x4 x5 x6\n"
		"End\n";

	struct SubmipRun
	{
		bool foundSolution = false;
		bool progressed = false;
		int ones = 0;
		int nodes = 0;
		bool solutionVerified = false;
	};

	/* Drive SubMipHeuristic directly, bypassing KernelPump, so outcome mapping can be asserted
	 * against models whose answer is known by hand. */
	SubmipRun runSubmip(const std::string &scratch, const std::string &cfgName,
						const std::string &extraKeys, const std::string &modelPath,
						double timeLimit)
	{
		std::string cfg = std::string("solver=cpx\nmethod=submip\nmipPresolve=0\n"
									  "multiThreading=0\nmodelLogging=0\nseed=1\n") +
						  extraKeys;
		std::string cfgPath = writeFile(scratch, cfgName, cfg);
		gConfig().load(cfgPath, /*merge=*/false);

		MIPModelPtr model(new CPXModel());
		model->intParam(IntParam::Threads, 1);
		model->readModel(modelPath);

		SubmipRun r;
		SubMipHeuristic h;
		h.readConfig();
		if (!h.init(model))
			return r;
		auto outcome = h.run(timeLimit, HeuristicStart(), false);
		r.foundSolution = outcome.foundSolution;
		r.progressed = outcome.progressed;
		r.nodes = h.getIterations();
		if (outcome.foundSolution)
		{
			std::vector<double> x;
			h.getSolution(x);
			r.solutionVerified = (x.size() == (size_t)model->ncols()) && model->isSolutionFeasible(x);
			std::vector<char> t(model->ncols());
			model->ctypes(&t[0]);
			for (int j = 0; j < model->ncols(); ++j)
				if (t[j] == 'B' && x[j] > 0.5)
					++r.ones;
		}
		return r;
	}

} // namespace

int main(int argc, char const *argv[])
{
	if (argc < 2)
	{
		std::cout << "usage: heuristic_interface_test <scratch dir>\n";
		return 2;
	}
	const std::string scratch = argv[1];
	std::error_code ec;
	std::filesystem::create_directories(scratch, ec); // harmless if it already exists

	HeuristicsFactory::getInstance().registerClass<RecordingHeuristic>("test-recording");
	HeuristicsFactory::getInstance().registerClass<ScriptedHeuristic>("test-scripted");

	const std::string modelPath = writeFile(scratch, "heur_iface_test.lp", kModelLp);

	try
	{
		/* ============================================================== *
		 * 1. the factory
		 * ============================================================== */
		section("factory");
		{
			std::vector<std::string> ids;
			HeuristicsFactory::getInstance().getIDs(std::back_inserter(ids));
			check(std::find(ids.begin(), ids.end(), "feaspump") != ids.end(),
				  "\"feaspump\" is registered by libkp's static recorder");

			FeasibilityHeuristicPtr h(HeuristicsFactory::getInstance().create("feaspump"));
			check(h != nullptr && h->name() == "feaspump", "create(\"feaspump\") yields a usable object");
			check(HeuristicsFactory::getInstance().create("no-such-heuristic") == nullptr,
				  "an unregistered name yields null (KP turns this into a throw)");
		}

		/* ============================================================== *
		 * 1b. name resolution -- what main.cpp uses to dispatch method=
		 *     and to validate kp.heuristic
		 * ============================================================== */
		section("name resolution (method= and kp.heuristic)");
		{
			check(resolveHeuristicName("feaspump") == "feaspump", "exact name resolves");
			// main.cpp used to upper-case method= before comparing, so every casing must keep working
			check(resolveHeuristicName("FEASPUMP") == "feaspump", "upper case resolves (method=FEASPUMP)");
			check(resolveHeuristicName("FeasPump") == "feaspump", "mixed case resolves (method=FeasPump)");
			check(resolveHeuristicName("nope").empty(), "unknown name resolves to empty");
			check(resolveHeuristicName("").empty(), "empty name resolves to empty");

			// a heuristic registered by this test is dispatchable exactly like feaspump, which is
			// the whole point of M6: a new heuristic needs no change to main.cpp or KernelPump
			check(resolveHeuristicName("test-scripted") == "test-scripted",
				  "a newly registered heuristic resolves too (no main.cpp change needed)");
			check(resolveHeuristicName("TEST-SCRIPTED") == "test-scripted",
				  "...and does so case-insensitively");

			const std::string names = registeredHeuristicNames();
			check(names.find("feaspump") != std::string::npos &&
					  names.find("test-scripted") != std::string::npos,
				  "registeredHeuristicNames() lists them for error messages");
		}

		/* ============================================================== *
		 * 2. interface defaults match FeasibilityPump::pump()'s defaults
		 * ============================================================== */
		section("defaults");
		{
			HeuristicStart s;
			check(s.empty() && s.frac.empty(), "HeuristicStart::frac defaults to empty");
			check(s.dist == INFBOUND, "HeuristicStart::dist defaults to INFBOUND (pump's xStartDist)");
			check(!s.primalFeasible, "HeuristicStart::primalFeasible defaults to false (pump's pFeas)");
			HeuristicOutcome o;
			check(!o.foundSolution && !o.progressed, "HeuristicOutcome defaults to {false,false}");
		}

		/* ============================================================== *
		 * 3. lifecycle, observed through the real FP
		 * ============================================================== */
		section("lifecycle (test-recording over FeasibilityPump)");
		gLog.clear();
		RunResult rec = runKp(scratch, "cfg_recording.cfg", baseConfig("test-recording", true), modelPath);

		check(rec.initOk, "KernelPump::Init succeeded on the fixture model");
		check(rec.numBuckets == 2, "fixture yields 2 buckets (6 binaries, maxBucketSize=2, flat mode)");
		check(gLog.constructions == 1, "exactly ONE heuristic instance constructed per KP run");
		check(gLog.instances.size() == 1 ||
				  std::all_of(gLog.instances.begin(), gLog.instances.end(),
							  [&](const void *p) { return p == gLog.instances.front(); }),
			  "the same instance is reused for every bucket");
		check(gLog.readConfigCalls == 1, "readConfig() called exactly once");
		check(gLog.readConfigBeforeFirstInit, "readConfig() called before the first init()");
		check(gLog.initAlwaysPrecededRun, "every run() was preceded by its own init()");
		check(gLog.initCalls == rec.lastBucketVisited + 1,
			  "init() called once per visited subproblem (lastBucketVisited+1)");
		check(!gLog.runs.empty(), "run() was called at least once");

		// warm start: empty on the first invocation, populated afterwards
		if (!gLog.runs.empty())
		{
			check(gLog.runs.front().startEmpty, "first run() gets an empty warm start");
			check(gLog.runs.front().startDist == INFBOUND, "first run() gets dist == INFBOUND");
			check(!gLog.runs.front().startPrimalFeasible, "first run() gets primalFeasible == false");
		}

		// stopOnNoImprovement is false for the last bucket only
		{
			bool flagsOk = true;
			for (size_t i = 0; i < gLog.runs.size(); ++i)
			{
				// KP computes this as (curr_bucket_index != total_num_buckets - 1); invocation i
				// corresponds to curr_bucket_index == (int)i - 1
				bool expected = ((int)i - 1) != rec.numBuckets - 1;
				if (gLog.runs[i].stopOnNoImprovement != expected)
					flagsOk = false;
			}
			check(flagsOk, "stopOnNoImprovement is false for the last bucket only");
		}

		// all non-last invocations share one budget: min_time_per_bucket, computed once
		if (gLog.runs.size() >= 2)
		{
			bool equalBudgets = true;
			for (size_t i = 1; i + 1 < gLog.runs.size(); ++i)
				if (gLog.runs[i].timeLimit != gLog.runs.front().timeLimit)
					equalBudgets = false;
			check(equalBudgets, "every non-last invocation gets the same time budget");
			check(gLog.runs.back().timeLimit > 0.0, "the last invocation gets a positive budget");
		}

		/* ============================================================== *
		 * 4. absorption branches, driven by canned outcomes
		 * ============================================================== */
		section("absorption branches (test-scripted)");

		// (a) init() == false  ->  whole reference kernel absorbed, loop continues
		{
			gScript.clear();
			gScript.initResults = {false};
			RunResult r = runKp(scratch, "cfg_s_initfalse.cfg", baseConfig("test-scripted", true), modelPath);
			check(!r.foundSolution, "init()==false: no solution reported");
			check(r.lastBucketVisited == r.numBuckets,
				  "init()==false: KP still walks every bucket");
			check(r.firstBucketToIterPump == -1,
				  "init()==false: firstBucketToIterPump stays -1 (never iterated)");
			check(r.varsInKernel == 6, "init()==false: every binary ends up absorbed into the kernel");
		}

		// (b) {false,false}  ->  whole reference kernel absorbed permanently
		{
			gScript.clear();
			gScript.outcomes = {HeuristicOutcome{false, false}};
			RunResult r = runKp(scratch, "cfg_s_notprog.cfg", baseConfig("test-scripted", true), modelPath);
			check(!r.foundSolution, "{false,false}: no solution reported");
			check(r.varsInKernel == 6, "{false,false}: whole reference kernel absorbed each time");
			check(r.firstBucketToIterPump == -1,
				  "{false,false}: firstBucketToIterPump stays -1");
		}

		// (c) {false,true} with forceBucketVarsIntoKernel=1  ->  whole reference kernel absorbed
		{
			gScript.clear();
			gScript.outcomes = {HeuristicOutcome{false, true}};
			gScript.closestDists = {INFBOUND}; // no improvement over KP's running best
			RunResult r = runKp(scratch, "cfg_s_force1.cfg", baseConfig("test-scripted", true), modelPath);
			check(!r.foundSolution, "{false,true} force=1: no solution reported");
			check(r.varsInKernel == 6, "{false,true} force=1: whole reference kernel absorbed");
			check(r.firstBucketToIterPump == 0,
				  "{false,true}: firstBucketToIterPump records the initial kernel as 0");
		}

		// (d) {false,true}, no improvement, force=0  ->  nothing absorbed
		{
			gScript.clear();
			gScript.outcomes = {HeuristicOutcome{false, true}};
			gScript.closestDists = {INFBOUND};
			RunResult r = runKp(scratch, "cfg_s_force0_noimpr.cfg", baseConfig("test-scripted", false), modelPath);
			check(!r.foundSolution, "{false,true} force=0 no improvement: no solution reported");
			check(r.varsInKernel == 2,
				  "{false,true} force=0 no improvement: kernel unchanged (bucket vars leave again)");
		}

		// (e) {false,true}, improving closest point, force=0  ->  only frac>0 binaries absorbed
		{
			gScript.clear();
			gScript.outcomes = {HeuristicOutcome{false, true}};
			gScript.closestDists = {9.0, 8.0, 7.0}; // strictly improving, so the branch is taken
			gScript.fracSupport = {0, 1, 2};		// x1,x2,x3 are "used" by the closest point
			RunResult r = runKp(scratch, "cfg_s_force0_impr.cfg", baseConfig("test-scripted", false), modelPath);
			check(!r.foundSolution, "{false,true} force=0 improving: no solution reported");
			// the kernel may only gain the 3 support columns on top of the 2 it started with,
			// so 5 is the ceiling; wholesale absorption (6) would mean the branch was not taken
			check(r.varsInKernel >= 3 && r.varsInKernel <= 5,
				  "{false,true} force=0 improving: only support binaries absorbed (3..5, never all 6)");
		}

		// (f) {true,true}  ->  solution recorded and the loop breaks early
		{
			gScript.clear();
			gScript.outcomes = {HeuristicOutcome{true, true}};
			gScript.closestDists = {1.0};
			RunResult r = runKp(scratch, "cfg_s_found.cfg", baseConfig("test-scripted", true), modelPath);
			check(r.foundSolution, "{true,true}: solution reported");
			check(r.lastBucketVisited == 0,
				  "{true,true} on the initial kernel: lastBucketVisited == 0 and the loop breaks");
			check(r.varsInKernel == 2, "{true,true}: reference kernel absorbed, nothing more");
		}

		/* ============================================================== *
		 * 5. iteration accounting is NOT reset per bucket
		 * ============================================================== */
		section("shared iteration budget");
		{
			// test-recording forwards to ONE FeasibilityPump, whose iteration counter accumulates
			// across buckets. A per-bucket instance, or a resetTotal() between buckets, would make
			// this sequence drop back towards zero instead of climbing.
			bool nonDecreasing = true;
			for (size_t i = 1; i < gLog.iterationsAfterRun.size(); ++i)
				if (gLog.iterationsAfterRun[i] < gLog.iterationsAfterRun[i - 1])
					nonDecreasing = false;
			check(nonDecreasing,
				  "the iteration counter never decreases between buckets (budget is shared, not reset)");
			check(!gLog.iterationsAfterRun.empty() &&
					  rec.iterations == gLog.iterationsAfterRun.back(),
				  "KernelPump::getIterations() reports the heuristic's latest accumulated total");

			std::cout << "         iterations after each bucket:";
			for (int it : gLog.iterationsAfterRun)
				std::cout << " " << it;
			std::cout << "\n";
			if (!gLog.iterationsAfterRun.empty() && gLog.iterationsAfterRun.back() == 0)
				std::cout << "         NOTE: FP did not need to iterate on this toy fixture, so the two\n"
							 "               checks above are weak here. FP's shared-budget semantics are\n"
							 "               covered by the golden-file suite, which includes a case ending\n"
							 "               at exactly fp.iterLimit.\n";

			// forwarding, proven with a sentinel rather than a real pump
			gScript.clear();
			gScript.outcomes = {HeuristicOutcome{false, true}};
			gScript.closestDists = {INFBOUND};
			RunResult fwd = runKp(scratch, "cfg_s_iters.cfg", baseConfig("test-scripted", true), modelPath);
			check(fwd.iterations == 1234,
				  "KernelPump::getIterations() forwards the heuristic's value verbatim (sentinel 1234)");
		}

		/* ============================================================== *
		 * 5b. the submip heuristic
		 * ============================================================== */
		section("submip: registration and dispatch");
		{
			std::vector<std::string> ids;
			HeuristicsFactory::getInstance().getIDs(std::back_inserter(ids));
			check(std::find(ids.begin(), ids.end(), "submip") != ids.end(),
				  "\"submip\" is registered by libkp's static recorder");
			check(resolveHeuristicName("submip") == "submip", "submip resolves for method=/kp.heuristic");
			check(resolveHeuristicName("SubMIP") == "submip", "...case-insensitively");
			FeasibilityHeuristicPtr h(HeuristicsFactory::getInstance().create("submip"));
			check(h != nullptr && h->name() == "submip", "create(\"submip\") yields a usable object");
			check(h && h->getIterations() == 0, "a fresh submip reports 0 nodes");
		}

		section("submip: outcome mapping");
		{
			const std::string feasPath = modelPath; // 6 binaries, sum = 5: satisfiable
			const std::string infPath = writeFile(scratch, "submip_infeasible.lp", kInfeasibleLp);

			// (a) feasible -> {true, true}, and the point must survive the row check
			auto a = runSubmip(scratch, "cfg_submip_feas.cfg", "submip.nodeLimit=500\n", feasPath, 30.0);
			check(a.foundSolution && a.progressed, "feasible sub-MIP -> {true, true}");
			check(a.solutionVerified, "...and the returned point passes isSolutionFeasible()");

			// (b) proven infeasible -> {false, false}; this is what makes KP absorb the bucket
            auto b = runSubmip(scratch, "cfg_submip_inf.cfg", "submip.nodeLimit=500\n", infPath, 30.0);
			check(!b.foundSolution, "infeasible sub-MIP -> no solution");
			check(!b.progressed, "proven infeasible -> progressed == false (KP absorbs the bucket)");

			// (c) budget exhausted must NEVER be reported as infeasible. This is the case that
			// would differ between backends if status() were used instead of isProvenInfeasible().
			auto c = runSubmip(scratch, "cfg_submip_tiny.cfg", "submip.nodeLimit=1\n", feasPath, 0.0);
			check(!(!c.foundSolution && !c.progressed),
				  "a limit-stopped run is never misreported as infeasible");
			std::cout << (c.foundSolution
							  ? "         NOTE: solved even at a 0s/1-node budget, so the limit path\n"
								"               was not exercised here; (b) still pins the other direction.\n"
							  : "         (limit path exercised: no solution, progressed == true)\n");
		}

		section("submip: the objective override follows FP's initial-LP rules");
		{
			/* Asserted on the objective itself rather than on which solution the solver returns:
			 * with SolutionLimit=1 the first feasible point is whatever CPLEX reaches first, which
			 * on a small model is the same regardless of objective. (The override does change
			 * results on real instances -- e.g. supportcase18 under config16+submip sets 163
			 * binaries with MIN, 170 with the original objective and 174 with MAX.) */
			const std::string choicePath = writeFile(scratch, "submip_choice.lp", kChoiceLp);

			auto objectiveOf = [&](bool null_, bool sum_, bool max_, bool rev_) {
				MIPModelPtr m(new CPXModel());
				m->readModel(choicePath);
				SubMipHeuristic::applyObjectiveOverride(m, null_, sum_, max_, rev_);
				int n = m->ncols();
				std::vector<double> obj(n);
				std::vector<char> t(n);
				m->objcoefs(&obj[0]);
				m->ctypes(&t[0]);
				double binSum = 0.0, nonBinSum = 0.0;
				for (int j = 0; j < n; ++j)
					(t[j] == 'B' ? binSum : nonBinSum) += obj[j];
				return std::make_tuple(binSum, nonBinSum, m->objSense(), m->objOffset());
			};

			int nbin = 0;
			{
				MIPModelPtr m(new CPXModel());
				m->readModel(choicePath);
				std::vector<char> t(m->ncols());
				m->ctypes(&t[0]);
				for (int j = 0; j < m->ncols(); ++j)
					if (t[j] == 'B')
						++nbin;
			}

			auto mn = objectiveOf(false, true, false, false);
			check(std::get<0>(mn) == (double)nbin && std::get<1>(mn) == 0.0,
				  "sum-of-binaries: coefficient 1 on every binary, 0 elsewhere");
			check(std::get<2>(mn) == ObjSense::MIN, "...MaxSense=0 gives MIN (towards all-zeros)");
			check(std::get<3>(mn) == 0.0, "...and the objective offset is cleared");

			auto mx = objectiveOf(false, true, true, false);
			check(std::get<2>(mx) == ObjSense::MAX, "MaxSense=1 gives MAX (towards all-ones)");

			auto nl = objectiveOf(true, false, false, false);
			check(std::get<0>(nl) == 0.0 && std::get<1>(nl) == 0.0,
				  "null objective: every coefficient zeroed");

			auto both = objectiveOf(true, true, true, false);
			check(std::get<0>(both) == 0.0,
				  "null takes precedence over sum-of-binaries, as in FP");

			auto rev = objectiveOf(false, true, false, true);
			check(std::get<2>(rev) == ObjSense::MAX,
				  "reverseObjectiveFunction flips the sense last (MIN -> MAX)");

			auto orig = objectiveOf(false, false, false, false);
			check(std::get<0>(orig) == 0.0 && std::get<2>(orig) == ObjSense::MIN,
				  "no keys set: the model's own objective is left untouched");
		}

		section("submip: no budget shared across invocations");
		{
			// Unlike FP, whose fp.iterLimit is decremented across buckets, each submip invocation
			// gets a fresh node limit. Two runs of the same model must do comparable work.
			auto r1 = runSubmip(scratch, "cfg_submip_b1.cfg", "submip.nodeLimit=500\n", modelPath, 30.0);
			auto r2 = runSubmip(scratch, "cfg_submip_b2.cfg", "submip.nodeLimit=500\n", modelPath, 30.0);
			check(r1.foundSolution == r2.foundSolution,
				  "a second invocation is not starved by the first (no shared budget)");
		}

		section("submip: validateConfig rejects selective absorption");
		{
			std::string cfg = baseConfig("submip", /*forceBucketVars=*/false);
			std::string cfgPath = writeFile(scratch, "cfg_submip_force0.cfg", cfg);
			gConfig().load(cfgPath, /*merge=*/false);
			FeasibilityHeuristicPtr h(HeuristicsFactory::getInstance().create("submip"));
			bool threw = false;
			try
			{
				h->validateConfig(/*drivenByKernelPump=*/true);
			}
			catch (const std::exception &)
			{
				threw = true;
			}
			check(threw, "submip + kp.forceBucketVarsIntoKernel=0 is refused when driven by KP");

			// method=submip is refused outright: standalone it is method=solver capped at
			// submip.nodeLimit, i.e. the same search but artificially limited
			bool threwStandalone = false;
			std::string standaloneMsg;
			try
			{
				h->validateConfig(/*drivenByKernelPump=*/false);
			}
			catch (const std::exception &e)
			{
				threwStandalone = true;
				standaloneMsg = e.what();
			}
			check(threwStandalone, "method=submip is refused: it would duplicate method=solver");
			check(standaloneMsg.find("method=solver") != std::string::npos,
				  "...and the message points at method=solver instead");

			FeasibilityHeuristicPtr fp(HeuristicsFactory::getInstance().create("feaspump"));
			bool fpThrew = false;
			try
			{
				fp->validateConfig(true);
			}
			catch (const std::exception &)
			{
				fpThrew = true;
			}
			check(!fpThrew, "...and feaspump imposes no such constraint");
		}

		section("submip: driven by KernelPump end to end");
		{
			std::string cfg = baseConfig("submip", /*forceBucketVars=*/true);
			RunResult r = runKp(scratch, "cfg_kp_submip.cfg", cfg, modelPath);
			check(r.initOk, "KernelPump::Init succeeds with kp.heuristic=submip");
			check(r.foundSolution, "KP+submip finds a solution on the fixture");
			check(r.lastBucketVisited >= 0 && r.lastBucketVisited <= r.numBuckets,
				  "lastBucketVisited stays in range");
		}

		/* ============================================================== *
		 * 6. Init() failure: the heuristic is never created/used
		 * ============================================================== */
		section("Init() failure");
		{
			// MIP-infeasible model; with mipPresolve=1 CPLEX's presolve proves it and Init bails
			// out before the heuristic is ever asked to do anything.
			const char *infeasibleLp =
				"Maximize\n"
				" obj: x1 + x2\n"
				"Subject To\n"
				" c1: x1 + x2 = 5\n"
				"Binaries\n"
				" x1 x2\n"
				"End\n";
			std::string badPath = writeFile(scratch, "heur_iface_infeasible.lp", infeasibleLp);
			std::string cfg = baseConfig("test-scripted", true);
			// enable presolve so infeasibility is detected inside Init()
			auto pos = cfg.find("mipPresolve=0");
			if (pos != std::string::npos)
				cfg.replace(pos, std::string("mipPresolve=0").size(), "mipPresolve=1");

			gScript.clear();
			RunResult r = runKp(scratch, "cfg_s_infeasible.cfg", cfg, badPath);
			check(!r.initOk, "Init() returns false on a MIP-infeasible model");
			check(r.iterations == 0,
				  "getIterations() is 0 (not a crash) when Init() failed before creating state");
			check(!r.foundSolution, "no solution is reported when Init() failed");
		}
	}
	catch (const std::exception &e)
	{
		std::cout << "\n  FAIL  unexpected exception: " << e.what() << "\n";
		++failures;
		++checks;
	}

	std::cout << "\n" << (checks - failures) << "/" << checks << " checks passed\n";
	if (failures)
		std::cout << "RESULT: " << failures << " FAILURE(S)\n";
	else
		std::cout << "RESULT: all checks passed\n";
	return failures ? 1 : 0;
}
