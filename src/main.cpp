/**
 * @file main.cpp
 * @brief Main App
 *
 * @author Domenico Salvagnin <dominiqs at gmail dot com>
 * @author Gioni Mexi <gionimexi at gmail dot com>
 * 2023
 */

#include <iostream>
#include <string>
#include <algorithm>
#include <iomanip>
#include <fstream>

#include <utils/args_parser.h>
#include <utils/fileconfig.h>
#include <utils/path.h>
#include <utils/floats.h>
#include <utils/timer.h>
#include <utils/randgen.h>
#include <utils/consolelog.h>
#include <utils/path.h>

#include "kernelpump/feaspump.h"
#include "kernelpump/version.h"
#include "kernelpump/kernelpump.h"
#include "kernelpump/solution.h"

#ifdef HAS_CPLEX
#include "kernelpump/cpxmodel.h"
#endif
#ifdef HAS_XPRESS
#include "kernelpump/xprsmodel.h"
#endif
#ifdef HAS_SCIP
#include "kernelpump/scipmodel.h"
#endif
#if defined(HAS_SCIP) && defined(HAS_ORTOOLS)
#include "kernelpump/pdlpmodel.h"
#endif
#include <fmt/format.h>

// macro type savers
#define LOG_ITEM(name, value) consoleLog("{} = {}", name, value)

using namespace dominiqs;

static const uint64_t DEF_SEED = 0;

int main(int argc, char const *argv[])
{
	// config/options
	ArgsParser args;
	args.parse(argc, argv);
	if (args.input.size() < 1)
	{
		consoleError("usage: kp prob_file");
		return -1;
	}
	mergeConfig(args, gConfig());
	std::string solution_folder = gConfig().get("solutionFolder", std::string("../solutions/test/"));
	std::string runName = gConfig().get("runName", std::string("default"));
	std::string testset = gConfig().get("testset", std::string("unknown"));
	std::string solver = gConfig().get("solver", std::string("cpx"));
	bool mipPresolve = gConfig().get("mipPresolve", true);
	bool mipFeasEmphasis = gConfig().get("mipFeasEmphasis", false);
	std::string method = gConfig().get("method", std::string(""));
	bool modelLogging = gConfig().get("modelLogging", true);
	bool multiThreading = gConfig().get("multiThreading", 0);
	bool printSol = gConfig().get("printSol", false);
	double timeLimit = gConfig().get("timeLimit", 1e+20);
	double pdlpTol = gConfig().get("fp.pdlpTol", 1.0e-6);
	double pdlpTolDecreaseFactor = gConfig().get("fp.pdlpTolDecreaseFactor", 1.0);
	double pdlpWarmStart = gConfig().get("fp.pdlpWarmStart", 0);

	std::string probName = getProbName(Path(args.input[0]).getBasename());
	// logger
	consoleInfo("Timestamp: {}", currentDateTime());
	consoleInfo("[config]");
	LOG_ITEM("method", method);
	LOG_ITEM("probName", probName);
	LOG_ITEM("solutionFolder", solution_folder);
	LOG_ITEM("testset", testset);
	LOG_ITEM("solver", solver);
	LOG_ITEM("runName", runName);
	LOG_ITEM("presolve", mipPresolve);
	LOG_ITEM("mipFeasEmphasis", mipFeasEmphasis);
	LOG_ITEM("multiThreading", multiThreading);
	LOG_ITEM("gitHash", KP_GIT_HASH);
	LOG_ITEM("kpVersion", KP_VERSION);
	LOG_ITEM("printSol", printSol);
	LOG_ITEM("timeLimit", timeLimit);
	// seed
	uint64_t seed = gConfig().get<uint64_t>("seed", DEF_SEED);
	LOG_ITEM("seed", seed);
	// seed = generateSeed(seed);
	gConfig().set<uint64_t>("seed", seed);

	bool solveOriginalMIP = false, solveKernelPump = false, solveFeasPump = false;

	std::string upper_method = method;
	std::transform(upper_method.begin(), upper_method.end(), upper_method.begin(), ::toupper);
	if (upper_method.empty())
		throw std::runtime_error(fmt::format("No method selected"));
	else if (upper_method == "SOLVER")
		solveOriginalMIP = true;
	else if (upper_method == "FEASPUMP")
		solveFeasPump = true;
	else if (upper_method == "KERNELPUMP")
		solveKernelPump = true;
	else
		throw std::runtime_error(fmt::format("Selected invalid method {}", method));

	MIPModelPtr model;
	StopWatch watch;
	watch.start();
#ifdef HAS_CPLEX
	if (solver == "cpx")
		model = MIPModelPtr(new CPXModel());
#else
	if (solver == "cpx")
		throw std::runtime_error(fmt::format("Did not compile support for solver {}", solver));
#endif
#ifdef HAS_XPRESS
	if (solver == "xprs")
		model = MIPModelPtr(new XPRSModel());
#else
	if (solver == "xprs")
		throw std::runtime_error(fmt::format("Did not compile support for solver {}", solver));
#endif
#ifdef HAS_SCIP
	if (solver == "scip")
		model = MIPModelPtr(new SCIPModel());
#else
	if (solver == "scip")
		throw std::runtime_error(fmt::format("Did not compile support for solver {}", solver));
#endif
#if defined(HAS_SCIP) && defined(HAS_ORTOOLS)
	if (solver == "pdlp")
		model = MIPModelPtr(new PDLPModel());
#else
	if (solver == "pdlp")
		throw std::runtime_error(fmt::format("Did not compile support for solver {}", solver));
#endif

	if (!model)
		throw std::runtime_error("No solver available");

	DOMINIQS_ASSERT(model);
	auto defaultTimeLimit = model->dblParam(DblParam::TimeLimit);

	double integralityEps = model->dblParam(DblParam::IntegralityTolerance);
	gConfig().set("fp.integralityEps", integralityEps);
#ifndef SILENT_EXEC
	model->logging(modelLogging);
#endif //< SILENT_EXEC
	if (!multiThreading)
		model->intParam(IntParam::Threads, 1); // IntParam::Threads is solver-agnostic!
	try
	{
		model->readModel(args.input[0]);

		// std::vector<double> incumbent{0.3, 0.4, 0.5, 10.0001};
		// auto result = model->computeRelativeIntegralityGap(incumbent);
		// consoleInfo("gap = {}", result);
		// return 0;
		// consoleLog("originalProblem: #rows={} #cols={} #nnz={}",
		// 		   model->nrows(), model->ncols(), model->nnz());

		// this would be to see if there exists a trivial solution (all zero or all one).
		// int num_vars = model->ncols();
		// std::vector<char> xType(num_vars);
		// model->ctypes(&xType[0]);
		// for (int i = 0; i < num_vars; ++i)
		// {
		// 	if (xType[i] == 'B')
		// 	{
		// 		model->ub(i, 0);
		// 	}
		// }
		// model->mipopt();

		std::vector<double> x;
		bool foundSolution = false;
		Solution solution;

		// model->dblParam(DblParam::WorkMem, 100000.0);

		// std::vector<double> frac{0.00000001, 2.0, 0.3, 0.4};
		// std::tie(solution.real_integrality_gap_, solution.num_frac_) = model->computeIntegralityGap(frac);
		// std::cout << solution.real_integrality_gap_ << " " << solution.num_frac_ << std::endl;
		// return 0;

		if (solveOriginalMIP)
		{
			auto defaultNumSolLimit = model->intParam(IntParam::SolutionLimit);
			auto timeLeft = std::max(timeLimit - watch.getElapsed(), 0.0);
			model->dblParam(DblParam::TimeLimit, timeLeft);
			model->intParam(IntParam::SolutionLimit, 1); // stop if/when find first feasible solution!
			if (mipFeasEmphasis)
				model->intParam(IntParam::Emphasis, 1); // emphasis in finding feasible solution.
			model->mipopt();
			foundSolution = !(model->isInfeasibleOrTimeReached());
			x.resize(model->ncols());
			if (foundSolution)
				model->sol(&(x[0]));
			watch.stop();

			// reset configs.
			model->intParam(IntParam::SolutionLimit, defaultNumSolLimit);

			solution.is_feasible_ = foundSolution;
			if (foundSolution)
			{
				solution.projection_integrality_gap_ = 0.0;
				solution.real_integrality_gap_ = 0.0;
			}
			solution.total_time_spent_ = watch.getTotal();
		}
		else if (solveKernelPump)
		{
			// kernel pump
			KernelPump kp;
			kp.readConfig();

			auto result = kp.Init(model);
			if (result)
			{
				auto timeLeft = std::max(timeLimit - watch.getElapsed(), 0.0);
				kp.Run(timeLeft);
				if (kp.foundSolution())
				{
					foundSolution = true;
					kp.getSolution(x);
				}
			}
			watch.stop();
			solution.num_buckets_ = kp.getNumBuckets();				   // does not consider the initial kernel.
			solution.last_bucket_visited_ = kp.getLastBucketVisited(); // 0 == initial kernel.
			solution.is_feasible_ = foundSolution;
			solution.num_iterations_ = kp.getIterations();
			solution.first_bucket_to_iter_pump_ = kp.getFirstBucketToIterPump();
			solution.projection_integrality_gap_ = kp.getClosestDist();
			solution.num_binary_vars_added_ = kp.getNumVarsInKernel();
			solution.num_binary_vars_with_value_one_ = kp.num_binary_vars_with_value_1_in_solution();
			if (foundSolution)
			{
				std::tie(solution.real_integrality_gap_, solution.num_frac_) = model->computeIntegralityGap(x, 0.001);
				// DOMINIQS_ASSERT(equal(solution.real_integrality_gap_, 0.0, integralityEps));
				DOMINIQS_ASSERT(solution.num_frac_ == 0);
				// consoleInfo("gap = {}", solution.real_integrality_gap_);
			}
			else
			{
				std::vector<double> closest_frac;
				kp.getClosestFrac(closest_frac);
				std::tie(solution.real_integrality_gap_, solution.num_frac_) = model->computeIntegralityGap(closest_frac, 0.001);
			}
			std::cout << "gap = " << solution.real_integrality_gap_ << " | num frac= " << solution.num_frac_ << std::endl;
			solution.total_time_spent_ = watch.getTotal();
			solution.time_spent_building_kernel_buckets_ = kp.getTimeSpentBuildingKernelBuckets();

			// std::cout << "terminou tudo com tempo " << solution.total_time_spent_ << std::endl;

			kp.Reset();
		}
		else if (solveFeasPump)
		{
			// set tolerance for PDLP and warm start parameter
			if (solver == "pdlp")
			{
				model->dblParam(DblParam::PdlpTolerance, pdlpTol);
				model->dblParam(DblParam::PdlpToleranceDecreaseFactor, pdlpTolDecreaseFactor);
				model->intParam(IntParam::PdlpWarmStart, pdlpWarmStart);
			}
			// feaspump
			FeasibilityPump fp;
			fp.readConfig();

			auto result = fp.init(model);
			if (result)
			{
				auto timeLeft = std::max(timeLimit - watch.getElapsed(), 0.0);
				fp.pump(timeLeft, false);
				if (fp.foundSolution())
				{
					foundSolution = true;
					fp.getSolution(x);
				}
			}
			watch.stop();

			solution.is_feasible_ = fp.foundSolution();
			solution.num_iterations_ = fp.getIterations();
			solution.projection_integrality_gap_ = fp.getClosestDist();
			if (foundSolution)
			{
				std::tie(solution.real_integrality_gap_, solution.num_frac_) = model->computeIntegralityGap(x, 0.001);
				// DOMINIQS_ASSERT(equal(solution.real_integrality_gap_, 0.0, integralityEps));
				DOMINIQS_ASSERT(solution.num_frac_ == 0);
			}
			else
			{
				std::vector<double> closest_frac;
				fp.getClosestFrac(closest_frac);
				std::tie(solution.real_integrality_gap_, solution.num_frac_) = model->computeIntegralityGap(closest_frac, 0.001);
			}
			std::cout << "gap = " << solution.real_integrality_gap_ << " | num frac= " << solution.num_frac_ << std::endl;
			solution.total_time_spent_ = watch.getTotal();

			fp.resetTotal();
		}

		// check and print solution found.
		if (foundSolution)
		{
			consoleInfo("[Feasible solution found]");
			consoleLog("Total time spent: {:.2f}", watch.getTotal());
			// compute objective in original space
			int n = model->ncols();
			std::vector<double> obj(n);
			model->objcoefs(&obj[0]);
			double objValue = model->objOffset();
			objValue += dotProduct(&obj[0], &x[0], n);

			// check solution for feasibility
			int m = model->nrows();
			std::vector<std::string> rNames(m, "");
			model->rowNames(rNames);
			auto rows = model->rows();
			for (int i = 0; i < m; i++)
			{
				const ConstraintPtr &c = (*rows)[i];
				// model->row(i, c->row, c->sense, c->rhs, c->range);
				if (c->sense == 'N')
					continue;
				if (!c->satisfiedBy(&x[0], 0.001))
					throw std::runtime_error(fmt::format("Constraint {} violated by {}", rNames[i], c->violation(&x[0])));
			}
			consoleLog("Double check feasibility done.");

			std::vector<char> xType(n);
			model->ctypes(&xType[0]);
			double bound = 0;
			// std::vector<double> lbs(n, 0), ubs(n, 0);
			// model->lbs(&(lbs[0]));
			// model->ubs(&(ubs[0]));

			for (int i = 0; i < n; ++i)
			{
				if ((xType[i] == 'B') || (xType[i] == 'I'))
				{
					bound = x[i];
					// if (lessThan(bound, lbs[i]) || greaterThan(bound, ubs[i]))
					// if (!equal(x[i], round(x[i])))
					// 	consoleWarn("{} <= {} <= {}", lbs[i], x[i], ubs[i]);
					// if (!isInteger(x[i], 1.0e-5))
					// 	consoleError("{} <= {} <= {}", lbs[i], x[i], ubs[i]);
					model->lb(i, bound);
					model->ub(i, bound);
				}
			}

			// compute optimal value of solution found (by fixing binary and integers found and re-optimizing over the continuous variables).

			model->dblParam(DblParam::TimeLimit, defaultTimeLimit);
			model->dblParam(DblParam::FeasibilityTolerance, 1.0e-5); // increase a bit the tolerance to avoid infeasibility.
			model->switchToLP();
			model->lpopt('S', false, false);
			double ReOptimizedObjValue = model->objval();

			std::cout << "Solution:" << std::endl;
			std::cout << "=obj= " << std::setprecision(15) << objValue << " | reoptimized= " << ReOptimizedObjValue << std::endl;

			solution.value_ = objValue;
			solution.reopt_value_ = ReOptimizedObjValue;
			if (printSol)
			{
				// print solution
				std::vector<std::string> xNames;
				model->colNames(xNames);
				DOMINIQS_ASSERT(xNames.size() == x.size());
				for (unsigned int i = 0; i < x.size(); i++)
				{
					if (isNotNull(x[i], integralityEps))
						std::cout << xNames[i] << " " << std::setprecision(15) << x[i] << std::endl;
				}
			}
		}

		solution.WriteToFile(solution_folder, runName, probName, seed);
	}
	catch (std::exception &e)
	{
		consoleError(e.what());
	}
	return 0;
}
