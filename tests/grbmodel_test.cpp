/**
 * @file grbmodel_test.cpp
 * @brief Interface tests for GRBModel on tiny models (milestone M3 of the Gurobi backend spec).
 *
 * Usage: grbmodel_test <license file> <scratch dir>
 *
 * The license file is passed explicitly to GRBModel (never the user environment's GRB_LICENSE_FILE).
 * The toy models are written into <scratch dir>; all of them fit the size-limited pip license.
 * Expected values are computed by hand in the comments next to each check.
 */

#include "kernelpump/grbmodel.h"
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <utils/floats.h>

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

	bool approx(double a, double b, double tol = 1e-6) { return std::fabs(a - b) <= tol; }

	std::string writeFile(const std::string &dir, const std::string &name, const std::string &content)
	{
		std::string path = dir + "/" + name;
		std::ofstream(path) << content;
		return path;
	}

	/* max x + y  s.t. x + 2y <= 4, 3x + y <= 6, 0 <= x,y <= 10 integer.
	 * LP optimum: x = 1.6, y = 1.2, obj 2.8.  MIP optimum: obj 2 (e.g. x=1,y=1 or x=2,y=0). */
	const char *TINY_MIP =
		"Maximize\n obj: x + y\nSubject To\n c1: x + 2 y <= 4\n c2: 3 x + y <= 6\n"
		"Bounds\n 0 <= x <= 10\n 0 <= y <= 10\nGeneral\n x y\nEnd\n";

	/* x + y >= 3 and x + y <= 1: infeasible LP. */
	const char *INFEASIBLE_LP =
		"Minimize\n obj: x + y\nSubject To\n c1: x + y >= 3\n c2: x + y <= 1\n"
		"Bounds\n 0 <= x <= 10\n 0 <= y <= 10\nEnd\n";

	/* b1 + b2 >= 1 with b1, b2 in [0,1]; KP-style: fixing both upper bounds to 0 makes it infeasible. */
	const char *CONFLICT_LP =
		"Minimize\n obj: b1 + b2 + z\nSubject To\n cover: b1 + b2 >= 1\n link: z - b1 >= 0\n"
		"Bounds\n 0 <= b1 <= 1\n 0 <= b2 <= 1\n 0 <= z <= 5\nEnd\n";

	/* ranged row r1: 2 <= x + y <= 4 (MPS L row with RANGES 2). min x + y -> 2. */
	const char *RANGED_MPS =
		"NAME          RNG\n"
		"ROWS\n N  obj\n L  r1\n"
		"COLUMNS\n    x         obj       1.0          r1        1.0\n"
		"    y         obj       1.0          r1        1.0\n"
		"RHS\n    rhs       r1        4.0\n"
		"RANGES\n    rng       r1        2.0\n"
		"BOUNDS\n UP bnd       x         10.0\n UP bnd       y         10.0\n"
		"ENDATA\n";
} // namespace

int main(int argc, char *argv[])
{
	if (argc < 3)
	{
		std::cerr << "usage: grbmodel_test <license file> <scratch dir>\n";
		return 2;
	}
	const std::string license = argv[1];
	const std::string dir = argv[2];

	try
	{
		std::cout << "[license guard]\n";
		{
			// Must run before any successful construction (the env is created once per process).
			bool threw = false;
			try
			{
				GRBModel guarded("", true);
			}
			catch (const std::exception &)
			{
				threw = true;
			}
			check(threw, "GRBModel(\"\", requireLicenseFile=true) throws before starting Gurobi");
		}

		std::string tinyPath = writeFile(dir, "tiny.lp", TINY_MIP);
		std::string infPath = writeFile(dir, "infeasible.lp", INFEASIBLE_LP);
		std::string conflictPath = writeFile(dir, "conflict.lp", CONFLICT_LP);
		std::string rangedPath = writeFile(dir, "ranged.mps", RANGED_MPS);

		std::cout << "[params set before readModel survive]\n";
		GRBModel m(license, true);
		m.intParam(IntParam::Threads, 1);
		m.dblParam(DblParam::TimeLimit, 30.0);
		m.readModel(tinyPath);
		check(m.intParam(IntParam::Threads) == 1, "Threads == 1 after readModel");
		check(approx(m.dblParam(DblParam::TimeLimit), 30.0), "TimeLimit == 30 after readModel");
		check(approx(m.dblParam(DblParam::FeasibilityTolerance), 1e-6), "FeasibilityTolerance default 1e-6");
		check(approx(m.dblParam(DblParam::IntegralityTolerance), 1e-5), "IntegralityTolerance default 1e-5");
		m.dblParam(DblParam::TimeLimit, INFBOUND);
		check(m.dblParam(DblParam::TimeLimit) == INFBOUND, "TimeLimit INFBOUND round-trips");

		std::cout << "[model data]\n";
		check(m.nrows() == 2 && m.ncols() == 2 && m.nnz() == 4, "2 rows, 2 cols, 4 nonzeros");
		check(m.objSense() == ObjSense::MAX, "objective sense MAX");
		std::vector<char> sen(2), typ(2);
		std::vector<double> rhs(2), rng(2), lb(2), ub(2), obj(2);
		m.sense(sen.data());
		m.rhs(rhs.data());
		m.range(rng.data());
		m.ctypes(typ.data());
		m.lbs(lb.data());
		m.ubs(ub.data());
		m.objcoefs(obj.data());
		check(sen[0] == 'L' && sen[1] == 'L', "senses L, L (translated from '<')");
		check(approx(rhs[0], 4) && approx(rhs[1], 6), "rhs 4, 6");
		check(rng[0] == 0.0 && rng[1] == 0.0, "range 0, 0");
		check(typ[0] == 'I' && typ[1] == 'I', "types I, I");
		check(lb[0] == 0 && lb[1] == 0 && approx(ub[0], 10) && approx(ub[1], 10), "bounds [0,10]");
		check(approx(obj[0], 1) && approx(obj[1], 1), "objective coefficients 1, 1");
		std::vector<std::string> cn, rn;
		m.colNames(cn);
		m.rowNames(rn);
		check(cn.size() == 2 && cn[0] == "x" && cn[1] == "y", "column names x, y");
		check(rn.size() == 2 && rn[0] == "c1" && rn[1] == "c2", "row names c1, c2");

		dominiqs::SparseMatrix rows;
		m.rows(rows);
		bool csrOk = rows.k == 2 && rows.nnz == 4 && rows.matbeg[0] == 0 && rows.matbeg[1] == 2 &&
					 rows.matind[0] == 0 && approx(rows.matval[0], 1) && rows.matind[1] == 1 && approx(rows.matval[1], 2) &&
					 rows.matind[2] == 0 && approx(rows.matval[2], 3) && rows.matind[3] == 1 && approx(rows.matval[3], 1);
		check(csrOk, "rows(): CSR [1 2; 3 1]");
		dominiqs::SparseVector rv;
		char rs;
		double rr, rg;
		m.row(1, rv, rs, rr, rg);
		check(rv.size() == 2 && rs == 'L' && approx(rr, 6) && rg == 0.0, "row(1): 3x + y <= 6");
		dominiqs::SparseVector cv;
		char ct;
		double clb, cub, cobj;
		m.col(1, cv, ct, clb, cub, cobj);
		check(cv.size() == 2 && ct == 'I' && approx(cub, 10) && approx(cobj, 1), "col(1): y in both rows");
		dominiqs::SparseMatrix colsM;
		m.cols(colsM);
		check(colsM.k == 2 && colsM.nnz == 4, "cols(): 2 columns, 4 nonzeros");

		std::cout << "[MIP solve]\n";
		m.intParam(IntParam::SolutionLimit, 2000000000);
		m.mipopt();
		check(!m.isInfeasibleOrTimeReached(), "MIP not infeasible / not time-limited");
		check(m.isPrimalFeas(), "MIP primal feasible");
		check(approx(m.objval(), 2.0), "MIP optimum 2");
		std::vector<double> x(2);
		m.sol(x.data());
		check(approx(x[0] + x[1], 2.0), "MIP solution sums to 2");
		check(approx(m.dblAttr(DblAttr::MIPDualBound), 2.0, 1e-4), "MIP dual bound 2");

		std::cout << "[presolve v1]\n";
		check(m.presolve(), "presolve() returns true");
		check(m.presolvedModel() == nullptr, "presolvedModel() == nullptr (no reductions)");
		check(m.postsolveSolution({1.0, 1.0}) == std::vector<double>({1.0, 1.0}), "postsolveSolution is identity");

		std::cout << "[clone + switchToLP + LP methods]\n";
		auto lp = m.clone();
		check(lp->intParam(IntParam::Threads) == 1, "clone keeps Threads == 1");
		lp->switchToLP();
		lp->ctypes(typ.data());
		check(typ[0] == 'C' && typ[1] == 'C', "switchToLP: all columns continuous");
		m.ctypes(typ.data());
		check(typ[0] == 'I' && typ[1] == 'I', "original model still integer");
		for (char method : {'S', 'P', 'D', 'B', 'A'})
		{
			bool ok = lp->lpopt(method, false, true);
			double v = ok ? lp->objval() : -1;
			check(ok && approx(v, 2.8, method == 'A' ? 1e-4 : 1e-6), std::string("lpopt('") + method + "') -> 2.8");
		}
		lp->lpopt('D', false, true);
		lp->sol(x.data());
		check(approx(x[0], 1.6) && approx(x[1], 1.2), "LP solution (1.6, 1.2)");
		std::vector<double> dj(2);
		lp->reduced_costs(dj.data());
		check(approx(dj[0], 0) && approx(dj[1], 0), "reduced costs of basic vars are 0");
		check(lp->isPrimalFeas(), "LP primal feasible");
		int lpMethodAfter = 0;
		GRBgetintparam(GRBgetenv(lp->getModel()), "Method", &lpMethodAfter);
		check(lpMethodAfter == -1, "lpopt restores Method to its previous value");

		std::cout << "[lazy updates]\n";
		lp->ub(0, 1.0);
		lp->ubs(ub.data());
		check(approx(ub[0], 1.0), "ub(0,1) visible immediately in ubs()");
		lp->lpopt('D', false, false);
		check(approx(lp->objval(), 2.5), "with x <= 1 LP optimum is 2.5");
		lp->fixCol(1, 0.5);
		lp->lbs(lb.data());
		lp->ubs(ub.data());
		check(approx(lb[1], 0.5) && approx(ub[1], 0.5), "fixCol sets both bounds");

		std::cout << "[infinity translation + data modifications]\n";
		lp->addEmptyCol("z", 'C', 0.0, INFBOUND, 0.0);
		check(lp->ncols() == 3, "addEmptyCol adds a column");
		lp->ubs(ub.data(), 2, 2);
		check(ub[0] == INFBOUND, "getter maps GRB_INFINITY back to INFBOUND");
		double rawUb = 0;
		GRBgetdblattrelement(lp->getModel(), "UB", 2, &rawUb);
		check(rawUb >= GRB_INFINITY, "INFBOUND stored as GRB_INFINITY inside Gurobi");
		int idx[2] = {0, 2};
		double val[2] = {1.0, -1.0};
		lp->addRow("d1", idx, val, 2, 'L', 0.5);
		check(lp->nrows() == 3, "addRow adds a row");
		lp->row(2, rv, rs, rr, rg);
		check(rv.size() == 2 && rs == 'L' && approx(rr, 0.5), "row(2) reads back x - z <= 0.5");
		lp->delRow(2);
		lp->delCol(2);
		check(lp->nrows() == 2 && lp->ncols() == 2, "delRow/delCol restore the shape");
		lp->objOffset(7.5);
		check(approx(lp->objOffset(), 7.5), "objOffset round-trips");
		lp->objSense(ObjSense::MIN);
		check(lp->objSense() == ObjSense::MIN, "objSense MIN round-trips");
		boost::dynamic_bitset<> enter(2);
		enter[0] = 1;
		lp->updateModelVarBounds(enter, std::nullopt);
		lp->ubs(ub.data());
		check(approx(ub[0], 1.0), "updateModelVarBounds: entering var gets ub 1");
		bool rangedThrew = false;
		try
		{
			lp->addRow("r", idx, val, 2, 'R', 1.0, 1.0);
		}
		catch (const std::exception &)
		{
			rangedThrew = true;
		}
		check(rangedThrew, "addRow with sense R throws");

		std::cout << "[infeasible LP]\n";
		GRBModel inf(license, true);
		inf.readModel(infPath);
		bool ok = inf.lpopt('D', false, true);
		check(ok, "lpopt returns true (solve ran) on an infeasible LP");
		check(inf.isInfeasibleOrTimeReached(), "isInfeasibleOrTimeReached() true");
		check(!inf.isPrimalFeas(), "isPrimalFeas() false");

		std::cout << "[findSetOfConflictingVariables]\n";
		GRBModel cf(license, true);
		cf.readModel(conflictPath);
		cf.ub(0, 0.0);
		cf.ub(1, 0.0);
		cf.lpopt('D', false, true);
		check(cf.isInfeasibleOrTimeReached(), "with b1,b2 <= 0 the LP is infeasible");
		boost::dynamic_bitset<> inactive(3);
		inactive[0] = 1;
		inactive[1] = 1;
		std::vector<int> conflictRows, conflictVars;
		int colsBefore = cf.ncols(), rowsBefore = cf.nrows();
		cf.findSetOfConflictingVariables(inactive, conflictRows, conflictVars, true, 30.0);
		bool subset = !conflictVars.empty();
		for (int v : conflictVars)
			subset = subset && (v == 0 || v == 1);
		check(subset, "feasrelax: conflicting vars non-empty and within {b1, b2}");
		std::cout << "        (feasrelax picked " << conflictVars.size() << " var(s))\n";
		check(cf.ncols() == colsBefore && cf.nrows() == rowsBefore, "feasrelax leaves the model unchanged");
		cf.ubs(ub.data());
		check(ub[0] == 0.0 && ub[1] == 0.0, "bounds untouched by feasrelax");
		cf.findSetOfConflictingVariables(inactive, conflictRows, conflictVars, false, 30.0);
		bool iisOk = !conflictRows.empty() && conflictRows[0] == 0 && !conflictVars.empty();
		check(iisOk, "IIS path: row 'cover' and a bounded var in the conflict");

		std::cout << "[ranged rows from MPS]\n";
		GRBModel rg2(license, true);
		rg2.readModel(rangedPath);
		std::vector<std::string> rcn;
		rg2.colNames(rcn);
		std::cout << "        ranged MPS read as " << rg2.nrows() << " row(s), " << rg2.ncols() << " col(s):";
		for (auto &n : rcn)
			std::cout << " " << n;
		std::cout << "\n";
		rg2.lpopt('D', false, true);
		check(approx(rg2.objval(), 2.0), "ranged row 2 <= x + y <= 4 respected (min x + y = 2)");
	}
	catch (const std::exception &e)
	{
		std::cout << "  FAIL  unexpected exception: " << e.what() << "\n";
		++failures;
	}

	std::cout << "\n"
			  << (checks - failures) << "/" << checks << " checks passed" << (failures ? ", FAILURES: " + std::to_string(failures) : "") << "\n";
	return failures ? 1 : 0;
}
