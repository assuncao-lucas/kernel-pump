/**
 * @file cpxmodel.cpp
 * @brief Implementation of MIPModelI for CPLEX
 *
 * @author Domenico Salvagnin <dominiqs at gmail dot com>
 * @author Gioni Mexi <gionimexi at gmail dot com>
 * 2023
 */

#include "kernelpump/cpxmodel.h"
#include <signal.h>
#include <cstring>
#include <stdexcept>
#include <utils/consolelog.h>

int CPXModel_UserBreak = 0;

static void userSignalBreak(int signum)
{
	CPXModel_UserBreak = 1;
}

static void throwCplexError(CPXCENVptr env, int status)
{
	const unsigned int BUF_SIZE = 4096;
	char errmsg[BUF_SIZE];
	CPXgeterrorstring(env, status, errmsg);
	int trailer = std::strlen(errmsg) - 1;
	if (trailer >= 0)
		errmsg[trailer] = '\0';
	throw std::runtime_error(errmsg);
}

/* Make a call to a Cplex API function checking its return status */
template <typename Func, typename... Args>
void CPX_CALL(Func cpxfunc, CPXENVptr env, Args &&...args)
{
	int status = cpxfunc(env, std::forward<Args>(args)...);
	if (status)
		throwCplexError(env, status);
}

template <typename Func, typename... Args>
int CPX_CALL_SILENT(Func cpxfunc, CPXENVptr env, Args &&...args) // never throws exception!
{
	return cpxfunc(env, std::forward<Args>(args)...);
}

CPXModel::CPXModel()
{
	int status = 0;

	env = CPXopenCPLEX(&status);
	if (status)
		throwCplexError(nullptr, status);

	lp = CPXcreateprob(env, &status, "");
	if (status)
	{
		CPXcloseCPLEX(&env);
		throwCplexError(nullptr, status);
	}
}

CPXModel::CPXModel(CPXENVptr _env, CPXLPptr _lp, bool _ownEnv, bool _ownLP) : env(_env), lp(_lp), ownEnv(_ownEnv), ownLP(_ownLP)
{
	DOMINIQS_ASSERT(env && lp);
}

CPXModel::~CPXModel()
{
	if (restoreSignalHandler)
		handleCtrlC(false);
	if (ownLP)
		CPXfreeprob(env, &lp);
	if (ownEnv)
		CPXcloseCPLEX(&env);
}

/* Read/Write */
void CPXModel::readModel(const std::string &filename)
{
	DOMINIQS_ASSERT(env && lp);
	CPX_CALL(CPXreadcopyprob, env, lp, filename.c_str(), nullptr);
}

void CPXModel::writeModel(const std::string &filename, const std::string &format) const
{
	DOMINIQS_ASSERT(env && lp);
	CPX_CALL(CPXwriteprob, env, lp, filename.c_str(), nullptr);
}

void CPXModel::writeSol(const std::string &filename) const
{
	DOMINIQS_ASSERT(env && lp);
	CPX_CALL(CPXsolwrite, env, lp, filename.c_str());
}

int CPXModel::status() const
{
	DOMINIQS_ASSERT(env && lp);
	return CPXgetstat(env, lp);
}

/* Solve */
bool CPXModel::lpopt(char method, bool decrease_tol, bool initial)
{
	// every call to solve must be SILENT to avoid throwing error when sub-problems are infeasible!
	DOMINIQS_ASSERT(env && lp);
	int status = 0;

	switch (method)
	{
	case 'S':
		status = CPX_CALL_SILENT(CPXlpopt, env, lp);
		break;
	case 'P':
		status = CPX_CALL_SILENT(CPXprimopt, env, lp);
		break;
	case 'D':
		status = CPX_CALL_SILENT(CPXdualopt, env, lp);
		break;
	case 'B':
		CPX_CALL(CPXsetintparam, env, CPX_PARAM_PREIND, CPX_OFF);
		status = CPX_CALL_SILENT(CPXbaropt, env, lp);
		CPX_CALL(CPXsetintparam, env, CPX_PARAM_PREIND, CPX_ON);
		break;
	case 'A':
	{
		// for the analytic point
		CPX_CALL(CPXsetintparam, env, CPXPARAM_SolutionType, CPX_NONBASIC_SOLN);
		CPX_CALL(CPXsetintparam, env, CPX_PARAM_PREIND, CPX_OFF);
		status = CPX_CALL_SILENT(CPXbaropt, env, lp);
		CPX_CALL(CPXsetintparam, env, CPXPARAM_SolutionType, CPX_BASIC_SOLN);
		CPX_CALL(CPXsetintparam, env, CPX_PARAM_PREIND, CPX_ON);
		break;
	}
	default:
		throw std::runtime_error("Unexpected method for lpopt");
	}

	// consoleError("status: {}", status);
	if (status) // if optimization failed.
		return false;

	// auto solve_status = this->status();
	// // if (method != 'D' && solve_status == CPX_STAT_OPTIMAL_INFEAS) // if found optimal solution, but, when rescaling, the solution values become infeasible, rerun with a more stable method (Dual simplex).
	// // 	status = CPX_CALL_SILENT(CPXdualopt, env, lp);
	// // solve_status = this->status();

	// int primalFeas = 0;
	// int solmethod = 0;
	// int soltype = 0;
	// CPX_CALL(CPXsolninfo, env, lp, &solmethod, &soltype, &primalFeas, nullptr);

	// consoleError("status: {} | sol type: {} | sol method: {} | pfeas: {}", solve_status, soltype, solmethod, primalFeas);

	// consoleError("solve status: {}", solve_status);
	// if (solve_status == CPX_STAT_INFEASIBLE || solve_status == CPX_STAT_INForUNBD)
	// 	return false;
	// consoleError("obj= {:.4f}", this->objval());
	return true;
}

bool CPXModel::isInfeasibleOrTimeReached()
{
	auto solve_status = status();
	// consoleError("solve status: {}", solve_status);
	if (solve_status == CPX_STAT_INFEASIBLE || solve_status == CPX_STAT_INForUNBD || solve_status == CPX_STAT_ABORT_DETTIME_LIM || solve_status == CPX_STAT_ABORT_TIME_LIM || solve_status == CPXMIP_TIME_LIM_INFEAS || solve_status == CPXMIP_DETTIME_LIM_INFEAS || solve_status == CPXMIP_INFEASIBLE || solve_status == CPXMIP_INForUNBD)
		return true;
	return false;
}

bool CPXModel::mipopt()
{
	// every call to solve must be SILENT to avoid throwing error when sub-problems are infeasible, but, in this case, should only be called in stage 3, when problem is known to be feasible already!
	DOMINIQS_ASSERT(env && lp);
	// CPX_CALL(CPXsetintparam, env, CPXPARAM_MIP_Strategy_FPHeur, 1);
	auto status = CPX_CALL_SILENT(CPXmipopt, env, lp);
	auto solve_status = this->status();
	// consoleError("status: {}", solve_status);
	return true;
}

bool CPXModel::presolve()
{
	DOMINIQS_ASSERT(env && lp);
	int status = CPX_CALL_SILENT(CPXpresolve, env, lp, CPX_ALG_NONE);
	if (status)
		return false;
	return true;
}

void CPXModel::postsolve()
{
	DOMINIQS_ASSERT(env && lp);
	// no-op in CPLEX
}

std::vector<double> CPXModel::postsolveSolution(const std::vector<double> &preX) const
{
	// uncrush solution
	int n = ncols();
	std::vector<double> origX(n, 0.0);
	CPX_CALL(CPXuncrushx, env, lp, &origX[0], &preX[0]);
	return origX;
}

std::vector<double> CPXModel::presolveSolution(const std::vector<double> &origX) const
{
	// crush solution
	int n = presolvedModel()->ncols();
	std::vector<double> preX(n, 0.0);
	CPX_CALL(CPXcrushx, env, lp, &origX[0], &preX[0]);
	return preX;
}

/* Get solution */
double CPXModel::objval() const
{
	DOMINIQS_ASSERT(env && lp);
	double ret;
	CPX_CALL(CPXgetobjval, env, lp, &ret);
	return ret;
}

void CPXModel::sol(double *x, int first, int last) const
{
	DOMINIQS_ASSERT(env && lp);
	DOMINIQS_ASSERT((first >= 0) && (first < ncols()));
	if (last == -1)
		last = ncols() - 1;
	DOMINIQS_ASSERT((last >= 0) && (last < ncols()));
	CPX_CALL(CPXgetx, env, lp, x, first, last);
}

void CPXModel::reduced_costs(double *x, int first, int last) const
{
	DOMINIQS_ASSERT(env && lp);
	DOMINIQS_ASSERT((first >= 0) && (first < ncols()));
	if (last == -1)
		last = ncols() - 1;
	DOMINIQS_ASSERT((last >= 0) && (last < ncols()));
	CPX_CALL(CPXgetdj, env, lp, x, first, last);
}

bool CPXModel::isPrimalFeas() const
{
	DOMINIQS_ASSERT(env && lp);
	int primalFeas = 0;
	CPX_CALL(CPXsolninfo, env, lp, nullptr, nullptr, &primalFeas, nullptr);
	return (primalFeas > 0);
}

/* Parameters */
void CPXModel::handleCtrlC(bool flag)
{
	if (flag)
	{
		CPXModel_UserBreak = 0;
		previousHandler = ::signal(SIGINT, userSignalBreak);
		restoreSignalHandler = true;
		CPX_CALL(CPXsetterminate, env, &CPXModel_UserBreak);
	}
	else
	{
		if (restoreSignalHandler)
		{
			::signal(SIGINT, previousHandler);
			restoreSignalHandler = false;
			CPX_CALL(CPXsetterminate, env, nullptr);
		}
	}
}

bool CPXModel::aborted() const
{
	return CPXModel_UserBreak;
}

void CPXModel::seed(int seed)
{
	DOMINIQS_ASSERT(env);
	CPX_CALL(CPXsetintparam, env, CPX_PARAM_RANDOMSEED, seed);
}

void CPXModel::logging(bool log)
{
	DOMINIQS_ASSERT(env);
	if (log)
		CPX_CALL(CPXsetintparam, env, CPX_PARAM_SCRIND, CPX_ON);
	else
		CPX_CALL(CPXsetintparam, env, CPX_PARAM_SCRIND, CPX_OFF);
}

int CPXModel::intParam(IntParam which) const
{
	DOMINIQS_ASSERT(env);
	int value;

	switch (which)
	{
	case IntParam::Threads:
		CPX_CALL(CPXgetintparam, env, CPX_PARAM_THREADS, &value);
		break;
	case IntParam::SolutionLimit:
		CPX_CALL(CPXgetintparam, env, CPX_PARAM_INTSOLLIM, &value);
		break;
	case IntParam::NodeLimit:
		CPX_CALL(CPXgetintparam, env, CPX_PARAM_NODELIM, &value);
		break;
	case IntParam::IterLimit:
		CPX_CALL(CPXgetintparam, env, CPX_PARAM_ITLIM, &value);
		break;
	case IntParam::PdlpWarmStart:
		break;
	case IntParam::Presolve:
		CPX_CALL(CPXgetintparam, env, CPXPARAM_Preprocessing_Presolve, &value);
		break;
	case IntParam::FeasOptMode:
		CPX_CALL(CPXgetintparam, env, CPXPARAM_Feasopt_Mode, &value);
		break;
	case IntParam::Emphasis:
		CPX_CALL(CPXgetintparam, env, CPXPARAM_Emphasis_MIP, &value);
		break;
	default:
		throw std::runtime_error("Unknown integer parameter");
	}

	return (int)value;
}

void CPXModel::intParam(IntParam which, int value)
{
	DOMINIQS_ASSERT(env);

	switch (which)
	{
	case IntParam::Threads:
		CPX_CALL(CPXsetintparam, env, CPX_PARAM_THREADS, value);
		break;
	case IntParam::SolutionLimit:
		CPX_CALL(CPXsetintparam, env, CPX_PARAM_INTSOLLIM, value);
		break;
	case IntParam::NodeLimit:
		CPX_CALL(CPXsetintparam, env, CPX_PARAM_NODELIM, value);
		break;
	case IntParam::IterLimit:
		CPX_CALL(CPXsetintparam, env, CPX_PARAM_ITLIM, value);
		break;
	case IntParam::PdlpWarmStart:
		break;
	case IntParam::Presolve:
		CPX_CALL(CPXsetintparam, env, CPXPARAM_Preprocessing_Presolve, value);
		break;
	case IntParam::FeasOptMode:
		CPX_CALL(CPXsetintparam, env, CPXPARAM_Feasopt_Mode, value);
		break;
	case IntParam::Emphasis:
		CPX_CALL(CPXsetintparam, env, CPXPARAM_Emphasis_MIP, value);
		break;
	default:
		throw std::runtime_error("Unknown integer parameter ");
	}
}

double CPXModel::dblParam(DblParam which) const
{
	DOMINIQS_ASSERT(env);
	double value;

	switch (which)
	{
	case DblParam::TimeLimit:
		CPX_CALL(CPXgetdblparam, env, CPX_PARAM_TILIM, &value);
		break;
	case DblParam::FeasibilityTolerance:
		CPX_CALL(CPXgetdblparam, env, CPX_PARAM_EPRHS, &value);
		break;
	case DblParam::IntegralityTolerance:
		CPX_CALL(CPXgetdblparam, env, CPX_PARAM_EPINT, &value);
		break;
	case DblParam::WorkMem:
		CPX_CALL(CPXgetdblparam, env, CPX_PARAM_WORKMEM, &value);
		break;
	default:
		throw std::runtime_error("Unknown double parameter");
	}

	return value;
}

void CPXModel::dblParam(DblParam which, double value)
{
	DOMINIQS_ASSERT(env);
	switch (which)
	{
	case DblParam::TimeLimit:
		CPX_CALL(CPXsetdblparam, env, CPX_PARAM_TILIM, value);
		break;
	case DblParam::FeasibilityTolerance:
		CPX_CALL(CPXsetdblparam, env, CPX_PARAM_EPRHS, value);
		break;
	case DblParam::IntegralityTolerance:
		CPX_CALL(CPXsetdblparam, env, CPX_PARAM_EPINT, value);
		break;
	case DblParam::WorkMem:
		CPX_CALL(CPXsetdblparam, env, CPX_PARAM_WORKMEM, value);
		break;
	default:
		throw std::runtime_error("Unknown double parameter");
	}
}

int CPXModel::intAttr(IntAttr which) const
{
	DOMINIQS_ASSERT(env && lp);
	int value = 0;

	switch (which)
	{
	case IntAttr::Nodes:
		value = CPXgetnodecnt(env, lp);
		break;
	case IntAttr::NodesLeft:
		value = CPXgetnodeleftcnt(env, lp);
		break;
	case IntAttr::BarrierIterations:
		value = CPXgetbaritcnt(env, lp);
		break;
	case IntAttr::SimplexIterations:
		value = CPXgetitcnt(env, lp);
		break;
	case IntAttr::PDLPIterations:
		break;
	default:
		throw std::runtime_error("Unknown integer attribute");
	}

	return value;
}

double CPXModel::dblAttr(DblAttr which) const
{
	DOMINIQS_ASSERT(env && lp);
	double value;

	switch (which)
	{
	case DblAttr::MIPDualBound:
		CPX_CALL(CPXgetbestobjval, env, lp, &value);
		break;
	default:
		throw std::runtime_error("Unknown double attribute");
	}

	return value;
}

void CPXModel::terminationReason(std::string &reason)
{
	// ToDo
	reason = "-";
}

/* Access model data */
int CPXModel::nrows() const
{
	DOMINIQS_ASSERT(env && lp);
	return CPXgetnumrows(env, lp);
}

int CPXModel::ncols() const
{
	DOMINIQS_ASSERT(env && lp);
	return CPXgetnumcols(env, lp);
}

int CPXModel::nnz() const
{
	int nnz;
	int tmp = 0;
	CPXgetrows(env, lp, &tmp, nullptr, nullptr, nullptr, 0, &nnz, 0, nrows() - 1);
	DOMINIQS_ASSERT(nnz <= 0);
	return -nnz;
}

double CPXModel::objOffset() const
{
	DOMINIQS_ASSERT(env && lp);
	double objOffset = 0.0;
	CPX_CALL(CPXgetobjoffset, env, lp, &objOffset);
	return objOffset;
}

ObjSense CPXModel::objSense() const
{
	DOMINIQS_ASSERT(env && lp);
	int cpxobjsen = CPXgetobjsen(env, lp);
	return (cpxobjsen > 0) ? ObjSense::MIN : ObjSense::MAX;
}

void CPXModel::lbs(double *lb, int first, int last) const
{
	DOMINIQS_ASSERT(env && lp);
	DOMINIQS_ASSERT((first >= 0) && (first < ncols()));
	if (last == -1)
		last = ncols() - 1;
	DOMINIQS_ASSERT((last >= 0) && (last < ncols()));
	DOMINIQS_ASSERT(first <= last);
	CPX_CALL(CPXgetlb, env, lp, lb, first, last);
}

void CPXModel::ubs(double *ub, int first, int last) const
{
	DOMINIQS_ASSERT(env && lp);
	DOMINIQS_ASSERT((first >= 0) && (first < ncols()));
	if (last == -1)
		last = ncols() - 1;
	DOMINIQS_ASSERT((last >= 0) && (last < ncols()));
	DOMINIQS_ASSERT(first <= last);
	CPX_CALL(CPXgetub, env, lp, ub, first, last);
}

void CPXModel::objcoefs(double *obj, int first, int last) const
{
	DOMINIQS_ASSERT(env && lp);
	DOMINIQS_ASSERT((first >= 0) && (first < ncols()));
	if (last == -1)
		last = ncols() - 1;
	DOMINIQS_ASSERT((last >= 0) && (last < ncols()));
	DOMINIQS_ASSERT(first <= last);
	CPX_CALL(CPXgetobj, env, lp, obj, first, last);
}

void CPXModel::ctypes(char *ctype, int first, int last) const
{
	DOMINIQS_ASSERT(env && lp);
	DOMINIQS_ASSERT((first >= 0) && (first < ncols()));
	if (last == -1)
		last = ncols() - 1;
	DOMINIQS_ASSERT((last >= 0) && (last < ncols()));
	DOMINIQS_ASSERT(first <= last);
	CPX_CALL(CPXgetctype, env, lp, ctype, first, last);
}

void CPXModel::sense(char *sense, int first, int last) const
{
	DOMINIQS_ASSERT(env && lp);
	DOMINIQS_ASSERT((first >= 0) && (first < nrows()));
	if (last == -1)
		last = nrows() - 1;
	DOMINIQS_ASSERT((last >= 0) && (last < nrows()));
	DOMINIQS_ASSERT(first <= last);
	CPX_CALL(CPXgetsense, env, lp, sense, first, last);
}

void CPXModel::rhs(double *rhs, int first, int last) const
{
	DOMINIQS_ASSERT(env && lp);
	DOMINIQS_ASSERT((first >= 0) && (first < nrows()));
	if (last == -1)
		last = nrows() - 1;
	DOMINIQS_ASSERT((last >= 0) && (last < nrows()));
	DOMINIQS_ASSERT(first <= last);
	CPX_CALL(CPXgetrhs, env, lp, rhs, first, last);
}

void CPXModel::range(double *range, int first, int last) const
{
	DOMINIQS_ASSERT(env && lp);
	DOMINIQS_ASSERT((first >= 0) && (first < nrows()));
	if (last == -1)
		last = nrows() - 1;
	DOMINIQS_ASSERT((last >= 0) && (last < nrows()));
	DOMINIQS_ASSERT(first <= last);
	CPX_CALL(CPXgetrngval, env, lp, range, first, last);
}

void CPXModel::row(int ridx, dominiqs::SparseVector &row, char &sense, double &rhs, double &rngval) const
{
	DOMINIQS_ASSERT(env && lp);
	DOMINIQS_ASSERT((ridx >= 0) && (ridx < nrows()));
	// get row nz
	int tmp = 0;
	int size;
	CPXgetrows(env, lp, &tmp, &tmp, 0, 0, 0, &size, ridx, ridx);
	// get coef
	size = -size;
	if (size)
	{
		row.resize(size);
		CPX_CALL(CPXgetrows, env, lp, &tmp, &tmp, row.idx(), row.coef(), size, &tmp, ridx, ridx);
	}
	else
	{
		row.clear();
	}
	// here to correctly handle empty constraints
	// get rhs
	CPX_CALL(CPXgetrhs, env, lp, &rhs, ridx, ridx);
	// get sense
	CPX_CALL(CPXgetsense, env, lp, &sense, ridx, ridx);
	// get rngval
	CPX_CALL(CPXgetrngval, env, lp, &rngval, ridx, ridx);
	// CPLEX treats ranged rows considering the constraint satisfied
	// if the linear expression is in the range [rhs, rhs+rngval].
	// However, we interpret ranged rows differently, and the allowed
	// range is [rhs-rngval,rhs] (in both cases, rngval >= 0).
	// So we might have to update the rhs
	if (sense == 'R')
	{
		DOMINIQS_ASSERT(rngval >= 0.0);
		rhs += rngval;
	}
}

void CPXModel::rows(dominiqs::SparseMatrix &matrix) const
{
	DOMINIQS_ASSERT(env && lp);
	// get nnz
	int tmp = 0;
	int beg = 0;
	int end = CPXgetnumrows(env, lp);
	int size;
	matrix.matbeg.resize(end);
	matrix.k = end;
	CPXgetrows(env, lp, &tmp, &(matrix.matbeg[0]), 0, 0, 0, &size, beg, end - 1);
	// get coef
	size = -size;
	if (size)
	{
		matrix.nnz = size;
		matrix.matind.resize(matrix.nnz);
		matrix.matval.resize(matrix.nnz);
		CPX_CALL(CPXgetrows, env, lp, &tmp,
				 &(matrix.matbeg[0]), &(matrix.matind[0]), &(matrix.matval[0]),
				 size, &tmp, beg, end - 1);
	}
	else
	{
		matrix.nnz = 0;
		matrix.matind.clear();
		matrix.matval.clear();
	}
}

void CPXModel::col(int cidx, dominiqs::SparseVector &col, char &type, double &lb, double &ub, double &obj) const
{
	DOMINIQS_ASSERT(env && lp);
	DOMINIQS_ASSERT((cidx >= 0) && (cidx < ncols()));
	// get col nz
	int tmp = 0;
	int size;
	CPXgetcols(env, lp, &tmp, &tmp, 0, 0, 0, &size, cidx, cidx);
	// get coef
	size = -size;
	if (size)
	{
		col.resize(size);
		CPX_CALL(CPXgetcols, env, lp, &tmp, &tmp, col.idx(), col.coef(), size, &tmp, cidx, cidx);
	}
	else
	{
		col.clear();
	}
	// here to correctly handle empty vars
	// get bounds
	CPX_CALL(CPXgetlb, env, lp, &lb, cidx, cidx);
	CPX_CALL(CPXgetub, env, lp, &ub, cidx, cidx);
	// get obj
	CPX_CALL(CPXgetobj, env, lp, &obj, cidx, cidx);
	// get type
	int status = CPXgetctype(env, lp, &type, cidx, cidx);
	if (status)
		type = 'C'; // cannot call CPXgetctype on an LP
}

void CPXModel::cols(dominiqs::SparseMatrix &matrix) const
{
	DOMINIQS_ASSERT(env && lp);
	// get nnz
	int tmp = 0;
	int beg = 0;
	int end = CPXgetnumcols(env, lp);
	int size;
	matrix.matbeg.resize(end);
	matrix.k = end;
	CPXgetcols(env, lp, &tmp, &(matrix.matbeg[0]), 0, 0, 0, &size, beg, end - 1);
	// get coef
	size = -size;
	if (size)
	{
		matrix.nnz = size;
		matrix.matind.resize(matrix.nnz);
		matrix.matval.resize(matrix.nnz);
		CPX_CALL(CPXgetcols, env, lp, &tmp,
				 &(matrix.matbeg[0]), &(matrix.matind[0]), &(matrix.matval[0]),
				 size, &tmp, beg, end - 1);
	}
	else
	{
		matrix.nnz = 0;
		matrix.matind.clear();
		matrix.matval.clear();
	}
}

void CPXModel::colNames(std::vector<std::string> &names, int first, int last) const
{
	DOMINIQS_ASSERT(env && lp);
	DOMINIQS_ASSERT((first >= 0) && (first < ncols()));
	if (last == -1)
		last = ncols() - 1;
	DOMINIQS_ASSERT((last >= 0) && (last < ncols()));
	DOMINIQS_ASSERT(first <= last);
	names.clear();
	int count = last - first + 1;
	std::vector<char> buffer;
	std::vector<char *> cnames(count, nullptr);
	int surplus;
	CPXgetcolname(env, lp, &cnames[0], nullptr, 0, &surplus, first, last);
	if (surplus)
	{
		buffer.resize(-surplus);
		CPX_CALL(CPXgetcolname, env, lp, &cnames[0], &buffer[0], buffer.size(), &surplus, first, last);
		for (int i = 0; i < count; i++)
			names.push_back(std::string(cnames[i]));
	}
	else
	{
		// no names
		for (int i = 0; i < count; i++)
			names.push_back("");
	}
}

void CPXModel::rowNames(std::vector<std::string> &names, int first, int last) const
{
	DOMINIQS_ASSERT(env && lp);
	DOMINIQS_ASSERT((first >= 0) && (first < nrows()));
	if (last == -1)
		last = nrows() - 1;
	DOMINIQS_ASSERT((last >= 0) && (last < nrows()));
	DOMINIQS_ASSERT(first <= last);
	names.clear();
	int count = last - first + 1;
	std::vector<char> buffer;
	std::vector<char *> rnames(count, 0);
	int surplus;
	CPXgetrowname(env, lp, &rnames[0], nullptr, 0, &surplus, first, last);
	if (surplus)
	{
		buffer.resize(-surplus);
		CPX_CALL(CPXgetrowname, env, lp, &rnames[0], &buffer[0], buffer.size(), &surplus, first, last);
		for (int i = 0; i < count; i++)
			names.push_back(std::string(rnames[i]));
	}
	else
	{
		// no names
		for (int i = 0; i < count; i++)
			names.push_back("");
	}
}

/* Data modifications */
void CPXModel::addEmptyCol(const std::string &name, char ctype, double lb, double ub, double obj)
{
	DOMINIQS_ASSERT(env && lp);
	char *cname = (char *)(name.c_str());
	const char *ctypeptr = (ctype == 'C') ? nullptr : &ctype; //< do not risk turning the model into a MIP
	CPX_CALL(CPXnewcols, env, lp, 1, &obj, &lb, &ub, ctypeptr, &cname);
}

void CPXModel::addCol(const std::string &name, const int *idx, const double *val, int cnt, char ctype, double lb, double ub, double obj)
{
	DOMINIQS_ASSERT(env && lp);
	int matbeg = 0;
	char *cname = (char *)(name.c_str());
	if (cnt > 0)
	{
		DOMINIQS_ASSERT(idx && val);
		CPX_CALL(CPXaddcols, env, lp, 1, cnt, &obj, &matbeg, idx, val, &lb, &ub, &cname);
		if (ctype != 'C')
		{
			int last = ncols() - 1;
			CPX_CALL(CPXchgctype, env, lp, 1, &last, &ctype);
		}
	}
	else
		CPX_CALL(CPXnewcols, env, lp, 1, &obj, &lb, &ub, &ctype, &cname);
}

void CPXModel::addRow(const std::string &name, const int *idx, const double *val, int cnt, char sense, double rhs, double rngval)
{
	DOMINIQS_ASSERT(env && lp);
	int matbeg = 0;
	char *rname = (char *)(name.c_str());
	if (sense == 'R')
	{
		DOMINIQS_ASSERT(rngval >= 0.0);
		// for ranged rows, we assue [rhs-rngval,rhs] while CPLEX uses [rhs, rhs+rngval]
		rhs -= rngval;
	}
	CPX_CALL(CPXaddrows, env, lp, 0, 1, cnt, &rhs, &sense, &matbeg, idx, val, nullptr, &rname);
	if (sense == 'R')
	{
		DOMINIQS_ASSERT(rngval >= 0.0);
		int ridx = nrows() - 1;
		DOMINIQS_ASSERT(ridx >= 0);
		CPX_CALL(CPXchgrngval, env, lp, 1, &ridx, &rngval);
	}
}

void CPXModel::delRow(int ridx)
{
	DOMINIQS_ASSERT(env && lp);
	CPX_CALL(CPXdelrows, env, lp, ridx, ridx);
}

void CPXModel::delCol(int cidx)
{
	DOMINIQS_ASSERT(env && lp);
	CPX_CALL(CPXdelcols, env, lp, cidx, cidx);
}

void CPXModel::delRows(int first, int last)
{
	DOMINIQS_ASSERT(env && lp);
	DOMINIQS_ASSERT((first >= 0) && (first < nrows()));
	DOMINIQS_ASSERT((last >= 0) && (last < nrows()));
	DOMINIQS_ASSERT(first <= last);
	CPX_CALL(CPXdelrows, env, lp, first, last);
}

void CPXModel::delCols(int first, int last)
{
	DOMINIQS_ASSERT(env && lp);
	DOMINIQS_ASSERT((first >= 0) && (first < ncols()));
	DOMINIQS_ASSERT((last >= 0) && (last < ncols()));
	DOMINIQS_ASSERT(first <= last);
	CPX_CALL(CPXdelcols, env, lp, first, last);
}

void CPXModel::objSense(ObjSense objsen)
{
	DOMINIQS_ASSERT(env && lp);
	CPXchgobjsen(env, lp, static_cast<int>(objsen));
}

void CPXModel::objOffset(double val)
{
	DOMINIQS_ASSERT(env && lp);
	CPX_CALL(CPXchgobjoffset, env, lp, val);
}

void CPXModel::lb(int cidx, double val)
{
	DOMINIQS_ASSERT(env && lp);
	DOMINIQS_ASSERT((cidx >= 0) && (cidx < ncols()));
	char lu = 'L';
	CPX_CALL(CPXchgbds, env, lp, 1, &cidx, &lu, &val);
}

void CPXModel::lbs(int cnt, const int *cols, const double *values)
{
	DOMINIQS_ASSERT(env && lp);
	std::vector<char> lu(cnt, 'L');
	CPX_CALL(CPXchgbds, env, lp, cnt, cols, &lu[0], values);
}

void CPXModel::ub(int cidx, double val)
{
	DOMINIQS_ASSERT(env && lp);
	DOMINIQS_ASSERT((cidx >= 0) && (cidx < ncols()));
	char lu = 'U';
	CPX_CALL(CPXchgbds, env, lp, 1, &cidx, &lu, &val);
}

void CPXModel::ubs(int cnt, const int *cols, const double *values)
{
	DOMINIQS_ASSERT(env && lp);
	std::vector<char> lu(cnt, 'U');
	CPX_CALL(CPXchgbds, env, lp, cnt, cols, &lu[0], values);
}

void CPXModel::fixCol(int cidx, double val)
{
	DOMINIQS_ASSERT(env && lp);
	DOMINIQS_ASSERT((cidx >= 0) && (cidx < ncols()));
	char lu = 'B';
	CPX_CALL(CPXchgbds, env, lp, 1, &cidx, &lu, &val);
}

void CPXModel::objcoef(int cidx, double val)
{
	DOMINIQS_ASSERT(env && lp);
	DOMINIQS_ASSERT((cidx >= 0) && (cidx < ncols()));
	CPX_CALL(CPXchgobj, env, lp, 1, &cidx, &val);
}

void CPXModel::objcoefs(int cnt, const int *cols, const double *values)
{
	DOMINIQS_ASSERT(env && lp);
	CPX_CALL(CPXchgobj, env, lp, cnt, cols, values);
}

void CPXModel::ctype(int cidx, char val)
{
	DOMINIQS_ASSERT(env && lp);
	DOMINIQS_ASSERT((cidx >= 0) && (cidx < ncols()));
	DOMINIQS_ASSERT((val == 'B') || (val == 'I') || (val == 'C'));
	CPX_CALL(CPXchgctype, env, lp, 1, &cidx, &val);
}

void CPXModel::ctypes(int cnt, const int *cols, const char *values)
{
	DOMINIQS_ASSERT(env && lp);
	CPX_CALL(CPXchgctype, env, lp, cnt, cols, values);
}

void CPXModel::switchToLP()
{
	DOMINIQS_ASSERT(env && lp);
	CPX_CALL(CPXchgprobtype, env, lp, CPXPROB_LP);
}

void CPXModel::switchToMIP()
{
}

/* Private interface */
CPXModel *CPXModel::clone_impl() const
{
	DOMINIQS_ASSERT(env && lp);
	int status = 0;
	CPXLPptr cloned = CPXcloneprob(env, lp, &status);
	if (status)
		throwCplexError(env, status);
	return new CPXModel(env, cloned, false, true);
}

CPXModel *CPXModel::presolvedmodel_impl() const
{
	int preStat;
	CPX_CALL(CPXgetprestat, env, lp, &preStat, nullptr, nullptr, nullptr, nullptr);

	// consoleError("Presolve status: {}", preStat);
	if ((preStat == 0) /* no presolve */ || (preStat == 2) /* reduced to empty problem*/)
	{
		if (preStat == 2) // reduced to empty problem
			consoleWarn("Solved in presolve (yields empty presolved problem). Presolve undone.");
		// nothing to do
		return nullptr;
	}
	else // return clone of presolved problem
	{
		CPXCLPptr redlp = nullptr;
		CPX_CALL(CPXgetredlp, env, lp, &redlp);
		int status = 0;
		CPXLPptr cloned = CPXcloneprob(env, redlp, &status);
		if (status)
			throwCplexError(env, status);
		CPXModel *premodel = new CPXModel(env, cloned, false, true);
		return premodel;
	}
	DOMINIQS_ASSERT(false);
	return nullptr;
}

void CPXModel::updateModelVarBounds(std::optional<boost::dynamic_bitset<>> vars_entering_problem, std::optional<boost::dynamic_bitset<>> vars_leaving_problem)
{
	// activate new variables to kernel.
	if (vars_entering_problem)
	{
		// std::cout << "entering: ";
		for (int var_index = vars_entering_problem->find_first(); var_index != boost::dynamic_bitset<>::npos; var_index = vars_entering_problem->find_next(var_index))
		{
			ub(var_index, 1.0);
			// curr_active_binary_vars_[var_index] = 1;
			// std::cout << var_index << " ";
		}
		// std::cout << std::endl;
	}

	// deactivate variables leaving kernel.
	if (vars_leaving_problem)
	{
		// std::cout << "leaving: ";
		for (int var_index = vars_leaving_problem->find_first(); var_index != boost::dynamic_bitset<>::npos; var_index = vars_leaving_problem->find_next(var_index))
		{
			ub(var_index, 0.0);
			// curr_active_binary_vars_[var_index] = 0;
			// std::cout << (*vars_)[var_index] << " ";
		}
		// std::cout << std::endl;
	}
	// cplex_->exportModel("update.lp");
}

void CPXModel::findSetOfConflictingVariables(boost::dynamic_bitset<> inactive_binary_vars, std::vector<int> &conflicting_constraints, std::vector<int> &conflicting_vars, bool optimize_set, double time_left)
{
	DOMINIQS_ASSERT(env && lp);
	conflicting_vars.clear();
	conflicting_constraints.clear();
	int numCols = ncols();
	dblParam(DblParam::TimeLimit, time_left);

	consoleLog(" * Attempt to make kernel LP feasible");
	if (optimize_set) // if you need to optimize the subset of conflicting constraints/vars (e.g., make it as small as possible).
	{
		intParam(IntParam::FeasOptMode, CPX_FEASOPT_MIN_SUM); // maybe also test CPX_FEASOPT_OPT_SUM. All other options make the problem an MIP! Too expensive.
		// only allow to relax (increase) the upper bounds of the binary variables that DO NOT belong to the kernel yet.
		// forbit the relaxation of the constraints bounds.
		std::vector<double> candidate_vars_with_upper_bounds_to_be_relaxed(numCols, 0);
		for (int var_index = inactive_binary_vars.find_first(); var_index != boost::dynamic_bitset<>::npos; var_index = inactive_binary_vars.find_next(var_index))
			candidate_vars_with_upper_bounds_to_be_relaxed[var_index] = 1;

		consoleLog(" * {} inactive binary vars considered for entering kernel (i.e., being activated: ub -> 1)", inactive_binary_vars.count());

		auto result = CPXfeasopt(env, lp, NULL, NULL, NULL, &(candidate_vars_with_upper_bounds_to_be_relaxed[0]));
		if (result == 0) // success
		{
			// if found a way to making the problem feasible.
			auto status = CPXgetstat(env, lp);
			// consoleWarn(" ststus: {}", status);
			if (status == CPX_STAT_FEASIBLE_RELAXED_SUM || status == CPX_STAT_OPTIMAL_RELAXED_SUM)
			{
				// consoleInfo("Found possible LP feasible model");
				// std::vector<double> solution(numCols);
				// sol(&(solution[0]));
				std::vector<double> infeasout(numCols);

				status = CPXgetcolinfeas(env, lp, NULL, &(infeasout[0]), 0, numCols - 1);
				for (int var_index = inactive_binary_vars.find_first(); var_index != boost::dynamic_bitset<>::npos; var_index = inactive_binary_vars.find_next(var_index))
					if (greaterThan(infeasout[var_index], 0))
						conflicting_vars.push_back(var_index);
			}
		}
	}
	else // if only need to find an arbitrary minimal set.
	{
		intParam(IntParam::Presolve, 0); // disable presolve to be able to get infeasibility info.
		// find a minimal conflict subset of constraints/variables.
		auto result = CPXrefineconflictext(env, lp, 0, 0, NULL, NULL, NULL, NULL);
		if (result == 0) // success
		{
			auto status = CPXgetstat(env, lp);
			// consoleWarn(" ststus: {}", status);
			if (status == CPX_STAT_CONFLICT_FEASIBLE || status == CPX_STAT_CONFLICT_MINIMAL)
			{
				int num_conflicting_cols = 0, num_conflicting_rows = 0;
				// retrieve the mininal conflict subset found.
				CPXgetconflict(env, lp, 0, NULL, NULL, &num_conflicting_rows, NULL, NULL, &num_conflicting_cols);

				conflicting_constraints.resize(num_conflicting_rows);
				std::vector<int> all_conflicting_vars(num_conflicting_cols);
				std::vector<int> row_stat(num_conflicting_rows);
				std::vector<int> col_stat(num_conflicting_cols);

				CPXgetconflict(env, lp, 0, &(conflicting_constraints[0]), &(row_stat[0]), &num_conflicting_rows, &(all_conflicting_vars[0]), &(col_stat[0]), &num_conflicting_cols);

				std::vector<std::string> col_names;
				colNames(col_names);
				int num_vars = this->ncols();
				std::vector<double> ubs(num_vars), lbs(num_vars);
				this->lbs(&(lbs[0]));
				this->ubs(&(ubs[0]));
				// filter out the columns that don't necessarily have an "upper bound" conflict.
				for (int i = 0; i < num_conflicting_cols; ++i)
				{
					if (col_stat[i] == CPX_CONFLICT_UB || col_stat[i] == CPX_CONFLICT_MEMBER || col_stat[i] == CPX_CONFLICT_POSSIBLE_UB || col_stat[i] == CPX_CONFLICT_POSSIBLE_MEMBER) // CPX_CONFLICT_POSSIBLE_MEMBER stands for a conflict of both lower and upper bounds!
					{
						conflicting_vars.push_back(all_conflicting_vars[i]);
						// consoleWarn("{}: {} [{},{}]", col_names[i], col_stat[i], lbs[i], ubs[i]);
					}
					// else
					// 	consoleWarn("{}: {} [{},{}]", col_names[i], col_stat[i], lbs[i], ubs[i]);
				}
			}
		}
	}
}