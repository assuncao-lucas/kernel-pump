/**
 * @file grbmodel.cpp
 * @brief Implementation of MIPModelI for Gurobi (C API)
 *
 * See include/kernelpump/grbmodel.h and the Gurobi backend spec for the design.
 */

#include "kernelpump/grbmodel.h"
#include <signal.h>
#include <algorithm>
#include <climits>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <unordered_map>
#include <fmt/format.h>
#include <utils/consolelog.h>
#include <utils/floats.h>

namespace
{
	volatile sig_atomic_t GRBModel_UserBreak = 0;
	GRBmodel *volatile GRBModel_ActiveModel = nullptr;

	void userSignalBreak(int)
	{
		GRBModel_UserBreak = 1;
		if (GRBModel_ActiveModel)
			GRBterminate(GRBModel_ActiveModel); // documented as safe to call asynchronously
	}

	[[noreturn]] void throwGurobiError(GRBenv *env, int error, const char *what)
	{
		const char *msg = env ? GRBgeterrormsg(env) : nullptr;
		throw std::runtime_error(fmt::format("Gurobi error {} in {}: {}", error, what, (msg && *msg) ? msg : "(no message)"));
	}

/* Call a Gurobi API function and throw on a non-zero return code. */
#define GRB_CALL(env, call)                         \
	do                                              \
	{                                               \
		int grb_err_ = (call);                      \
		if (grb_err_)                               \
			throwGurobiError((env), grb_err_, #call); \
	} while (0)

	/* Infinity translation: the code base uses INFBOUND (1e20), Gurobi GRB_INFINITY (1e100). */
	inline double toGrb(double v)
	{
		if (v >= INFBOUND)
			return GRB_INFINITY;
		if (v <= -INFBOUND)
			return -GRB_INFINITY;
		return v;
	}

	inline double fromGrb(double v)
	{
		if (v >= INFBOUND)
			return INFBOUND;
		if (v <= -INFBOUND)
			return -INFBOUND;
		return v;
	}

	inline char senseFromGrb(char s)
	{
		switch (s)
		{
		case GRB_LESS_EQUAL:
			return 'L';
		case GRB_GREATER_EQUAL:
			return 'G';
		case GRB_EQUAL:
			return 'E';
		default:
			throw std::runtime_error(fmt::format("Unexpected Gurobi constraint sense '{}'", s));
		}
	}

	inline char senseToGrb(char s)
	{
		switch (s)
		{
		case 'L':
			return GRB_LESS_EQUAL;
		case 'G':
			return GRB_GREATER_EQUAL;
		case 'E':
			return GRB_EQUAL;
		default:
			throw std::runtime_error(fmt::format("Unsupported constraint sense '{}' for Gurobi", s));
		}
	}

	void setProcessEnv(const char *name, const std::string &value)
	{
#ifdef _WIN32
		if (_putenv_s(name, value.c_str()) != 0)
#else
		if (setenv(name, value.c_str(), 1) != 0)
#endif
			throw std::runtime_error(fmt::format("Could not set environment variable {}", name));
	}

	/* One Gurobi environment per process, created on first use and shared by every model. */
	std::shared_ptr<GRBenv> acquireMasterEnv(const std::string &licenseFile, bool requireLicenseFile)
	{
		static std::weak_ptr<GRBenv> cached;
		if (auto existing = cached.lock())
			return existing;

		if (licenseFile.empty())
		{
			if (requireLicenseFile)
				throw std::runtime_error(
					"solver=grb requires grb.licenseFile=<path to the license to use>. "
					"Refusing to start Gurobi with whatever license the user environment points to "
					"(set grb.requireLicenseFile=0 to allow that).");
		}
		else
		{
			// Must happen before GRBstartenv: Gurobi reads GRB_LICENSE_FILE when the env starts.
			setProcessEnv("GRB_LICENSE_FILE", licenseFile);
		}

		GRBenv *raw = nullptr;
		int error = GRBemptyenv(&raw);
		if (error || !raw)
			throwGurobiError(raw, error ? error : -1, "GRBemptyenv");
		std::shared_ptr<GRBenv> env(raw, [](GRBenv *e)
									{ GRBfreeenv(e); });
		// quiet by default, like CPLEX (screen output off); logging(true) turns it on.
		GRB_CALL(env.get(), GRBsetintparam(env.get(), "OutputFlag", 0));
		GRB_CALL(env.get(), GRBstartenv(env.get()));

		int major = 0, minor = 0, tech = 0;
		GRBversion(&major, &minor, &tech);
		const char *lic = std::getenv("GRB_LICENSE_FILE");
		consoleLog("Gurobi {}.{}.{} started; license file: {}", major, minor, tech, lic ? lic : "(Gurobi default search)");

		cached = env;
		return env;
	}
} // namespace

GRBModel::GRBModel(const std::string &licenseFile, bool requireLicenseFile)
	: env(acquireMasterEnv(licenseFile, requireLicenseFile))
{
	GRB_CALL(env.get(), GRBnewmodel(env.get(), &model, "", 0, nullptr, nullptr, nullptr, nullptr, nullptr));
}

GRBModel::GRBModel(std::shared_ptr<GRBenv> _env, GRBmodel *_model) : env(std::move(_env)), model(_model)
{
	DOMINIQS_ASSERT(env && model);
}

GRBModel::~GRBModel()
{
	if (restoreSignalHandler)
		handleCtrlC(false);
	if (model)
		GRBfreemodel(model);
}

/* Helpers */
void GRBModel::flush() const
{
	if (dirty)
	{
		GRB_CALL(modelEnv(), GRBupdatemodel(model));
		dirty = false;
	}
}

int GRBModel::checkRange(int &first, int &last, int size) const
{
	DOMINIQS_ASSERT((first >= 0) && (first < size));
	if (last == -1)
		last = size - 1;
	DOMINIQS_ASSERT((last >= 0) && (last < size));
	DOMINIQS_ASSERT(first <= last);
	return last - first + 1;
}

void GRBModel::setIntParamEverywhere(const char *name, int value)
{
	// CPLEX params live on the env shared by all models; Gurobi copies params into each model.
	// Setting both keeps params set before readModel() (main.cpp does this) and in clones.
	GRB_CALL(modelEnv(), GRBsetintparam(modelEnv(), name, value));
	GRB_CALL(env.get(), GRBsetintparam(env.get(), name, value));
}

void GRBModel::setDblParamEverywhere(const char *name, double value)
{
	GRB_CALL(modelEnv(), GRBsetdblparam(modelEnv(), name, value));
	GRB_CALL(env.get(), GRBsetdblparam(env.get(), name, value));
}

int GRBModel::getIntAttrOr(const char *name, int fallback) const
{
	flush();
	int value = 0;
	return GRBgetintattr(model, name, &value) ? fallback : value;
}

double GRBModel::getDblAttrOr(const char *name, double fallback) const
{
	flush();
	double value = 0.0;
	return GRBgetdblattr(model, name, &value) ? fallback : value;
}

void GRBModel::getDblArray(const char *attr, double *out, int first, int last, int size, bool mapInf) const
{
	flush();
	int count = checkRange(first, last, size);
	GRB_CALL(modelEnv(), GRBgetdblattrarray(model, attr, first, count, out));
	if (mapInf)
		for (int i = 0; i < count; ++i)
			out[i] = fromGrb(out[i]);
}

void GRBModel::setDblList(const char *attr, int cnt, const int *cols, const double *values, bool mapInf)
{
	if (cnt <= 0)
		return;
	std::vector<int> ind(cols, cols + cnt);
	std::vector<double> val(values, values + cnt);
	if (mapInf)
		for (auto &v : val)
			v = toGrb(v);
	GRB_CALL(modelEnv(), GRBsetdblattrlist(model, attr, cnt, ind.data(), val.data()));
	markDirty();
}

void GRBModel::getNames(const char *attr, std::vector<std::string> &names, int first, int last, int size) const
{
	flush();
	int count = checkRange(first, last, size);
	std::vector<char *> raw(count, nullptr);
	GRB_CALL(modelEnv(), GRBgetstrattrarray(model, attr, first, count, raw.data()));
	names.clear();
	names.reserve(count);
	for (int i = 0; i < count; ++i)
		names.emplace_back(raw[i] ? raw[i] : "");
}

void GRBModel::getMatrix(bool byRow, dominiqs::SparseMatrix &matrix) const
{
	flush();
	int len = byRow ? nrows() : ncols();
	auto query = byRow ? GRBgetconstrs : GRBgetvars;
	matrix.k = len;
	matrix.matbeg.assign(len, 0);
	if (len == 0)
	{
		matrix.nnz = 0;
		matrix.matind.clear();
		matrix.matval.clear();
		return;
	}
	int numnz = 0;
	GRB_CALL(modelEnv(), query(model, &numnz, nullptr, nullptr, nullptr, 0, len));
	matrix.nnz = numnz;
	matrix.matind.resize(numnz);
	matrix.matval.resize(numnz);
	GRB_CALL(modelEnv(), query(model, &numnz, matrix.matbeg.data(),
							   numnz ? matrix.matind.data() : nullptr,
							   numnz ? matrix.matval.data() : nullptr, 0, len));
}

void GRBModel::getVector(bool byRow, int index, dominiqs::SparseVector &vec) const
{
	flush();
	auto query = byRow ? GRBgetconstrs : GRBgetvars;
	int numnz = 0;
	GRB_CALL(modelEnv(), query(model, &numnz, nullptr, nullptr, nullptr, index, 1));
	if (numnz > 0)
	{
		int beg = 0;
		vec.resize(numnz);
		GRB_CALL(modelEnv(), query(model, &numnz, &beg, vec.idx(), vec.coef(), index, 1));
	}
	else
		vec.clear();
}

/* Read/Write */
void GRBModel::readModel(const std::string &filename)
{
	DOMINIQS_ASSERT(env);
	GRBmodel *loaded = nullptr;
	// the new model copies the master env params, which the param setters keep in sync
	GRB_CALL(env.get(), GRBreadmodel(env.get(), filename.c_str(), &loaded));
	if (model)
		GRBfreemodel(model);
	model = loaded;
	dirty = false;
	constraints = nullptr; // cached data of the previous model, if any
	dependency = nullptr;
	int numRangeVars = 0;
	// Gurobi models ranged rows (MPS RANGES) as an equality plus an extra bounded column named RgXXX.
	std::vector<std::string> names;
	if (ncols() > 0)
	{
		colNames(names);
		for (const auto &n : names)
			if (n.rfind("Rg", 0) == 0)
				++numRangeVars;
	}
	if (numRangeVars)
		consoleWarn("Gurobi: model may contain {} range slack column(s) (names starting with 'Rg')", numRangeVars);
}

void GRBModel::writeModel(const std::string &filename, const std::string &format) const
{
	flush();
	// Gurobi picks the format from the extension; an explicit format is appended as extension.
	std::string target = format.empty() ? filename : filename + "." + format;
	GRB_CALL(modelEnv(), GRBwrite(model, target.c_str()));
}

void GRBModel::writeSol(const std::string &filename) const
{
	flush();
	GRB_CALL(modelEnv(), GRBwrite(model, filename.c_str()));
}

int GRBModel::status() const
{
	return getIntAttrOr("Status", GRB_LOADED);
}

/* Solve */
bool GRBModel::lpopt(char method, bool decrease_tol, bool initial)
{
	// Like CPXModel: never throw on solve failures, subproblems may well be infeasible.
	DOMINIQS_ASSERT(model);
	GRBenv *menv = modelEnv();
	int oldMethod = -1, oldPresolve = -1, oldCrossover = -1;
	GRB_CALL(menv, GRBgetintparam(menv, "Method", &oldMethod));
	GRB_CALL(menv, GRBgetintparam(menv, "Presolve", &oldPresolve));
	GRB_CALL(menv, GRBgetintparam(menv, "Crossover", &oldCrossover));

	switch (method)
	{
	case 'S': // automatic
		GRB_CALL(menv, GRBsetintparam(menv, "Method", -1));
		break;
	case 'P': // primal simplex
		GRB_CALL(menv, GRBsetintparam(menv, "Method", 0));
		break;
	case 'D': // dual simplex
		GRB_CALL(menv, GRBsetintparam(menv, "Method", 1));
		break;
	case 'B': // barrier (with crossover), presolve off as CPXModel does
		GRB_CALL(menv, GRBsetintparam(menv, "Method", 2));
		GRB_CALL(menv, GRBsetintparam(menv, "Presolve", 0));
		break;
	case 'A': // interior barrier point, no crossover (counterpart of CPX_NONBASIC_SOLN)
		GRB_CALL(menv, GRBsetintparam(menv, "Method", 2));
		GRB_CALL(menv, GRBsetintparam(menv, "Presolve", 0));
		GRB_CALL(menv, GRBsetintparam(menv, "Crossover", 0));
		break;
	default:
		throw std::runtime_error("Unexpected method for lpopt");
	}

	int error = GRBoptimize(model);
	dirty = false;

	// restoring params does not discard the solution in Gurobi
	GRB_CALL(menv, GRBsetintparam(menv, "Method", oldMethod));
	GRB_CALL(menv, GRBsetintparam(menv, "Presolve", oldPresolve));
	GRB_CALL(menv, GRBsetintparam(menv, "Crossover", oldCrossover));

	return error == 0;
}

bool GRBModel::isInfeasibleOrTimeReached()
{
	int st = status();
	if (st == GRB_INFEASIBLE || st == GRB_INF_OR_UNBD)
		return true;
	if (st == GRB_TIME_LIMIT || st == GRB_WORK_LIMIT)
	{
		// CPLEX: an LP stopped by time counts; a MIP stopped by time counts only without a solution.
		if (getIntAttrOr("IsMIP", 0) == 0)
			return true;
		return getIntAttrOr("SolCount", 0) == 0;
	}
	return false;
}

/* Proven infeasible, as opposed to stopped by a limit; see MIPModelI::isProvenInfeasible. */
bool GRBModel::isProvenInfeasible() const
{
	int st = status();
	return st == GRB_INFEASIBLE || st == GRB_INF_OR_UNBD;
}

bool GRBModel::mipopt()
{
	DOMINIQS_ASSERT(model);
	GRBoptimize(model); // silent, as in CPXModel
	dirty = false;
	return true;
}

/* Presolve/postsolve: not supported yet (v1), see header. */
bool GRBModel::presolve()
{
	static bool warned = false;
	if (!warned)
	{
		consoleWarn("Gurobi: MIP presolve is not supported yet (no uncrush in the C API); running on the original model");
		warned = true;
	}
	return true;
}

void GRBModel::postsolve()
{
	// no-op
}

std::vector<double> GRBModel::postsolveSolution(const std::vector<double> &preX) const
{
	return preX; // no presolve, solution spaces coincide
}

std::vector<double> GRBModel::presolveSolution(const std::vector<double> &origX) const
{
	return origX;
}

/* Get solution */
double GRBModel::objval() const
{
	flush();
	double ret = 0.0;
	GRB_CALL(modelEnv(), GRBgetdblattr(model, "ObjVal", &ret));
	return ret;
}

void GRBModel::sol(double *x, int first, int last) const
{
	getDblArray("X", x, first, last, ncols(), false);
}

void GRBModel::reduced_costs(double *x, int first, int last) const
{
	getDblArray("RC", x, first, last, ncols(), false);
}

bool GRBModel::isPrimalFeas() const
{
	if (getIntAttrOr("SolCount", 0) <= 0)
		return false;
	int st = status();
	if (st == GRB_OPTIMAL || st == GRB_SUBOPTIMAL)
		return true;
	// stopped early with a point available: check its violations against the tolerance
	double tol = dblParam(DblParam::FeasibilityTolerance);
	double constrVio = getDblAttrOr("ConstrVio", GRB_INFINITY);
	double boundVio = getDblAttrOr("BoundVio", GRB_INFINITY);
	return (constrVio <= tol) && (boundVio <= tol);
}

/* Parameters */
void GRBModel::handleCtrlC(bool flag)
{
	if (flag)
	{
		GRBModel_UserBreak = 0;
		GRBModel_ActiveModel = model;
		previousHandler = ::signal(SIGINT, userSignalBreak);
		restoreSignalHandler = true;
	}
	else if (restoreSignalHandler)
	{
		::signal(SIGINT, previousHandler);
		restoreSignalHandler = false;
		if (GRBModel_ActiveModel == model)
			GRBModel_ActiveModel = nullptr;
	}
}

bool GRBModel::aborted() const
{
	return GRBModel_UserBreak;
}

void GRBModel::seed(int seed)
{
	setIntParamEverywhere("Seed", seed);
}

void GRBModel::logging(bool log)
{
	setIntParamEverywhere("OutputFlag", log ? 1 : 0);
}

int GRBModel::intParam(IntParam which) const
{
	GRBenv *menv = modelEnv();
	int ivalue = 0;
	double dvalue = 0.0;

	switch (which)
	{
	case IntParam::Threads:
		GRB_CALL(menv, GRBgetintparam(menv, "Threads", &ivalue));
		return ivalue;
	case IntParam::SolutionLimit:
		GRB_CALL(menv, GRBgetintparam(menv, "SolutionLimit", &ivalue));
		return ivalue;
	case IntParam::NodeLimit: // double in Gurobi
		GRB_CALL(menv, GRBgetdblparam(menv, "NodeLimit", &dvalue));
		return dvalue >= static_cast<double>(INT_MAX) ? INT_MAX : static_cast<int>(dvalue);
	case IntParam::IterLimit: // double in Gurobi
		GRB_CALL(menv, GRBgetdblparam(menv, "IterationLimit", &dvalue));
		return dvalue >= static_cast<double>(INT_MAX) ? INT_MAX : static_cast<int>(dvalue);
	case IntParam::PdlpWarmStart:
		return 0;
	case IntParam::Presolve: // report CPLEX-style on/off
		GRB_CALL(menv, GRBgetintparam(menv, "Presolve", &ivalue));
		return ivalue == 0 ? 0 : 1;
	case IntParam::FeasOptMode:
		return feasOptMode;
	case IntParam::Emphasis:
		GRB_CALL(menv, GRBgetintparam(menv, "MIPFocus", &ivalue));
		return ivalue;
	default:
		throw std::runtime_error("Unknown integer parameter");
	}
}

void GRBModel::intParam(IntParam which, int value)
{
	switch (which)
	{
	case IntParam::Threads:
		setIntParamEverywhere("Threads", value);
		break;
	case IntParam::SolutionLimit:
		setIntParamEverywhere("SolutionLimit", value);
		break;
	case IntParam::NodeLimit:
		setDblParamEverywhere("NodeLimit", value >= INT_MAX ? GRB_INFINITY : static_cast<double>(value));
		break;
	case IntParam::IterLimit:
		setDblParamEverywhere("IterationLimit", value >= INT_MAX ? GRB_INFINITY : static_cast<double>(value));
		break;
	case IntParam::PdlpWarmStart:
		break;
	case IntParam::Presolve: // CPLEX on/off -> Gurobi off/automatic
		setIntParamEverywhere("Presolve", value == 0 ? 0 : -1);
		break;
	case IntParam::FeasOptMode:
		feasOptMode = value;
		break;
	case IntParam::Emphasis: // 1 = feasibility emphasis in CPLEX
		setIntParamEverywhere("MIPFocus", value);
		break;
	default:
		throw std::runtime_error("Unknown integer parameter");
	}
}

double GRBModel::dblParam(DblParam which) const
{
	GRBenv *menv = modelEnv();
	double value = 0.0;

	switch (which)
	{
	case DblParam::TimeLimit:
		GRB_CALL(menv, GRBgetdblparam(menv, "TimeLimit", &value));
		return fromGrb(value);
	case DblParam::FeasibilityTolerance:
		GRB_CALL(menv, GRBgetdblparam(menv, "FeasibilityTol", &value));
		return value;
	case DblParam::IntegralityTolerance:
		GRB_CALL(menv, GRBgetdblparam(menv, "IntFeasTol", &value));
		return value;
	case DblParam::WorkMem: // CPLEX MB, Gurobi GB
		GRB_CALL(menv, GRBgetdblparam(menv, "SoftMemLimit", &value));
		return value >= GRB_INFINITY ? INFBOUND : value * 1024.0;
	default:
		throw std::runtime_error("Unknown double parameter");
	}
}

void GRBModel::dblParam(DblParam which, double value)
{
	switch (which)
	{
	case DblParam::TimeLimit:
		setDblParamEverywhere("TimeLimit", toGrb(value));
		break;
	case DblParam::FeasibilityTolerance:
		setDblParamEverywhere("FeasibilityTol", value);
		break;
	case DblParam::IntegralityTolerance:
		setDblParamEverywhere("IntFeasTol", value);
		break;
	case DblParam::WorkMem:
		setDblParamEverywhere("SoftMemLimit", value >= INFBOUND ? GRB_INFINITY : value / 1024.0);
		break;
	default:
		throw std::runtime_error("Unknown double parameter");
	}
}

int GRBModel::intAttr(IntAttr which) const
{
	switch (which)
	{
	case IntAttr::Nodes:
		return static_cast<int>(getDblAttrOr("NodeCount", 0.0));
	case IntAttr::NodesLeft: // only available inside callbacks in Gurobi
		return 0;
	case IntAttr::BarrierIterations:
		return getIntAttrOr("BarIterCount", 0);
	case IntAttr::SimplexIterations:
		return static_cast<int>(getDblAttrOr("IterCount", 0.0));
	case IntAttr::PDLPIterations:
		return 0;
	default:
		throw std::runtime_error("Unknown integer attribute");
	}
}

double GRBModel::dblAttr(DblAttr which) const
{
	switch (which)
	{
	case DblAttr::MIPDualBound:
	{
		flush();
		double value = 0.0;
		GRB_CALL(modelEnv(), GRBgetdblattr(model, "ObjBound", &value));
		return fromGrb(value);
	}
	default:
		throw std::runtime_error("Unknown double attribute");
	}
}

void GRBModel::terminationReason(std::string &reason)
{
	reason = "-";
}

/* Access model data */
int GRBModel::nrows() const
{
	flush();
	int value = 0;
	GRB_CALL(modelEnv(), GRBgetintattr(model, "NumConstrs", &value));
	return value;
}

int GRBModel::ncols() const
{
	flush();
	int value = 0;
	GRB_CALL(modelEnv(), GRBgetintattr(model, "NumVars", &value));
	return value;
}

int GRBModel::nnz() const
{
	flush();
	int value = 0;
	GRB_CALL(modelEnv(), GRBgetintattr(model, "NumNZs", &value));
	return value;
}

double GRBModel::objOffset() const
{
	flush();
	double value = 0.0;
	GRB_CALL(modelEnv(), GRBgetdblattr(model, "ObjCon", &value));
	return value;
}

ObjSense GRBModel::objSense() const
{
	flush();
	int value = GRB_MINIMIZE;
	GRB_CALL(modelEnv(), GRBgetintattr(model, "ModelSense", &value));
	return (value > 0) ? ObjSense::MIN : ObjSense::MAX;
}

void GRBModel::lbs(double *lb, int first, int last) const
{
	getDblArray("LB", lb, first, last, ncols(), true);
}

void GRBModel::ubs(double *ub, int first, int last) const
{
	getDblArray("UB", ub, first, last, ncols(), true);
}

void GRBModel::objcoefs(double *obj, int first, int last) const
{
	getDblArray("Obj", obj, first, last, ncols(), false);
}

void GRBModel::ctypes(char *ctype, int first, int last) const
{
	flush();
	int count = checkRange(first, last, ncols());
	GRB_CALL(modelEnv(), GRBgetcharattrarray(model, "VType", first, count, ctype));
}

void GRBModel::sense(char *sense, int first, int last) const
{
	flush();
	int count = checkRange(first, last, nrows());
	GRB_CALL(modelEnv(), GRBgetcharattrarray(model, "Sense", first, count, sense));
	for (int i = 0; i < count; ++i)
		sense[i] = senseFromGrb(sense[i]);
}

void GRBModel::rhs(double *rhs, int first, int last) const
{
	getDblArray("RHS", rhs, first, last, nrows(), true);
}

void GRBModel::range(double *range, int first, int last) const
{
	int count = checkRange(first, last, nrows());
	std::fill(range, range + count, 0.0); // Gurobi has no ranged rows (they become a row + slack column)
}

void GRBModel::row(int ridx, dominiqs::SparseVector &row, char &sense, double &rhs, double &rngval) const
{
	DOMINIQS_ASSERT((ridx >= 0) && (ridx < nrows()));
	getVector(true, ridx, row);
	char grbSense = 0;
	GRB_CALL(modelEnv(), GRBgetcharattrarray(model, "Sense", ridx, 1, &grbSense));
	sense = senseFromGrb(grbSense);
	GRB_CALL(modelEnv(), GRBgetdblattrarray(model, "RHS", ridx, 1, &rhs));
	rhs = fromGrb(rhs);
	rngval = 0.0;
}

void GRBModel::rows(dominiqs::SparseMatrix &matrix) const
{
	getMatrix(true, matrix);
}

void GRBModel::col(int cidx, dominiqs::SparseVector &col, char &type, double &lb, double &ub, double &obj) const
{
	DOMINIQS_ASSERT((cidx >= 0) && (cidx < ncols()));
	getVector(false, cidx, col);
	GRB_CALL(modelEnv(), GRBgetdblattrarray(model, "LB", cidx, 1, &lb));
	GRB_CALL(modelEnv(), GRBgetdblattrarray(model, "UB", cidx, 1, &ub));
	GRB_CALL(modelEnv(), GRBgetdblattrarray(model, "Obj", cidx, 1, &obj));
	GRB_CALL(modelEnv(), GRBgetcharattrarray(model, "VType", cidx, 1, &type));
	lb = fromGrb(lb);
	ub = fromGrb(ub);
}

void GRBModel::cols(dominiqs::SparseMatrix &matrix) const
{
	getMatrix(false, matrix);
}

void GRBModel::colNames(std::vector<std::string> &names, int first, int last) const
{
	getNames("VarName", names, first, last, ncols());
}

void GRBModel::rowNames(std::vector<std::string> &names, int first, int last) const
{
	getNames("ConstrName", names, first, last, nrows());
}

/* Data modifications */
void GRBModel::addEmptyCol(const std::string &name, char ctype, double lb, double ub, double obj)
{
	GRB_CALL(modelEnv(), GRBaddvar(model, 0, nullptr, nullptr, obj, toGrb(lb), toGrb(ub), ctype,
								   name.empty() ? nullptr : name.c_str()));
	markDirty();
}

void GRBModel::addCol(const std::string &name, const int *idx, const double *val, int cnt, char ctype, double lb, double ub, double obj)
{
	std::vector<int> ind(idx, idx + (cnt > 0 ? cnt : 0));
	std::vector<double> coef(val, val + (cnt > 0 ? cnt : 0));
	GRB_CALL(modelEnv(), GRBaddvar(model, cnt, cnt > 0 ? ind.data() : nullptr, cnt > 0 ? coef.data() : nullptr,
								   obj, toGrb(lb), toGrb(ub), ctype, name.empty() ? nullptr : name.c_str()));
	markDirty();
}

void GRBModel::addRow(const std::string &name, const int *idx, const double *val, int cnt, char sense, double rhs, double rngval)
{
	// A ranged row would add a slack column in Gurobi and shift column indices; no caller needs it.
	if (sense == 'R')
		throw std::runtime_error("GRBModel::addRow: ranged rows are not supported");
	std::vector<int> ind(idx, idx + cnt);
	std::vector<double> coef(val, val + cnt);
	GRB_CALL(modelEnv(), GRBaddconstr(model, cnt, cnt > 0 ? ind.data() : nullptr, cnt > 0 ? coef.data() : nullptr,
									  senseToGrb(sense), toGrb(rhs), name.empty() ? nullptr : name.c_str()));
	markDirty();
}

void GRBModel::delRow(int ridx)
{
	delRows(ridx, ridx);
}

void GRBModel::delCol(int cidx)
{
	delCols(cidx, cidx);
}

void GRBModel::delRows(int first, int last)
{
	int count = checkRange(first, last, nrows());
	std::vector<int> ind(count);
	for (int i = 0; i < count; ++i)
		ind[i] = first + i;
	GRB_CALL(modelEnv(), GRBdelconstrs(model, count, ind.data()));
	markDirty();
}

void GRBModel::delCols(int first, int last)
{
	int count = checkRange(first, last, ncols());
	std::vector<int> ind(count);
	for (int i = 0; i < count; ++i)
		ind[i] = first + i;
	GRB_CALL(modelEnv(), GRBdelvars(model, count, ind.data()));
	markDirty();
}

void GRBModel::objSense(ObjSense objsen)
{
	GRB_CALL(modelEnv(), GRBsetintattr(model, "ModelSense", static_cast<int>(objsen)));
	markDirty();
}

void GRBModel::objOffset(double val)
{
	GRB_CALL(modelEnv(), GRBsetdblattr(model, "ObjCon", val));
	markDirty();
}

void GRBModel::lb(int cidx, double val)
{
	lbs(1, &cidx, &val);
}

void GRBModel::lbs(int cnt, const int *cols, const double *values)
{
	setDblList("LB", cnt, cols, values, true);
}

void GRBModel::ub(int cidx, double val)
{
	ubs(1, &cidx, &val);
}

void GRBModel::ubs(int cnt, const int *cols, const double *values)
{
	setDblList("UB", cnt, cols, values, true);
}

void GRBModel::fixCol(int cidx, double val)
{
	lb(cidx, val);
	ub(cidx, val);
}

void GRBModel::objcoef(int cidx, double val)
{
	objcoefs(1, &cidx, &val);
}

void GRBModel::objcoefs(int cnt, const int *cols, const double *values)
{
	setDblList("Obj", cnt, cols, values, false);
}

void GRBModel::ctype(int cidx, char val)
{
	DOMINIQS_ASSERT((val == 'B') || (val == 'I') || (val == 'C'));
	ctypes(1, &cidx, &val);
}

void GRBModel::ctypes(int cnt, const int *cols, const char *values)
{
	if (cnt <= 0)
		return;
	std::vector<int> ind(cols, cols + cnt);
	std::vector<char> val(values, values + cnt);
	GRB_CALL(modelEnv(), GRBsetcharattrlist(model, "VType", cnt, ind.data(), val.data()));
	markDirty();
}

void GRBModel::switchToLP()
{
	int n = ncols();
	if (n == 0)
		return;
	std::vector<char> continuous(n, GRB_CONTINUOUS);
	GRB_CALL(modelEnv(), GRBsetcharattrarray(model, "VType", 0, n, continuous.data()));
	markDirty();
}

void GRBModel::switchToMIP()
{
	// no-op, exactly like CPXModel (feaspump.cpp relies on this behaviour)
}

/* Private interface */
GRBModel *GRBModel::clone_impl() const
{
	flush();
	GRBmodel *copy = GRBcopymodel(model);
	if (!copy)
		throwGurobiError(modelEnv(), -1, "GRBcopymodel");
	return new GRBModel(env, copy);
}

GRBModel *GRBModel::presolvedmodel_impl() const
{
	return nullptr; // v1: no presolve, callers treat this as "no reductions"
}

void GRBModel::updateModelVarBounds(std::optional<boost::dynamic_bitset<>> vars_entering_problem, std::optional<boost::dynamic_bitset<>> vars_leaving_problem)
{
	std::vector<int> ind;
	std::vector<double> val;
	if (vars_entering_problem)
		for (auto i = vars_entering_problem->find_first(); i != boost::dynamic_bitset<>::npos; i = vars_entering_problem->find_next(i))
		{
			ind.push_back(static_cast<int>(i));
			val.push_back(1.0);
		}
	if (vars_leaving_problem)
		for (auto i = vars_leaving_problem->find_first(); i != boost::dynamic_bitset<>::npos; i = vars_leaving_problem->find_next(i))
		{
			ind.push_back(static_cast<int>(i));
			val.push_back(0.0);
		}
	if (!ind.empty())
		ubs(static_cast<int>(ind.size()), ind.data(), val.data());
}

void GRBModel::findSetOfConflictingVariables(boost::dynamic_bitset<> inactive_binary_vars, std::vector<int> &conflicting_constraints, std::vector<int> &conflicting_vars, bool optimize_set, double time_left)
{
	conflicting_vars.clear();
	conflicting_constraints.clear();
	int numCols = ncols();
	int numRows = nrows();
	dblParam(DblParam::TimeLimit, time_left); // as CPXModel: the time limit stays set on this model

	// GRBfeasrelax and GRBcomputeIIS work on a copy: feasrelax adds columns/rows to its model,
	// while CPXfeasopt leaves the problem unchanged.
	flush();
	GRBmodel *work = GRBcopymodel(model);
	if (!work)
		throwGurobiError(modelEnv(), -1, "GRBcopymodel");
	struct Guard
	{
		GRBmodel *m;
		~Guard() { GRBfreemodel(m); }
	} guard{work};

	consoleLog(" * Attempt to make kernel LP feasible");
	if (optimize_set)
	{
		// CPLEX FeasOpt mode -> (relaxobjtype, minrelax): MIN/OPT x SUM(0)/INF(2)/QUAD(1)
		static const int relaxObjType[3] = {0, 2, 1};
		int relaxobjtype = relaxObjType[(feasOptMode / 2) % 3];
		int minrelax = feasOptMode % 2;

		std::vector<double> ubpen(numCols, GRB_INFINITY); // infinity = bound may not be relaxed
		for (auto i = inactive_binary_vars.find_first(); i != boost::dynamic_bitset<>::npos; i = inactive_binary_vars.find_next(i))
			ubpen[i] = 1.0;
		consoleLog(" * {} inactive binary vars considered for entering kernel (i.e., being activated: ub -> 1)", inactive_binary_vars.count());

		double feasobj = 0.0;
		if (GRBfeasrelax(work, relaxobjtype, minrelax, nullptr, ubpen.data(), nullptr, &feasobj) != 0)
			return;
		if (GRBoptimize(work) != 0)
			return;
		int st = 0, solCount = 0, totalCols = 0;
		GRBgetintattr(work, "Status", &st);
		GRBgetintattr(work, "SolCount", &solCount);
		GRBgetintattr(work, "NumVars", &totalCols);
		if (solCount <= 0 || totalCols <= numCols)
			return;

		// The artificial column relaxing the upper bound of variable v is named "ArtU_<name of v>".
		std::vector<std::string> names;
		colNames(names);
		std::unordered_map<std::string, int> indexOf;
		indexOf.reserve(names.size());
		for (int j = 0; j < numCols; ++j)
			indexOf.emplace(names[j], j);

		int numArt = totalCols - numCols;
		std::vector<char *> artNames(numArt, nullptr);
		std::vector<double> artValues(numArt, 0.0);
		GRB_CALL(GRBgetenv(work), GRBgetstrattrarray(work, "VarName", numCols, numArt, artNames.data()));
		GRB_CALL(GRBgetenv(work), GRBgetdblattrarray(work, "X", numCols, numArt, artValues.data()));
		for (int a = 0; a < numArt; ++a)
		{
			if (!artNames[a] || std::strncmp(artNames[a], "ArtU_", 5) != 0)
				continue;
			auto it = indexOf.find(artNames[a] + 5);
			if (it == indexOf.end())
				continue;
			int var_index = it->second;
			if (inactive_binary_vars[var_index] && greaterThan(artValues[a], 0.0))
				conflicting_vars.push_back(var_index);
		}
	}
	else
	{
		if (GRBcomputeIIS(work) != 0)
			return;
		std::vector<int> iisConstr(numRows, 0), iisUB(numCols, 0);
		if (numRows > 0 && GRBgetintattrarray(work, "IISConstr", 0, numRows, iisConstr.data()) == 0)
			for (int i = 0; i < numRows; ++i)
				if (iisConstr[i])
					conflicting_constraints.push_back(i);
		if (numCols > 0 && GRBgetintattrarray(work, "IISUB", 0, numCols, iisUB.data()) == 0)
			for (int j = 0; j < numCols; ++j)
				if (iisUB[j])
					conflicting_vars.push_back(j);
	}
}
