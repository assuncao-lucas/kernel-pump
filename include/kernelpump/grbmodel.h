/**
 * @file grbmodel.h
 * @brief Implementation of MIPModelI for Gurobi (C API)
 *
 * Mirrors CPXModel method for method. Design notes (see the Gurobi backend spec):
 *  - one GRBenv per process, shared by all models through a shared_ptr;
 *  - parameter setters write to the model's own env AND to the master env, so that
 *    parameters set before readModel() (as main.cpp does) survive, like CPLEX's env-wide params;
 *  - Gurobi applies model changes lazily, so getters flush pending changes first;
 *  - the code base uses INFBOUND (1e20) for infinity, Gurobi uses GRB_INFINITY (1e100):
 *    values are translated at this class boundary;
 *  - MIP presolve is not supported yet (Gurobi publishes no crush/uncrush), so presolve()
 *    reports "no reductions" and FP/KP run on the original model.
 */

#ifndef GRBMODEL_H
#define GRBMODEL_H

#include "mipmodel.h"
#include <gurobi_c.h>
#include <memory>
#include <string>
#include <vector>

using namespace dominiqs;

class GRBModel : public MIPModelI
{
public:
	/**
	 * Create the process-wide Gurobi environment (first call only) and an empty model.
	 *
	 * licenseFile: when non-empty, GRB_LICENSE_FILE is set to it for this process before the
	 * environment starts, overriding any user-level value.
	 * requireLicenseFile: when true and licenseFile is empty, throw before touching Gurobi.
	 * This guards against silently checking out whatever license the user environment points to.
	 */
	explicit GRBModel(const std::string &licenseFile = "", bool requireLicenseFile = true);
	GRBModel(std::shared_ptr<GRBenv> env, GRBmodel *model);
	~GRBModel() override;
	GRBModel(const GRBModel &) = delete;
	GRBModel &operator=(const GRBModel &) = delete;
	std::unique_ptr<GRBModel> clone() const { return std::unique_ptr<GRBModel>(this->clone_impl()); }
	/* Read/Write */
	void readModel(const std::string &filename) override;
	void writeModel(const std::string &filename, const std::string &format = "") const override;
	void writeSol(const std::string &filename) const override;
	/* Solve */
	bool lpopt(char method, bool decrease_tol, bool initial) override;
	int status() const override;
	bool mipopt() override;
	/* Presolve/postsolve */
	bool presolve() override;
	void postsolve() override;
	std::vector<double> postsolveSolution(const std::vector<double> &preX) const override;
	std::vector<double> presolveSolution(const std::vector<double> &origX) const override;
	/* Get solution */
	double objval() const override;
	void sol(double *x, int first = 0, int last = -1) const override;
	void reduced_costs(double *x, int first = 0, int last = -1) const override;
	bool isPrimalFeas() const override;
	/* Parameters */
	void handleCtrlC(bool flag) override;
	bool aborted() const override;
	void seed(int seed) override;
	void logging(bool log) override;
	int intParam(IntParam which) const override;
	void intParam(IntParam which, int value) override;
	double dblParam(DblParam which) const override;
	void dblParam(DblParam which, double value) override;
	int intAttr(IntAttr which) const override;
	double dblAttr(DblAttr which) const override;
	void terminationReason(std::string &reason) override;
	/* Access model data */
	int nrows() const override;
	int ncols() const override;
	int nnz() const override;
	double objOffset() const override;
	ObjSense objSense() const override;
	void lbs(double *lb, int first = 0, int last = -1) const override;
	void ubs(double *ub, int first = 0, int last = -1) const override;
	void objcoefs(double *obj, int first = 0, int last = -1) const override;
	void ctypes(char *ctype, int first = 0, int last = -1) const override;
	void sense(char *sense, int first = 0, int last = -1) const override;
	void rhs(double *rhs, int first = 0, int last = -1) const override;
	void range(double *range, int first = 0, int last = -1) const override;
	void row(int ridx, dominiqs::SparseVector &row, char &sense, double &rhs, double &rngval) const override;
	void rows(dominiqs::SparseMatrix &matrix) const override;
	void col(int cidx, dominiqs::SparseVector &col, char &type, double &lb, double &ub, double &obj) const override;
	void cols(dominiqs::SparseMatrix &matrix) const override;
	void colNames(std::vector<std::string> &names, int first = 0, int last = -1) const override;
	void rowNames(std::vector<std::string> &names, int first = 0, int last = -1) const override;
	/* Data modifications */
	void addEmptyCol(const std::string &name, char ctype, double lb, double ub, double obj) override;
	void addCol(const std::string &name, const int *idx, const double *val, int cnt, char ctype, double lb, double ub, double obj) override;
	void addRow(const std::string &name, const int *idx, const double *val, int cnt, char sense, double rhs, double rngval = 0.0) override;
	void delRow(int ridx) override;
	void delCol(int cidx) override;
	void delRows(int first, int last) override;
	void delCols(int first, int last) override;
	void objSense(ObjSense objsen) override;
	void objOffset(double val) override;
	void lb(int cidx, double val) override;
	void lbs(int cnt, const int *cols, const double *values) override;
	void ub(int cidx, double val) override;
	void ubs(int cnt, const int *cols, const double *values) override;
	void fixCol(int cidx, double val) override;
	void objcoef(int cidx, double val) override;
	void objcoefs(int cnt, const int *cols, const double *values) override;
	void ctype(int cidx, char val) override;
	void ctypes(int cnt, const int *cols, const char *values) override;
	void switchToLP() override;
	void switchToMIP() override;
	void updateModelVarBounds(std::optional<boost::dynamic_bitset<>> vars_entering_problem, std::optional<boost::dynamic_bitset<>> vars_leaving_problem) override;
	void findSetOfConflictingVariables(boost::dynamic_bitset<> inactive_binary_vars, std::vector<int> &conflicting_constraints, std::vector<int> &conflicting_vars, bool optimize_set, double time_limit) override;
	bool isInfeasibleOrTimeReached() override;
	/* Access to underlying Gurobi objects */
	GRBenv *getMasterEnv() const { return env.get(); }
	GRBmodel *getModel() const { return model; }

private:
	GRBModel *clone_impl() const override;
	GRBModel *presolvedmodel_impl() const override;

	/* helpers */
	GRBenv *modelEnv() const { return GRBgetenv(model); }
	void flush() const;	   // apply pending lazy updates before a query
	void markDirty() { dirty = true; }
	int checkRange(int &first, int &last, int size) const; // validates and normalizes [first,last], returns count
	void setIntParamEverywhere(const char *name, int value);
	void setDblParamEverywhere(const char *name, double value);
	int getIntAttrOr(const char *name, int fallback) const;
	double getDblAttrOr(const char *name, double fallback) const;
	void getDblArray(const char *attr, double *out, int first, int last, int size, bool mapInf) const;
	void setDblList(const char *attr, int cnt, const int *cols, const double *values, bool mapInf);
	void getNames(const char *attr, std::vector<std::string> &names, int first, int last, int size) const;
	void getMatrix(bool byRow, dominiqs::SparseMatrix &matrix) const;
	void getVector(bool byRow, int index, dominiqs::SparseVector &vec) const;

private:
	std::shared_ptr<GRBenv> env; // master environment (shared by all models of the process)
	GRBmodel *model = nullptr;	 // owned
	mutable bool dirty = false;
	int feasOptMode = 0; // CPLEX FeasOpt mode numbering (0 = MIN_SUM), mapped onto GRBfeasrelax arguments
	bool restoreSignalHandler = false;
	using SignalHandler = void (*)(int);
	SignalHandler previousHandler = nullptr;
};

#endif /* GRBMODEL_H */
