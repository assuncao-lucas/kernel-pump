/**
 * @file mipmodel.h
 * @brief MIP Model Interface
 *
 * @author Domenico Salvagnin <dominiqs at gmail dot com>
 * @author Gioni Mexi <gionimexi at gmail dot com>
 * 2023
 */

#ifndef MIPMODEL_H
#define MIPMODEL_H

#include <string>
#include <vector>
#include <unordered_set>
#include <memory>
#include <utils/maths.h>
#include <utils/asserter.h>
#include <utils/consolelog.h>
#include <boost/dynamic_bitset.hpp>

using namespace dominiqs;

/* Objective Sense values */
enum class ObjSense
{
	MIN = 1,
	MAX = -1
};

enum class IntParam
{
	Threads,
	SolutionLimit,
	NodeLimit,
	IterLimit,
	PdlpWarmStart,
	Presolve,
	FeasOptMode,
	Emphasis
};

enum class DblParam
{
	TimeLimit,
	FeasibilityTolerance,
	IntegralityTolerance,
	PdlpTolerance,
	PdlpToleranceDecreaseFactor,
	WorkMem
};

enum class IntAttr
{
	Nodes,
	NodesLeft,
	BarrierIterations,
	SimplexIterations,
	PDLPIterations
};

enum class DblAttr
{
	MIPDualBound
};

/* Interface for a MIP model and solver */
class MIPModelI
{
public:
	virtual ~MIPModelI() {}
	std::unique_ptr<MIPModelI> clone() const { return std::unique_ptr<MIPModelI>(this->clone_impl()); }
	/* Read/Write */
	virtual void readModel(const std::string &filename) = 0;
	virtual void writeModel(const std::string &filename, const std::string &format = "") const = 0;
	virtual void writeSol(const std::string &filename) const = 0;
	/* Solve */
	virtual bool lpopt(char method, bool decrease_tol, bool initial) = 0;
	virtual int status() const = 0;
	virtual bool mipopt() = 0;
	/* Presolve/Postsolve */
	virtual bool presolve() = 0;
	virtual void postsolve() = 0;
	std::unique_ptr<MIPModelI> presolvedModel() const { return std::unique_ptr<MIPModelI>(this->presolvedmodel_impl()); }
	virtual std::vector<double> postsolveSolution(const std::vector<double> &preX) const = 0; /* get solution vector in the original space */
	virtual std::vector<double> presolveSolution(const std::vector<double> &origX) const = 0; /* get solution vector in the presolved space */
	/* Get solution */
	virtual double objval() const = 0;
	virtual void sol(double *x, int first = 0, int last = -1) const = 0;
	virtual void reduced_costs(double *x, int first = 0, int last = -1) const = 0;
	virtual bool isPrimalFeas() const = 0;
	/* Parameters */
	virtual void handleCtrlC(bool flag) = 0;
	virtual bool aborted() const = 0;
	virtual void seed(int seed) = 0;
	virtual void logging(bool log) = 0;
	virtual int intParam(IntParam which) const = 0;
	virtual void intParam(IntParam which, int value) = 0;
	virtual double dblParam(DblParam which) const = 0;
	virtual void dblParam(DblParam which, double value) = 0;
	virtual int intAttr(IntAttr which) const = 0;
	virtual double dblAttr(DblAttr which) const = 0;
	virtual void terminationReason(std::string &reason) = 0;
	/* Access model data */
	virtual int nrows() const = 0;
	virtual int ncols() const = 0;
	virtual int nnz() const = 0;
	virtual double objOffset() const = 0;
	virtual ObjSense objSense() const = 0;
	virtual void lbs(double *lb, int first = 0, int last = -1) const = 0;
	virtual void ubs(double *ub, int first = 0, int last = -1) const = 0;
	virtual void objcoefs(double *obj, int first = 0, int last = -1) const = 0;
	virtual void ctypes(char *ctype, int first = 0, int last = -1) const = 0;
	virtual void sense(char *sense, int first = 0, int last = -1) const = 0;
	virtual void range(double *range, int first = 0, int last = -1) const = 0;
	virtual void rhs(double *rhs, int first = 0, int last = -1) const = 0;
	virtual void row(int ridx, dominiqs::SparseVector &row, char &sense, double &rhs, double &rngval) const = 0;
	virtual void rows(dominiqs::SparseMatrix &matrix) const = 0;
	std::tuple<double, int> computeIntegralityGap(const std::vector<double> &x, double integralityEps) const
	{
		int num_infeas_values = 0; // number of integer/binary variables outside their domain (in practice, should be the ones that assume fractional values).
		// but, for correctness, we also include the ones that assume integer values outside their given domain. For instance, x = 1, when x mist be within {3,4,5}, for instance.
		bool empty = x.empty();
		double integrality_gap = 0.0;
		int num_vars = ncols();
		// consoleInfo("{}", num_vars);
		std::vector<double> lbs(num_vars), ubs(num_vars);
		std::vector<char> x_type(num_vars);
		this->ctypes(&(x_type[0]));

		if (!empty)
		{
			this->lbs(&(lbs[0]));
			this->ubs(&(ubs[0]));
		}

		int num_int_vars = 0;
		double min_gap_var = 0.0, dist_lb = 0.0, dist_ub = 0.0;
		for (int i = 0; i < num_vars; ++i)
		{
			// consoleInfo("{}", x_type[i]);
			if (x_type[i] == 'B' || x_type[i] == 'I')
			{
				++num_int_vars;

				if (empty) // if empty solution, consider all integer/binary variables with gap = 1.0;
					min_gap_var = 1.0;
				else
				{

					if (greaterEqualThan(x[i], lbs[i], integralityEps) && lessEqualThan(x[i], ubs[i], integralityEps))
					{
						min_gap_var = std::fabs(x[i] - round(x[i]));
						// consoleInfo("frac {} | [{},{}] | {}", x[i], lbs[i], ubs[i], min_gap_var);
					}
					else
					{
						dist_lb = std::fabs(x[i] - lbs[i]);
						dist_ub = std::fabs(x[i] - ubs[i]);
						min_gap_var = std::min(dist_lb, dist_ub);
						// consoleInfo("frac {} | [{},{}] | min [{},{}] = {}", x[i], lbs[i], ubs[i], dist_lb, dist_ub, min_gap_var);
					}
				}

				if (greaterThan(min_gap_var, 0.0, integralityEps))
					++num_infeas_values;
				integrality_gap += min_gap_var;
			}
		}

		if (num_int_vars == 0)
			return {0.0, 0};

		// std::cout << "relative_integrality_gap: " << integrality_gap << "/" << num_int_vars << "= " << integrality_gap / num_int_vars << std::endl;

		return {integrality_gap, num_infeas_values};
	}

	std::shared_ptr<std::vector<ConstraintPtr>> rows()
	{
		if (!constraints)
			retrieveConstraints_impl();
		return constraints;
	}

	std::shared_ptr<std::vector<boost::dynamic_bitset<>>> colsDependency()
	{
		if (!dependency)
			retrieveDependency_impl();
		return dependency;
	}

	int numIntegerAndBinaryCols() const
	{
		int counter = 0;
		int num_vars = ncols();
		std::vector<char> types(num_vars);
		ctypes(&(types[0]));
		for (int i = 0; i < num_vars; ++i)
			if (types[i] == 'B' || types[i] == 'I')
				++counter;
		return counter;
	}

	int numBinaryCols() const
	{
		int counter = 0;
		int num_vars = ncols();
		std::vector<char> types(num_vars);
		ctypes(&(types[0]));
		for (int i = 0; i < num_vars; ++i)
			if (types[i] == 'B')
				++counter;
		return counter;
	}

	bool isSolutionFeasible(const std::vector<double> &x)
	{
		// std::vector<std::string> rNames(nrows(), "");
		// rowNames(rNames);
		// auto rows = this->rows();
		// for (int i = 0; i < rows->size(); i++)
		// {
		// 	const ConstraintPtr &c = (*rows)[i];
		// 	// model->row(i, c->row, c->sense, c->rhs, c->range);
		// 	if (c->sense == 'N')
		// 		continue;
		// 	if (!c->satisfiedBy(&x[0]))
		// 	{
		// 		if (!silent)
		// 		{
		// 			consoleError("Constraint {} violated by {}", rNames[i], c->violation(&x[0]));
		// 			consoleWarn("sol status: {}", status());
		// 			getchar();
		// 			getchar();
		// 			return false;
		// 		}
		// 	}
		// }

		const auto &rows = *(this->rows());
		for (auto c : rows)
		{
			if (!c->satisfiedBy(&x[0]))
				return false;
		}
		return true;
	}

	void printDependencies()
	{
		consoleInfo("[dependency matrix]");
		colsDependency();
		std::vector<std::string> xNames;
		colNames(xNames);
		for (int i = 0; i < dependency->size(); ++i)
		{
			std::string dep_list;
			for (int j = ((*dependency)[i]).find_first(); j != boost::dynamic_bitset<>::npos; j = ((*dependency)[i]).find_next(j))
			{
				dep_list += xNames[j] + " ";
			}
			consoleLog("col {}: {}", xNames[i], dep_list);
		}
	}
	virtual void findSetOfConflictingVariables(boost::dynamic_bitset<> inactive_binary_vars, std::vector<int> &conflicting_constraints, std::vector<int> &conflicting_vars, bool optimize_set, double time_limit) = 0;
	virtual void col(int cidx, dominiqs::SparseVector &col, char &type, double &lb, double &ub, double &obj) const = 0;
	virtual void cols(dominiqs::SparseMatrix &matrix) const = 0;
	virtual void colNames(std::vector<std::string> &names, int first = 0, int last = -1) const = 0;
	virtual void rowNames(std::vector<std::string> &names, int first = 0, int last = -1) const = 0;
	/* Data modifications */
	virtual void addEmptyCol(const std::string &name, char ctype, double lb, double ub, double obj) = 0;
	virtual void addCol(const std::string &name, const int *idx, const double *val, int cnt, char ctype, double lb, double ub, double obj) = 0;
	virtual void addRow(const std::string &name, const int *idx, const double *val, int cnt, char sense, double rhs, double rngval = 0.0) = 0;
	virtual void delRow(int ridx) = 0;
	virtual void delCol(int cidx) = 0;
	virtual void delRows(int first, int last) = 0;
	virtual void delCols(int first, int last) = 0;
	virtual void objSense(ObjSense objsen) = 0;
	virtual void objOffset(double val) = 0;
	virtual void lb(int cidx, double val) = 0;
	virtual void lbs(int cnt, const int *cols, const double *values) = 0;
	virtual void ub(int cidx, double val) = 0;
	virtual void ubs(int cnt, const int *cols, const double *values) = 0;
	virtual void fixCol(int cidx, double val) = 0;
	virtual void objcoef(int cidx, double val) = 0;
	virtual void objcoefs(int cnt, const int *cols, const double *values) = 0;
	virtual void ctype(int cidx, char val) = 0;
	virtual void ctypes(int cnt, const int *cols, const char *values) = 0;
	virtual void switchToLP() = 0;
	virtual void switchToMIP() = 0;
	virtual void updateModelVarBounds(std::optional<boost::dynamic_bitset<>> vars_entering_problem, std::optional<boost::dynamic_bitset<>> vars_leaving_problem) = 0;
	virtual bool isInfeasibleOrTimeReached() = 0;

private:
	virtual MIPModelI *clone_impl() const = 0;
	virtual MIPModelI *presolvedmodel_impl() const = 0;
	void retrieveDependency_impl()
	{
		int nCols = ncols();
		dependency = std::make_shared<std::vector<boost::dynamic_bitset<>>>(nCols, boost::dynamic_bitset<>(nCols, 0));
		const auto &rows = *(this->rows());
		// for (int rowIndex = 0; rowIndex < rows.size(); ++rowIndex)
		// {
		// 	(*dependency)[rowIndex] = std::list<int>();
		for (const auto &constraint : rows)
		{
			for (int i = 0; i < constraint->row.capacity(); ++i)
			{
				int var_i_index = constraint->row.idx()[i];
				for (int j = i + 1; j < constraint->row.capacity(); ++j)
				{
					int var_j_index = constraint->row.idx()[j];
					(*dependency)[var_i_index][var_j_index] = 1;
					(*dependency)[var_j_index][var_i_index] = 1;
				}
			}
		}
	}

	void retrieveConstraints_impl()
	{
		dominiqs::SparseMatrix matrix;
		rows(matrix);

		int beg = 0;
		int end = nrows();

		std::vector<char> sense(end, 0);
		std::vector<double> rhs(end, 0.0);
		std::vector<double> rngval(end, 0.0);

		// here to correctly handle empty constraints
		// get rhs
		this->rhs(&rhs[0], beg, end - 1);
		// get sense
		this->sense(&sense[0], beg, end - 1);
		// get rngval
		this->range(&rngval[0], beg, end - 1);

		constraints = std::make_unique<std::vector<ConstraintPtr>>(end);
		// std::cout << constraints->size() << std::endl;
		// for (int i = 0; i < constraints.size(); ++i)
		// 	std::cout << *(constraints)[i] << std::endl;
		// getchar();
		// getchar();
		for (int i = 0; i < end; ++i)
		{
			(*constraints)[i] = std::make_shared<Constraint>();
			// CPLEX treats ranged rows considering the constraint satisfied
			// if the linear expression is in the range [rhs, rhs+rngval].
			// However, we interpret ranged rows differently, and the allowed
			// range is [rhs-rngval,rhs] (in both cases, rngval >= 0).
			// So we might have to update the rhs
			if (sense[i] == 'R')
			{
				DOMINIQS_ASSERT(rngval[i] >= 0.0);
				(rhs[i]) += rngval[i];
			}

			// fill the actual Constraint objects with the previously filled structures.
			auto &curr_constraint = *((*constraints)[i]);
			curr_constraint.sense = sense[i];
			curr_constraint.rhs = rhs[i];
			curr_constraint.range = rngval[i];

			// compute number of non-zero entries in constraint coefficient vector.
			// if last row, the end of the sub-vector of its corresponding indexes/values is the end of the matind/matval (position matrix.nnz-1)
			int constraint_nnz = i < end - 1 ? matrix.matbeg[i + 1 - beg] - matrix.matbeg[i - beg] : matrix.nnz - matrix.matbeg[i - beg];
			// beg  ==0, and could be omitted.
			// just to be aware: when a given row is empty (has no non-zero coefficient), the corresponding matbeg (begin of its elements
			// in matval and matind corresponds to the same begin of the next - non-empty - constraint). So, it is safe to always consider, for a row i, the
			// number of non-zero elements as (matbeg[i+1] - matbeg[i]), even when i is empty! Case empty, matbeg[i+1] - matbeg[i] == 0,

			// std::cout << "constraint " << i << std::endl;
			// std::cout << "matbeg: " << matrix.matbeg[i] << std::endl;
			// if (constraint_nnz > 0)
			// 	std::cout << "range: " << "[" << matrix.matbeg[i] << "," << matrix.matbeg[i] + constraint_nnz - 1 << "]" << std::endl;
			// else
			// 	std::cout << "range: -" << std::endl;

			auto &row = curr_constraint.row;
			if (constraint_nnz > 0)
			{
				row.resize(constraint_nnz);
				for (int j = 0; j < constraint_nnz; ++j)
				{
					auto index = matrix.matbeg[i - beg] + j;
					row.idx()[j] = matrix.matind[index];
					row.coef()[j] = matrix.matval[index];
				}
			}
			else
				row.clear();
		}
	}

protected:
	std::shared_ptr<std::vector<ConstraintPtr>> constraints = nullptr;			// must be computed just once, at the first call of rows(void).
	std::shared_ptr<std::vector<boost::dynamic_bitset<>>> dependency = nullptr; // must be computed just once, at the first call of colsDependency(void).
};

using MIPModelPtr = std::shared_ptr<MIPModelI>;

#endif /* MIPMODEL_H */
