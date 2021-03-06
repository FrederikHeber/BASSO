/*
 * Solver.hpp
 *
 *  Created on: Jun 26, 2015
 *      Author: heber
 */

#ifndef MATRIXFACTORIZERBASE_SOLVERS_INRANGESOLVER_HPP_
#define MATRIXFACTORIZERBASE_SOLVERS_INRANGESOLVER_HPP_

#include "BassoConfig.h"

#include <Eigen/Dense>

#include "MatrixFactorizer/Database/InnerProblemDatabase.hpp"
#include "Minimizations/types.hpp"
#include "Solvers/AuxiliaryConstraints/AuxiliaryConstraints.hpp"

class CommandLineOptions;
class InnerProblemDatabase;

/** Solver combines projection onto range and inverse problem solving with
 * SESOP that requires right-hand side to be in range.
 *
 */
struct InRangeSolver
{
	/** Constructor for class InRangeSolver.
	 *
	 * @param _opts options defining the manner of solving.
	 * @param _overall_keys set of keys to accumulate from solver's overall
	 * 			table
	 * @param _projection_delta delta for projection the right-hand side (i.e.
	 * 		  the current column/row of the matrix factor that is modified)
	 * 		  onto the range of the operator (i.e. the respective fixed matrix
	 * 		   factor), e.g. 1e-8
	 */
	InRangeSolver(
			const CommandLineOptions &_opts,
			const InnerProblemDatabase::keys_t &_overall_keys,
			const double _projection_delta = 1e-8
			);

	/** Copy constructor for class InRangeSolver.
	 *
	 * We need this to be able to store solvers in a std::vector, e.g. for
	 * parallel runs when we want to preconstruct them.
	 *
	 * @param other other instance to copy from
	 */
	InRangeSolver(const InRangeSolver &other);

	/** Functor that projects onto range and solves.
	 *
	 * @param _matrix problem matrix
	 * @param _rhs right-hand side
	 * @param _solution_start initial value for solution
	 * @param _solution on return containing solution
	 * @param _dim unique index (for output statements only)
	 * @param _auxiliary_constraints auxiliary constraints to fulfill
	 * @return true - success
	 */
	bool operator()(
			const Eigen::MatrixXd &_matrix,
			const Eigen::VectorXd &_rhs,
			const Eigen::VectorXd &_solution_start,
			Eigen::VectorXd &_solution,
			const unsigned int _dim,
			const AuxiliaryConstraints::ptr_t&_auxiliary_constraints = AuxiliaryConstraints::ptr_t()
			);

	/** Insert accumulated values from projector problem into given \a _table.
	 *
	 * @param _table table to insert values to
	 * @param _suffix suffix to database column key for distinction
	 */
	void insertAccumulatedProjectorValues(
			Table &_table,
			const std::string &_suffix) const;

	/** Insert accumulated values from solver problem into given \a _table.
	 *
	 * @param _table table to insert values to
	 * @param _suffix suffix to database column key for distinction
	 */
	void insertAccumulatedSolverValues(
			Table &_table,
			const std::string &_suffix) const;

	/** Getter for the accumulated values of the projector inner problem.
	 *
	 * This is used by the parallel implementation to obtain the values
	 * for the the Slave nodes to transmit to Master node.
	 *
	 * @return ref to accumulated values from projection problem
	 */
	const AccumulatedValues& getAccumulatedProjectorValues() const
	{ return projectorDB.getAccumulatedValues(); }

	/** Getter for the accumulated values of the minimization inner problem.
	 *
	 * This is used by the parallel implementation to obtain the values
	 * for the the Slave nodes to transmit to Master node.
	 *
	 * @return ref to accumulated values from minimization problem
	 */
	const AccumulatedValues& getAccumulatedSolverValues() const
	{ return solverDB.getAccumulatedValues(); }

	/** Function to add values accumulated in other solvers, e.g. due to parallel
	 * runs, inside this solver's databases.
	 *
	 * @param _values_first start iterator
	 * @param _values_last end iterator
	 */
	void insertProjectorValues(
			const std::vector<AccumulatedValues>::const_iterator &_values_first,
			const std::vector<AccumulatedValues>::const_iterator &_values_last
			)
	{ projectorDB.insertValues(_values_first,_values_last); }

	/** Function to add values accumulated in other solvers, e.g. due to parallel
	 * runs, inside this solver's databases.
	 *
	 * @param _values_first start iterator
	 * @param _values_last end iterator
	 */
	void insertSolverValues(
			const std::vector<AccumulatedValues>::const_iterator &_values_first,
			const std::vector<AccumulatedValues>::const_iterator &_values_last
			)
	{ solverDB.insertValues(_values_first,_values_last); }

private:
	const CommandLineOptions &opts;
	//!> table database for projector iteration information
	Database_ptr_t projector_db;
	//!> table database for solver iteration information
	Database_ptr_t solver_db;

	InnerProblemDatabase &projectorDB;
	InnerProblemDatabase &solverDB;

	//!> replacement to use as delta when performing projection onto range
	const double projection_delta;
};



#endif /* MATRIXFACTORIZERBASE_SOLVERS_INRANGESOLVER_HPP_ */
