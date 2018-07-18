/*
 * RangeProjectionSolver.hpp
 *
 *  Created on: Jun 26, 2015
 *      Author: heber
 */

#ifndef RANGEPROJECTORBASE_RANGEPROJECTIONSOLVER_HPP_
#define RANGEPROJECTORBASE_RANGEPROJECTIONSOLVER_HPP_

#include "BassoConfig.h"

#include <Eigen/Dense>

#include "Minimizations/types.hpp"
#include "Minimizations/Minimizers/GeneralMinimizer.hpp"
#include "Minimizations/Minimizers/MinimizerFactory.hpp"
#include "Minimizations/Minimizers/StoppingCriteria/StoppingCriteriaFactory.hpp"
#include "Options/CommandLineOptions.hpp"
#include "Solvers/FeasibilityProblem.hpp"

/** This class contains all code for projecting a give vector onto the
 * range of a given matrix.
 */
struct RangeProjectionSolver : public FeasibilityProblem
{
	/** Constructor for class RangeProjectionSolver.
	 *
	 * @param _matrix matrix onto whose range to project
	 * @param _rhs vector to project
	 * @param _database ref to database
	 * @param _opts options defining manner of inverse problem solving
	 */
	RangeProjectionSolver(
			const Eigen::MatrixXd &_matrix,
			const Eigen::VectorXd &_rhs,
			Database_ptr_t &_database,
			const CommandLineOptions &_opts
			);

	virtual ~RangeProjectionSolver() {}

	/** Operator to solve this specific feasibility problem.
	 *
	 * @param _startingvalue start value or zero
	 * @return result container with success state and solutions
	 */
	GeneralMinimizer::ReturnValues operator()(
			const SpaceElement_ptr_t &_startingvalue)
	{ return operator()(_startingvalue, SpaceElement_ptr_t()); }

	/** Operator to solve this specific feasibility problem while
	 * comparing to a true solution.
	 *
	 * @param _startingvalue start value or zero
	 * @param _truesolution true solution for true error calculation
	 * @return result container with success state and solutions
	 */
	GeneralMinimizer::ReturnValues operator()(
			const SpaceElement_ptr_t &_startingvalue,
			const SpaceElement_ptr_t _truesolution);

	/** Returns a zero start value using the internal \a inverseproblem.
	 *
	 * @return zero start value
	 */
	SpaceElement_ptr_t getZeroStartvalue() const;

	/** Returns specific name of the problem for output.
	 *
	 * @return name of the problem
	 */
	const std::string& getName() const
	{ return name; }

	/** Resets the feasibility problem.
	 *
	 */
	void clear();

	/** Finalizes the feasibility problem.
	 *
	 */
	void finish();

private:
	//!> database reference for storing iteration information
	Database_ptr_t database;
	//!> storing options on the manner of solving
	const CommandLineOptions &opts;

private:
	// stuff constructed prior to operator()

	//!> right-hand side
	SpaceElement_ptr_t rhs;
	//!> dual of right-hand side
	SpaceElement_ptr_t dualrhs;
	//!> inverse problem
	InverseProblem_ptr_t inverseproblem;
	//!> pre-constructed stopping criteria
	StoppingCriterion::ptr_t stopping_criterion;
	//!> pre-constructed minimizer
	MinimizerFactory::instance_ptr_t minimizer;

	const std::string name;
};


#endif /* RANGEPROJECTORBASE_RANGEPROJECTIONSOLVER_HPP_ */
