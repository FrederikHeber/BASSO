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

#include "Options/CommandLineOptions.hpp"

/** This class contains all code for projecting a give vector onto the
 * range of a given matrix.
 */
struct RangeProjectionSolver
{
	/** Constructor for class RangeProjectionSolver.
	 *
	 * @param _database ref to database
	 * @param _opts options defining manner of inverse problem solving
	 */
	RangeProjectionSolver(
			Database_ptr_t &_database,
			const CommandLineOptions &_opts
			);

	/** Functor that projects given \a _rhs onto \a _matrix by solving
	 * the inverse problem.
	 *
	 * @param _matrix matrix onto whose range to project
	 * @param _rhs vector to project
	 * @param _resultingvalue projected vector
	 * @return true - success
	 */
	bool operator()(
			const Eigen::MatrixXd &_matrix,
			const Eigen::VectorXd &_rhs,
			Eigen::VectorXd &_resultingvalue
			);

private:
	//!> database reference for storing iteration information
	Database_ptr_t database;
	//!> storing options on the manner of solving
	const CommandLineOptions opts;
};


#endif /* RANGEPROJECTORBASE_RANGEPROJECTIONSOLVER_HPP_ */