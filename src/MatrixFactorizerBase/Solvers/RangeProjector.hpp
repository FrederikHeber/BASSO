/*
 * RangeProjector.hpp
 *
 *  Created on: Jun 26, 2015
 *      Author: heber
 */

#ifndef MATRIXFACTORIZERBASE_SOLVERS_RANGEPROJECTOR_HPP_
#define MATRIXFACTORIZERBASE_SOLVERS_RANGEPROJECTOR_HPP_

#include "BassoConfig.h"

#include <Eigen/Dense>

#include "Minimizations/types.hpp"

#include "MatrixFactorizerBase/Options/MatrixFactorizerOptions.hpp"

/** This class contains all code for projecting a give vector onto the
 * range of a given matrix.
 */
struct RangeProjector
{
	/** Constructor for class RangeProjector.
	 *
	 * @param _database ref to database
	 * @param _opts options defining manner of inverse problem solving
	 */
	RangeProjector(
			Database_ptr_t &_database,
			const MatrixFactorizerOptions &_opts
			);

	/** Functor that projects given \a _rhs onto \a _matrix by solving
	 * the inverse problem.
	 *
	 * @param _matrix matrix onto whose range to project
	 * @param _rhs vector to project
	 * @param _resultingvalue projected vector
	 * @param _nonnegative controls whether solution must be non-negative
	 * @return true - success
	 */
	bool operator()(
			const Eigen::MatrixXd &_matrix,
			const Eigen::VectorXd &_rhs,
			Eigen::VectorXd &_resultingvalue,
			const bool _nonnegative
			);

private:
	//!> database reference for storing iteration information
	Database_ptr_t database;
	//!> storing options on the manner of solving
	const MatrixFactorizerOptions opts;
};


#endif /* MATRIXFACTORIZERBASE_SOLVERS_RANGEPROJECTOR_HPP_ */
