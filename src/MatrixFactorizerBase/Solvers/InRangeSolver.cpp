/*
 * InRangeSolver.cpp
 *
 *  Created on: Jun 26, 2015
 *      Author: heber
 */


#include "BassoConfig.h"

#include "InRangeSolver.hpp"

#include "Log/Logging.hpp"

#include "Database/Database_mock.hpp"
#include "MatrixFactorizerBase/Options/MatrixFactorizerOptions.hpp"
#include "MatrixFactorizerBase/Solvers/InverseProblemSolver.hpp"
#include "MatrixFactorizerBase/Solvers/RangeProjector.hpp"

InRangeSolver::InRangeSolver(const MatrixFactorizerOptions &_opts) :
	opts(_opts),
	mock_db(new Database_mock)
{}


bool InRangeSolver::operator()(
		const Eigen::MatrixXd &_matrix,
		const Eigen::VectorXd &_rhs,
		const Eigen::VectorXd &_solution_start,
		Eigen::VectorXd &_solution,
		const unsigned int _dim,
		const unsigned int _loop_nr
		)
{
	Eigen::VectorXd projected_rhs(_rhs);
	BOOST_LOG_TRIVIAL(debug)
		<< "------------------------ n=" << _dim << " ------------------";
	// project right-hand side onto range of matrix
	BOOST_LOG_TRIVIAL(trace)
			<< "Initial y_" << _dim << " is "
			<< projected_rhs.transpose();
	{
		// project y onto image of K
		RangeProjector projector(
				mock_db,
				opts
				);
		if (!projector(
				_matrix,
				_rhs,
				projected_rhs,
				false))
			return false;
	}
	BOOST_LOG_TRIVIAL(trace)
			<< "Projected y_" << _dim << " is "
			<< projected_rhs.transpose();
	BOOST_LOG_TRIVIAL(trace)
			<< "Difference y_" << _dim << "-y'_" << _dim << " is "
			<< (_rhs-projected_rhs).transpose();
	BOOST_LOG_TRIVIAL(debug)
			<< "........................ n=" << _dim << " ..................";

	// solve inverse problem for projected right-hand side
	_solution = _solution_start;
	InverseProblemSolver solver(
			mock_db,
			opts,
			true /* true solution calculation */);
	if (!solver(
			_matrix,
			projected_rhs,
			_solution_start,
			_solution,
			_loop_nr >= 3))
		return false;
	BOOST_LOG_TRIVIAL(trace)
		<< "Resulting vector is " << _solution.transpose();

	return true;
}


