/*
 * InRangeSolver.cpp
 *
 *  Created on: Jun 26, 2015
 *      Author: heber
 */


#include "BassoConfig.h"

#include "InRangeSolver.hpp"

#include "Log/Logging.hpp"

#include "MatrixFactorizer/Helpers/detail.hpp"
#include "Options/CommandLineOptions.hpp"
#include "RangeProjector/RangeProjector/RangeProjectionSolver.hpp"
#include "Solvers/InverseProblemSolver.hpp"

InRangeSolver::InRangeSolver(
		const CommandLineOptions &_opts,
		const InnerProblemDatabase::keys_t &_overall_keys
		) :
	opts(_opts),
	projector_db(new InnerProblemDatabase),
	solver_db(new InnerProblemDatabase(_overall_keys)),
	projectorDB(static_cast<InnerProblemDatabase &>(*projector_db.get())),
	solverDB(static_cast<InnerProblemDatabase &>(*solver_db.get()))
{}

bool InRangeSolver::operator()(
		const Eigen::MatrixXd &_matrix,
		const Eigen::VectorXd &_rhs,
		const Eigen::VectorXd &_solution_start,
		Eigen::VectorXd &_solution,
		const unsigned int _dim
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
		RangeProjectionSolver projector(
				projector_db,
				opts
				);
		if (!projector(
				_matrix,
				_rhs,
				projected_rhs))
			return false;
		projector_db->clear();
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
			solver_db,
			opts,
			false /* true solution calculation */);
	if (!solver(
			_matrix,
			projected_rhs,
			_solution_start,
			_solution)) {
		BOOST_LOG_TRIVIAL(trace)
			<< "Resulting vector is " << _solution.transpose();
		return false;
	}
	solver_db->clear();


	return true;
}

void InRangeSolver::insertAccumulatedProjectorValues(
		Table &_table) const
{
	projectorDB.insertAccumulatedValues(_table);
}

void InRangeSolver::insertAccumulatedSolverValues(
		Table &_table) const
{
	solverDB.insertAccumulatedValues(_table);
}


