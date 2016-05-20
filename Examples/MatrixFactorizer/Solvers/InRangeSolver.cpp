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
#include "Minimizations/Elements/ElementCreator.hpp"
#include "Options/CommandLineOptions.hpp"
#include "Solvers/FeasibilityProblem.hpp"
#include "Solvers/InverseProblemSolver.hpp"
#include "Solvers/RangeProjectionSolver.hpp"
#include "Solvers/SolverFactory/SolverFactory.hpp"
#include "Solvers/SplitFeasibilitySolver.hpp"

InRangeSolver::InRangeSolver(
		const CommandLineOptions &_opts,
		const InnerProblemDatabase::keys_t &_overall_keys,
		const double _projection_delta
		) :
	opts(_opts),
	projector_db(new InnerProblemDatabase(_overall_keys)),
	solver_db(new InnerProblemDatabase(_overall_keys)),
	projectorDB(static_cast<InnerProblemDatabase &>(*projector_db.get())),
	solverDB(static_cast<InnerProblemDatabase &>(*solver_db.get())),
	projection_delta(_projection_delta)
{}

InRangeSolver::InRangeSolver(const InRangeSolver &other) :
	opts(other.opts),
	projector_db(other.projector_db),
	solver_db(other.solver_db),
	projectorDB(static_cast<InnerProblemDatabase &>(*projector_db.get())),
	solverDB(static_cast<InnerProblemDatabase &>(*solver_db.get())),
	projection_delta(other.projection_delta)
{}

bool InRangeSolver::operator()(
		const Eigen::MatrixXd &_matrix,
		const Eigen::VectorXd &_rhs,
		const Eigen::VectorXd &_solution_start,
		Eigen::VectorXd &_solution,
		const unsigned int _dim,
		const AuxiliaryConstraints::ptr_t &_auxiliary_constraints
		)
{
	Eigen::VectorXd projected_rhs;
	LOG(debug, "------------------------ n=" << _dim << " ------------------");
	// project right-hand side onto range of matrix
	BOOST_LOG_TRIVIAL(trace)
			<< "Initial y_" << _dim << " is "
			<< _rhs.transpose();
	{
		// use smaller delta for the projection and SESOP
		CommandLineOptions projection_opts(opts);
		projection_opts.delta = projection_delta;
		projection_opts.algorithm_name =
				MinimizerFactory::TypeNames[MinimizerFactory::sequentialsubspace];

		RangeProjectionSolver projector(
				_matrix,
				_rhs,
				projector_db,
				projection_opts
				);

		// zero start value
		SpaceElement_ptr_t dualy0 = projector.getZeroStartvalue();

		// solve
		projector_db->clear();
		GeneralMinimizer::ReturnValues result =
				projector(dualy0);
		projector_db->finish();

		// get projected rhs
		if (result.status == GeneralMinimizer::ReturnValues::finished)
			projected_rhs = RepresentationAdvocate::get(result.m_dual_solution);
		else
			return false;
	}
	BOOST_LOG_TRIVIAL(trace)
			<< "Projected y_" << _dim << " is "
			<< projected_rhs.transpose();
	BOOST_LOG_TRIVIAL(trace)
			<< "Difference y_" << _dim << "-y'_" << _dim << " is "
			<< (_rhs-projected_rhs).transpose();

	LOG(debug, "........................ n=" << _dim << " ..................");

	// prepare inverse problem
	InverseProblem_ptr_t inverseproblem =
			SolverFactory::createInverseProblem(
					opts, _matrix, projected_rhs);

	FeasibilityProblem::ptr_t solver;
	if (_auxiliary_constraints) {
		SplitFeasibilitySolver *SFP =
				new SplitFeasibilitySolver(opts);
		FeasibilityProblem::ptr_t IP(
				new InverseProblemSolver(
						inverseproblem,
						solver_db,
						opts,
						false /* true solution calculation */)
				);
		SFP->registerFeasibilityProblem(IP);
		SFP->registerAuxiliaryConstraints(_auxiliary_constraints);
		solver.reset(SFP);
	} else {
		solver.reset(
				new InverseProblemSolver(
						inverseproblem,
						solver_db,
						opts,
						false /* true solution calculation */)
				);
	}

	// solve inverse problem for projected right-hand side
	SpaceElement_ptr_t solution_start = ElementCreator::create(
			inverseproblem->x->getSpace(),
			_solution_start);
	solver_db->clear();
	GeneralMinimizer::ReturnValues result =
			(*solver)(solution_start);
	solver_db->finish();

	if (result.status == GeneralMinimizer::ReturnValues::finished) {
		_solution = RepresentationAdvocate::get(result.m_solution);
		LOG(trace, "Resulting vector is " << _solution.transpose());
		return true;
	} else
		return false;
}

void InRangeSolver::insertAccumulatedProjectorValues(
		Table &_table,
		const std::string &_suffix) const
{
	projectorDB.insertAccumulatedValues(_table, _suffix);
}

void InRangeSolver::insertAccumulatedSolverValues(
		Table &_table,
		const std::string &_suffix) const
{
	solverDB.insertAccumulatedValues(_table, _suffix);
}


