/*
 * InverseProblemSolver.cpp
 *
 *  Created on: Jun 26, 2015
 *      Author: heber
 */


#include "BassoConfig.h"

#include "InverseProblemSolver.hpp"

#include <cassert>

#include "Log/Logging.hpp"
#include "MatrixFactorizerBase/Helpers/detail.hpp"
#include "Minimizations/Elements/ElementCreator.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Mappings/LinearMapping.hpp"
#include "Minimizations/Mappings/SingularValueDecomposition.hpp"
#include "Minimizations/Minimizers/GeneralMinimizer.hpp"
#include "Minimizations/Minimizers/MinimizerFactory.hpp"
#include "Minimizations/Minimizers/MinimizationExceptions.hpp"
#include "Options/CommandLineOptions.hpp"
#include "Solvers/SolverFactory/SolverFactory.hpp"

InverseProblemSolver::InverseProblemSolver(
		Database_ptr_t &_database,
		const CommandLineOptions &_opts,
		const bool _checkTrueSolution
		) :
		database(_database),
		opts(_opts),
		checkTrueSolution(_checkTrueSolution)
{}

bool InverseProblemSolver::operator()(
		const Eigen::MatrixXd &_matrix,
		const Eigen::MatrixXd &_rhs,
		const Eigen::VectorXd &_startingvalue,
		Eigen::VectorXd &_solution,
		const bool _nonnegative
		)
{
	// prepare inverse problem
	InverseProblem_ptr_t inverseproblem =
			SolverFactory::createInverseProblem(
					opts, _matrix, _rhs);

	bool status = operator()(inverseproblem, _startingvalue);
	detail::setResultingVector(inverseproblem->x, _solution, _nonnegative);
	return status;
}

bool InverseProblemSolver::operator()(
		InverseProblem_ptr_t &_inverseproblem,
		const Eigen::VectorXd &_startingvalue
		)
{
	GeneralMinimizer::ReturnValues result;
	MinimizerFactory::instance_ptr_t minimizer =
			SolverFactory::createMinimizer(
					opts, _inverseproblem, database, opts.maxiter);

	if (minimizer == NULL) {
		BOOST_LOG_TRIVIAL(error)
				<< "Minimizer could not be constructed, exiting.";
		return false;
	}

	SpaceElement_ptr_t truesolution;
	if (checkTrueSolution) {
		const LinearMapping &A = static_cast<LinearMapping&>(*_inverseproblem->A);
		const SpaceElement_ptr_t &rhs = _inverseproblem->y;
		SingularValueDecomposition svd = A.getSVD();
		const SpaceElement_ptr_t truesolution = svd.solve(rhs);

		// empty or true solution from diagonalization
		BOOST_LOG_TRIVIAL(trace)
				<< "True solution is " << *truesolution
				<< " with norm "
				<< (A(truesolution) - rhs)->Norm()/rhs->Norm();
	} else {
		truesolution =
				_inverseproblem->x->getSpace()->createElement();
	}

	// prepare start value and dual solution
	SpaceElement_ptr_t x0 = ElementCreator::create(
			_inverseproblem->x->getSpace(),
			_startingvalue);
	*_inverseproblem->x = x0;
	if (x0->getSpace()->getDimension() < 10)
		BOOST_LOG_TRIVIAL(debug)
			<< "Starting at x0 = " << x0;
	SpaceElement_ptr_t dualx0;
	if (opts.type_spacex == "lp") {
		dualx0 = (*_inverseproblem->x->getSpace()->getDualityMapping())(x0);
	} else {
		dualx0 = _inverseproblem->x->getSpace()->getDualSpace()->createElement();
		dualx0->setZero();
	}

	// and minimize
	try{
		result = (*minimizer)(
						_inverseproblem,
						x0,
						dualx0,
						truesolution);
		minimizer->resetState();
	} catch (MinimizationIllegalValue_exception &e) {
		std::cerr << "Illegal value for "
				<< *boost::get_error_info<MinimizationIllegalValue_name>(e)
				<< std::endl;
		return false;
	}
	assert( *result.m_solution == *_inverseproblem->x );

	return true;
}

