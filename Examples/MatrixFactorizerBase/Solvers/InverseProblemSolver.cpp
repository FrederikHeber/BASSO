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
#include "Minimizations/Minimizers/GeneralMinimizer.hpp"
#include "Minimizations/Minimizers/MinimizerFactory.hpp"
#include "Minimizations/Minimizers/MinimizationExceptions.hpp"
#include "Options/CommandLineOptions.hpp"
#include "SolverFactory/SolverFactory.hpp"

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
	GeneralMinimizer::ReturnValues result;
	// prepare inverse problem
	InverseProblem_ptr_t inverseproblem =
			SolverFactory::createInverseProblem(
					opts, _matrix, _rhs);
	MinimizerFactory::instance_ptr_t minimizer =
			SolverFactory::createMinimizer(
					opts, inverseproblem, database, opts.maxiter);
	if (minimizer == NULL) {
		BOOST_LOG_TRIVIAL(error)
				<< "Minimizer could not be constructed, exiting.";
		return false;
	}

	SpaceElement_ptr_t truesolution;
	if (checkTrueSolution) {
		// SVD only gives true solution for l2 norm
		assert( inverseproblem->y->getSpace()->getNorm()->getPvalue() == 2. );

		// empty or true solution from diagonalization
		Eigen::MatrixXd copymatrix = _matrix;
		Eigen::JacobiSVD<Eigen::MatrixXd> svd =
				copymatrix.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
		const Eigen::VectorXd truesolution_vector =
				svd.solve(_rhs);
		BOOST_LOG_TRIVIAL(trace)
				<< "True solution is " << truesolution_vector.transpose()
				<< " with norm "
				<< (_matrix*truesolution_vector - _rhs).norm()/_rhs.norm();
		truesolution =
				ElementCreator::create(
						inverseproblem->x->getSpace(),
						truesolution_vector);
	} else {
		truesolution =
				inverseproblem->x->getSpace()->createElement();
	}

	// prepare start value and dual solution
	SpaceElement_ptr_t x0 = ElementCreator::create(
			inverseproblem->x->getSpace(),
			_startingvalue);
	*inverseproblem->x = x0;
	if (x0->getSpace()->getDimension() < 10)
		BOOST_LOG_TRIVIAL(debug)
			<< "Starting at x0 = " << x0;
	SpaceElement_ptr_t dualx0;
	if (opts.type_spacex == "lp") {
		dualx0 = (*inverseproblem->x->getSpace()->getDualityMapping())(x0);
	} else {
		dualx0 = inverseproblem->x->getSpace()->getDualSpace()->createElement();
		dualx0->setZero();
	}

	// and minimize
	try{
		result = (*minimizer)(
						inverseproblem,
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
	detail::setResultingVector(result.m_solution, _solution, _nonnegative);

	return true;
}

