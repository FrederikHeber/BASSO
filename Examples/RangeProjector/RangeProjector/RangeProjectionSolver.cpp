/*
 * RangeProjectionSolver.cpp
 *
 *  Created on: Jun 26, 2015
 *      Author: heber
 */


#include "BassoConfig.h"

#include <boost/assign.hpp>

#include "RangeProjectionSolver.hpp"

#include "Log/Logging.hpp"
#include "Math/Helpers.hpp"
#include "MatrixFactorizer/Helpers/detail.hpp"
#include "Minimizations/Elements/ElementCreator.hpp"
#include "Minimizations/Elements/RepresentationAdvocate.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/InverseProblems/InverseProblem.hpp"
#include "Minimizations/Mappings/LinearMappingFactory.hpp"
#include "Minimizations/Minimizers/GeneralMinimizer.hpp"
#include "Minimizations/Minimizers/MinimizerFactory.hpp"
#include "Minimizations/Minimizers/MinimizationExceptions.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"
#include "Minimizations/Spaces/NormedSpaceFactory.hpp"
#include "Minimizations/Minimizers/StoppingCriteria/StoppingCriteriaFactory.hpp"
#include "Solvers/SolverFactory/SolverFactory.hpp"

using namespace boost::assign;

RangeProjectionSolver::RangeProjectionSolver(
		Database_ptr_t &_database,
		const CommandLineOptions &_opts
		) :
		database(_database),
		opts(_opts)
{
	// use smaller delta for the projection and SESOP
	const_cast<CommandLineOptions &>(opts).delta = 1e-8;
	const_cast<CommandLineOptions &>(opts).algorithm_name =
			MinimizerFactory::TypeNames[MinimizerFactory::sequentialsubspace];
}

bool RangeProjectionSolver::operator()(
		const Eigen::MatrixXd &_matrix,
		const Eigen::VectorXd &_rhs,
		Eigen::VectorXd &_resultingvalue
		)
{
	// prepare right-hand side
	NormedSpace_ptr_t Y;
	NormedSpace_ptr_t Ys;
	{
		NormedSpaceFactory::args_t args;
		args += boost::any(opts.py), boost::any(opts.powery);
		Y = NormedSpaceFactory::create(
				_matrix.innerSize(), opts.type_spacey, args);
		Ys = Y->getDualSpace();
	}
	NormedSpace_ptr_t X;
	NormedSpace_ptr_t Xs;
	{
		NormedSpaceFactory::args_t args;
		args += boost::any(opts.px), boost::any(opts.powerx);
		X = NormedSpaceFactory::create(
				_matrix.outerSize(), opts.type_spacex, args);
		Xs = X->getDualSpace();
	}
	// and the LinearMapping
	Mapping_ptr_t As =
			LinearMappingFactory::createInstance(Ys,Xs,_matrix.transpose());

	SpaceElement_ptr_t rhs = ElementCreator::create(Y, _rhs);
	SpaceElement_ptr_t dualrhs =
			(*Y->getDualityMapping())((-1.)*rhs);
	SpaceElement_ptr_t dualmappedrhs = (-1.)*((*As)(dualrhs));

	// prepare inverse problem: Y^\ast \rightarrow X^\ast
	InverseProblem_ptr_t inverseproblem(
			new InverseProblem(As,Ys,Xs,dualmappedrhs) );

	// create stopping criterion
	StoppingCriteriaFactory stop_factory;
	StoppingCriterion::ptr_t stopping_criterion =
			stop_factory.create(opts.stopping_criteria, opts.stopping_args);

	// prepare minimizer
	MinimizerFactory::instance_ptr_t minimizer =
			SolverFactory::createMinimizer(
					opts, inverseproblem, database, stopping_criterion, opts.maxiter);
	if (minimizer == NULL) {
		BOOST_LOG_TRIVIAL(error)
				<< "Minimizer could not be constructed, exiting.";
		return false;
	}

	// empty true solution
	SpaceElement_ptr_t truesolution =
			inverseproblem->x->getSpace()->createElement();
	truesolution->setZero();

	// prepare start value and dual solution
	SpaceElement_ptr_t dualy0 =
			inverseproblem->x->getSpace()->createElement();
	dualy0->setZero();
	*inverseproblem->x = dualy0;
	if (dualy0->getSpace()->getDimension() < 10)
		BOOST_LOG_TRIVIAL(debug)
			<< "Starting at dualy0 = " << dualy0;
	SpaceElement_ptr_t y0;
	// only for smooth spaces we may use the duality mapping
	if (inverseproblem->x->getSpace()->getNorm()->isSmooth()) {
		y0 = (*inverseproblem->x->getSpace()->getDualityMapping())(dualy0);
	} else {
		y0 = inverseproblem->x->getSpace()->getDualSpace()->createElement();
		y0->setZero();
	}

	// and minimize
	{
		GeneralMinimizer::ReturnValues result;
		try{
			result = (*minimizer)(
							inverseproblem,
							dualy0,
							y0,
							truesolution);
			minimizer->resetState();
		} catch (MinimizationIllegalValue_exception &e) {
			std::cerr << "Illegal value for "
					<< *boost::get_error_info<MinimizationIllegalValue_name>(e)
					<< std::endl;
			return false;
		}
		*result.m_solution += dualrhs;
		SpaceElement_ptr_t projected_solution =
				rhs +
				(*inverseproblem->x->getSpace()->getDualityMapping())
					(result.m_solution);
		_resultingvalue = RepresentationAdvocate::get(projected_solution);
	}

	return true;
}

