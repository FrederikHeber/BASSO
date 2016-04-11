/*
 * RangeProjectionSolver.cpp
 *
 *  Created on: Jun 26, 2015
 *      Author: heber
 */


#include "BassoConfig.h"

#include <boost/assign.hpp>

#include "RangeProjectionSolver.hpp"

#include "Database/Database.hpp"
#include "Log/Logging.hpp"
#include "Math/Helpers.hpp"
#include "MatrixFactorizer/Helpers/detail.hpp"
#include "Minimizations/Elements/ElementCreator.hpp"
#include "Minimizations/Elements/RepresentationAdvocate.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/InverseProblems/InverseProblem.hpp"
#include "Minimizations/Mappings/LinearMappingFactory.hpp"
#include "Minimizations/Minimizers/GeneralMinimizer.hpp"
#include "Minimizations/Minimizers/MinimizationExceptions.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"
#include "Minimizations/Spaces/NormedSpaceFactory.hpp"
#include "Solvers/SolverFactory/SolverFactory.hpp"

using namespace boost::assign;

RangeProjectionSolver::RangeProjectionSolver(
		const Eigen::MatrixXd &_matrix,
		const Eigen::VectorXd &_rhs,
		Database_ptr_t &_database,
		const CommandLineOptions &_opts
		) :
		database(_database),
		opts(_opts),
		name("RangeProjection")
{
	// prepare spaces
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

	// prepare LinearMapping
	Mapping_ptr_t As =
			LinearMappingFactory::createInstance(
					Ys,Xs,_matrix.transpose());

	// prepare right-hand side
	rhs = ElementCreator::create(Y, _rhs);
	dualrhs = (*Y->getDualityMapping())((-1.)*rhs);
	SpaceElement_ptr_t dualmappedrhs = (-1.)*((*As)(dualrhs));

	// prepare inverse problem: Y^\ast \rightarrow X^\ast
	inverseproblem.reset(
			new InverseProblem(As,Ys,Xs,dualmappedrhs));

	// prepare minimizer
	minimizer = SolverFactory::createMinimizer(
					opts, inverseproblem, database);
}


GeneralMinimizer::ReturnValues RangeProjectionSolver::operator()(
		const SpaceElement_ptr_t &_startingvalue,
		const SpaceElement_ptr_t)
{
	GeneralMinimizer::ReturnValues result;
	result.residuum = 0.;

	if (minimizer == NULL) {
		BOOST_LOG_TRIVIAL(error)
				<< "Minimizer could not be constructed, exiting.";
		result.status = GeneralMinimizer::ReturnValues::error;
		return result;
	}

	// empty true solution
	SpaceElement_ptr_t truesolution =
			inverseproblem->x->getSpace()->createElement();
	truesolution->setZero();

	// init starting value for this loop
	result.m_solution = _startingvalue->getSpace()->createElement();
	*result.m_solution = _startingvalue;
	*inverseproblem->x = result.m_solution;
	if (result.m_solution->getSpace()->getDimension() < 10)
		BOOST_LOG_TRIVIAL(debug)
			<< "Starting at dualy0 = " << result.m_solution;

	// only for smooth spaces we may use the duality mapping
	if (inverseproblem->x->getSpace()->getNorm()->isSmooth()) {
		result.m_dual_solution =
				(*inverseproblem->x->getSpace()->getDualityMapping())(result.m_solution);
	} else {
		result.m_dual_solution =
				inverseproblem->x->getSpace()->getDualSpace()->createElement();
		result.m_dual_solution->setZero();
	}

	// and minimize
	{
		try{
			result = (*minimizer)(
							inverseproblem,
							result.m_solution,
							result.m_dual_solution,
							truesolution);
			minimizer->resetState();
		} catch (MinimizationIllegalValue_exception &e) {
			std::cerr << "Illegal value for "
					<< *boost::get_error_info<MinimizationIllegalValue_name>(e)
					<< std::endl;
			result.status = GeneralMinimizer::ReturnValues::error;
			return result;
		}
		*result.m_solution += dualrhs;
		SpaceElement_ptr_t projected_solution =
				rhs +
				(*inverseproblem->x->getSpace()->getDualityMapping())
					(result.m_solution);
		*result.m_dual_solution = projected_solution;
	}

	return result;
}

SpaceElement_ptr_t RangeProjectionSolver::getZeroStartvalue() const
{
	SpaceElement_ptr_t x0 =
			inverseproblem->x->getSpace()->createElement();
	x0->setZero();
	return x0;
}


void RangeProjectionSolver::clear()
{
	database->clear();
}

void RangeProjectionSolver::finish()
{
	database->finish();
}
