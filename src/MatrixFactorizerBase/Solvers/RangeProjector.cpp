/*
 * RangeProjector.cpp
 *
 *  Created on: Jun 26, 2015
 *      Author: heber
 */


#include "BassoConfig.h"

#include "RangeProjector.hpp"

#include "Log/Logging.hpp"
#include "Math/Helpers.hpp"
#include "MatrixFactorizerBase/Helpers/detail.hpp"
#include "Minimizations/Elements/ElementCreator.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/InverseProblems/InverseProblem.hpp"
#include "Minimizations/Mappings/LinearMappingFactory.hpp"
#include "Minimizations/Minimizers/GeneralMinimizer.hpp"
#include "Minimizations/Minimizers/MinimizerFactory.hpp"
#include "Minimizations/Minimizers/MinimizationExceptions.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"
#include "Minimizations/Spaces/NormedSpaceFactory.hpp"
#include "SolutionFactory/SolutionFactory.hpp"

RangeProjector::RangeProjector(
		Database_ptr_t &_database,
		const MatrixFactorizerOptions &_opts
		) :
		database(_database),
		opts(_opts)
{
	// use smaller delta for the projection and SESOP
	const_cast<MatrixFactorizerOptions &>(opts).delta = 1e-8;
	const_cast<MatrixFactorizerOptions &>(opts).algorithm_name =
			MinimizerFactory::TypeNames[MinimizerFactory::sequentialsubspace];
}

bool RangeProjector::operator()(
		const Eigen::MatrixXd &_matrix,
		const Eigen::VectorXd &_rhs,
		Eigen::VectorXd &_resultingvalue,
		const bool _nonnegative
		)
{
	// require dual values
	const double dualnormx =
			Helpers::ConjugateValue(opts.px);
	const double dualnormy =
			Helpers::ConjugateValue(opts.py);
	const double dualpowerx =
			Helpers::ConjugateValue(opts.powerx);
	const double dualpowery =
			Helpers::ConjugateValue(opts.powery);

	// prepare right-hand side
	NormedSpace_ptr_t Ys =
			NormedSpaceFactory::createLpInstance(
					_matrix.innerSize(), dualnormy, dualpowery);
	NormedSpace_ptr_t Xs =
			NormedSpaceFactory::createLpInstance(
					_matrix.outerSize(), dualnormx, dualpowerx);
	// and the LinearMapping
	Mapping_ptr_t As =
			LinearMappingFactory::createInstance(Ys,Xs,_matrix.transpose());
	const NormedSpace &Y = *Ys->getDualSpace();

	SpaceElement_ptr_t rhs = ElementCreator::create(Y, _rhs);
	SpaceElement_ptr_t dualrhs =
			(*Y.getDualityMapping())((-1.)*rhs);
	SpaceElement_ptr_t dualmappedrhs = (-1.)*((*As)(dualrhs));

	// prepare inverse problem: Y^\ast \rightarrow X^\ast
	InverseProblem_ptr_t inverseproblem(
			new InverseProblem(As,Ys,Xs,dualmappedrhs) );

	// prepare minimizer
	MinimizerFactory::instance_ptr_t minimizer =
			SolutionFactory::createMinimizer(
					opts, inverseproblem, database, opts.inner_iterations);
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
	if (opts.type_spacex == "lp") {
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
		detail::setResultingVector(projected_solution, _resultingvalue, _nonnegative);
	}

	return true;
}
