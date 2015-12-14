/*
 * VectorProjection.cpp
 *
 *  Created on: Sep 4, 2014
 *      Author: heber
 */

#include "BassoConfig.h"

#include "VectorProjection.hpp"

#include <utility>

#include "Log/Logging.hpp"
#include "Minimizations/Functions/BregmanDistance.hpp"
#include "Minimizations/Functions/Minimizers/FunctionalMinimizerFactory.hpp"
#include "Minimizations/Functions/VectorProjection_BregmanDistanceToLine.hpp"

VectorProjection::VectorProjection(
		const Norm &_lpnorm,
		const Mapping &_J_p,
		const double _p
		) :
		lpnorm(_lpnorm),
		J_p(_J_p),
		p(_p)
{}

const std::pair<double, double>
VectorProjection::operator()(
		const SpaceElement_ptr_t &_projectedonto,
		const SpaceElement_ptr_t &_tobeprojected,
		const double _Tol) const
{
	BregmanDistance distance(
			lpnorm,
			J_p,
			p
			);


	const unsigned int dim = 1;
	double tmin = 0.;
	double result = 0.;
	{

		// tmin=fminunc(@(t) BregmanFunctional(t,Jx,u,alpha+d,DualNormX,DualPowerX,TolX),t0,BregmanOptions);
		const BregmanDistanceToLine distancefunctional(
				distance,
				lpnorm,
				J_p,
				_projectedonto,
				_tobeprojected,
				p);
		FunctionalMinimizer<double>::ptr_t functionminimizer =
				FunctionalMinimizerFactory::create<double>(
						dim,
						distancefunctional,
						tmin);

//		const unsigned int inner_iterations =
				(*functionminimizer)(dim, _Tol, tmin);
		result = functionminimizer->getCurrentOptimumValue();

		BOOST_LOG_TRIVIAL(trace)
			<< "tmin is " << tmin;
	}

	return std::make_pair( result, tmin);
}
