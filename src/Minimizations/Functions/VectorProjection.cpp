/*
 * VectorProjection.cpp
 *
 *  Created on: Sep 4, 2014
 *      Author: heber
 */

#include "BassoConfig.h"

#include "VectorProjection.hpp"

#include <boost/log/trivial.hpp>

#include "Minimizations/Functions/BregmanDistance.hpp"
#include "Minimizations/Functions/FunctionMinimizer.hpp"
#include "Minimizations/Functions/VectorProjection_BregmanDistanceToLine.hpp"

VectorProjection::VectorProjection(
		const Norm &_lpnorm,
		const PowerTypeDualityMapping &_J_p,
		const double _p,
		const OperationCounter<
								Eigen::internal::scalar_product_traits<typename Eigen::internal::traits<Eigen::VectorXd>::Scalar, typename Eigen::internal::traits<Eigen::VectorXd>::Scalar>::ReturnType,
								const Eigen::MatrixBase<Eigen::VectorXd>&,
								const Eigen::MatrixBase<Eigen::VectorXd>&
								>& _ScalarVectorProduct
		) :
		lpnorm(_lpnorm),
		J_p(_J_p),
		p(_p),
		ScalarVectorProduct(_ScalarVectorProduct)
{}

const std::pair<double, double>
VectorProjection::operator()(
		const Eigen::VectorXd &_projectedonto,
		const Eigen::VectorXd &_tobeprojected,
		const double _Tol) const
{
	BregmanDistance distance(
			lpnorm,
			J_p,
			p,
			ScalarVectorProduct
			);


	double tmin = 0.;
	double result = 0.;
	{
		// tmin=fminunc(@(t) BregmanFunctional(t,Jx,u,alpha+d,DualNormX,DualPowerX,TolX),t0,BregmanOptions);
		BregmanDistanceToLine distancefunctional(
				distance,
				lpnorm,
				J_p,
				_projectedonto,
				_tobeprojected,
				p);

		FunctionMinimizer<double> minimizer(
				distancefunctional, tmin);

//		const unsigned int inner_iterations =
				minimizer(1, _Tol, tmin);
		result = distancefunctional(tmin);

		BOOST_LOG_TRIVIAL(trace)
			<< "tmin is " << tmin;
	}

	return std::make_pair( result, tmin);
}

// instantiate template functions
CONSTRUCT_FUNCTIONMINIMIZER(double)
