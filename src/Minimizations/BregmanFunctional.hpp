/*
 * BregmanFunctional.hpp
 *
 *  Created on: Mar 27, 2014
 *      Author: heber
 */

#ifndef BREGMANFUNCTIONAL_HPP_
#define BREGMANFUNCTIONAL_HPP_

#include <Eigen/Dense>
#include <utility>

#include "Minimizations/DualityMapping.hpp"

/** Functor to calculate BregmanFunctional functional/distance.
 *
 */
template <unsigned int p>
class BregmanFunctional
{
public:
	BregmanFunctional() {}
	~BregmanFunctional() {}

	/** Implements BregmanFunctional functional.
	 *
	 * \param _t
	 * \param _x vector
	 * \param _U
	 * \param _alpha
	 * \param _q power of the weight of the duality mapping
	 */
	std::pair<double,Eigen::VectorXd> operator()(
			const Eigen::VectorXd &_t,
			const Eigen::VectorXd &_x,
			const Eigen::MatrixXd &_U,
			const Eigen::VectorXd &_alpha,
			const unsigned int _q
			)
	{
		// x=x-U*t;
		const Eigen::VectorXd resx = _x - _U * _t;
		// fval=1/q*norm(x,p)^q+alpha'*t;
		const double fval =
				1./(double)_q * ::pow(resx.lpNorm<p>(), _q)
				+ _alpha.transpose() * _t;
		// gval=alpha-U'*DualityMapping(x,p,q,Tol);
		const DualityMapping<p> J_p(_q);
		const Eigen::VectorXd gval =
				_alpha -
				_U.transpose() * J_p(resx);

		return std::make_pair( fval, gval );
	}

};


#endif /* BREGMANFUNCTIONAL_HPP_ */
