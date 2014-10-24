/*
 * LpNorm.hpp
 *
 *  Created on: Mar 28, 2014
 *      Author: heber
 */

#ifndef LPNORM_HPP_
#define LPNORM_HPP_

#include "BassoConfig.h"

#include <cmath>
#include <Eigen/Dense>

#include "Minimizations/MinimizationExceptions.hpp"

class LpNorm
{
public:
	LpNorm(const double _p) :
		p(_p)
	{
		if (p < 0.)
			throw MinimizationIllegalValue_exception()
				<< MinimizationIllegalValue_name("p");

	}
	~LpNorm() {}

	double operator()(const Eigen::VectorXd &_x) const
	{
		if (p > 0. ) {
			double value = 0.;
			for (unsigned int i=0;i<_x.innerSize();++i)
				value += ::pow(fabs(_x[i]), p);
			return ::pow(value, 1./p);
		} else {
			// infinity norm
			return _x.lpNorm<Eigen::Infinity>();
		}
	}

	enum { Infinity = 0};

private:
	//!> p value for norm
	const double p;
};


#endif /* LPNORM_HPP_ */
