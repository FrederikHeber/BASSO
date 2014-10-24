/*
 * SmoothnessFunctional.hpp
 *
 *  Created on: Oct 24, 2014
 *      Author: heber
 */

#ifndef SMOOTHNESSFUNCTIONAL_HPP_
#define SMOOTHNESSFUNCTIONAL_HPP_

#include "BassoConfig.h"

class SmoothnessModulus;

/** This class is a functor evaluating a given smoothness modulus
 * in order to match with a given lambda.
 */
class SmoothnessFunctional
{
public:
	SmoothnessFunctional(
			const SmoothnessModulus &_modul,
			const double _lambda
			) :
		modul(_modul),
		lambda(_lambda)
	{}

	double operator()(double _arg) const
	{
		const double result = (modul)(_arg);
		const double norm = result/_arg - lambda;
		return norm*norm;
	}

private:
	const SmoothnessModulus &modul;
	const double lambda;
};




#endif /* SMOOTHNESSFUNCTIONAL_HPP_ */
