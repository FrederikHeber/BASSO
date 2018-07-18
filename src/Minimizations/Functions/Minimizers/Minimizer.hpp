/*
 * Minimizer.hpp
 *
 *  Created on: Nov 4, 2014
 *      Author: heber
 */

#ifndef MINIMIZER_HPP_
#define MINIMIZER_HPP_

#include "BassoConfig.h"

#include <boost/function.hpp>
#include <vector>

#include "Minimizations/Functions/Minimizers/FunctionMinimizer.hpp"

template <class T>
class Minimizer : public FunctionMinimizer
{
	typedef std::vector<double> array_type;

public:
	Minimizer(const unsigned int _N);
	virtual ~Minimizer();

	enum FunctionMinimizer::GradientStatus
	checkGradient(const double _tol) const;

	//!> typedef for the function to check when to stop the iteration
	typedef boost::function<
			enum FunctionMinimizer::GradientStatus (const double) > check_function_t;

	const unsigned int operator()(
			const double _Tol,
			array_type &_startvalue,
			const check_function_t &_checkfunction
			);

	void setMaxIterations(const unsigned int _maxiterations);

	const double getCurrentOptimumValue() const;

	const array_type &getCurrentOptimum() const;

	const array_type &getCurrentGradient() const;

	void convertToInternalType(
			const array_type & _t,
			T * const _x) const;

	void convertFromInternalType(
			const T * const _x,
			array_type &_t) const;
};

#include "Minimizer_gsl.hpp"
#include "Minimizer_nlopt.hpp"

#endif /* MINIMIZER_HPP_ */
