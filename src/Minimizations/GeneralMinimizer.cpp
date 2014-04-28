/*
 * GeneralMinimizer.cpp
 *
 *  Created on: Apr 28, 2014
 *      Author: heber
 */

#include "BassoConfig.h"

#include "Minimizations/GeneralMinimizer.hpp"

#include <boost/log/trivial.hpp>

GeneralMinimizer::GeneralMinimizer(
		const double _NormX,
		const double _NormY,
		const double _PowerX,
		const double _PowerY,
		const double _Delta,
		const double _C,
		const unsigned int _maxiter,
		const unsigned int _outputsteps
		) :
	val_NormX(_NormX),
	val_NormY(_NormY),
	val_DualNormX(val_NormX/(val_NormX - 1.)),
	PowerX(_PowerX),
	PowerY(_PowerY),
	DualPowerX(PowerX/(PowerX - 1.)),
	Delta(_Delta),
	MaxOuterIterations(_maxiter),
	TolX(1e-6),
	TolY(Delta),
	TolFun(1e-12),
	C(_C),
	outputsteps(_outputsteps),
	NormX(val_NormX),
	NormY(val_NormY),
	DualNormX(val_DualNormX),
	J_p(val_NormX),
	J_q(val_DualNormX),
	j_r(val_NormY)
{
	if ((C <= 0.)) // || ( C > 1.))
		throw MinimizationIllegalValue_exception()
			<< MinimizationIllegalValue_name("C");

	BOOST_LOG_TRIVIAL(debug)
		<< "p is " << val_NormX
		<< ", q is " << val_DualNormX
		<< ", r is " << val_NormY
		<< ", power of J_p is " <<  PowerX
		<< ", power of J_q is " <<  DualPowerX
		<< ", power of J_r is " <<  PowerY;

	// set tolerances values
	J_p.setTolerance(TolX);
	J_q.setTolerance(TolX);
	j_r.setTolerance(TolY);
}

