/*
 * Minimizer_nlopt.cpp
 *
 *  Created on: Nov 6, 2014
 *      Author: heber
 */


#include "BassoConfig.h"

#include "Minimizer_nlopt.hpp"

#include <sstream>

#include "Log/Logging.hpp"

Minimizer<NLopt_vector>::Minimizer(
		const unsigned int _N
		) :
		opt(nlopt::LD_MMA, _N),
		N(_N),
		iter(0),
		maxiterations(100),
		constant_positivity(0.),
		tolerance(1e-4),
		optimum(0.),
		tempoptimum(N),
		tempgradient(N)
{
	// prepare objective function
	opt.set_min_objective(&FunctionGradientCaller,
			static_cast<void *>(this));
}

enum Minimization::GradientStatus
Minimizer<NLopt_vector>::checkGradient(
		const double _Tol) const
{
	// we do nothing as nlopt::opt::set_xtol_rel does all the work instead
	if (iter >= maxiterations)
		return Minimization::gradient_success;
	else
		return Minimization::gradient_continue;
}

const unsigned int
Minimizer<NLopt_vector>::minimize(
		const double _Tol,
		array_type &_startvalue,
		const check_function_t &_checkfunction
		)
{
	checkfunction = _checkfunction;
	tolerance = _Tol;
	iter = 0;

	// set lower bound positivity constraint
	std::vector<double> lb(N, -HUGE_VAL);
	for (std::vector<unsigned int>::const_iterator it = PositivityBoundIndices.begin();
			it != PositivityBoundIndices.end(); ++it) {
		lb[*it] = constant_positivity;
		// make sure _startvalue adheres bounds
		if (_startvalue[*it] < constant_positivity)
			_startvalue[*it] = constant_positivity;
	}
	opt.set_lower_bounds(lb);

	// set tolerance
	opt.set_xtol_rel(_Tol);

	double minf;
	try {
//		nlopt::result result =
				opt.optimize(_startvalue, minf);
	} catch (nlopt::forced_stop &e) {
		BOOST_LOG_TRIVIAL(debug)
				<< "NLopt stopped forcedly by checkfunction, "
				<< "checkfunction says: "
				<< (_checkfunction(_Tol) == Minimization::gradient_success ?
						"valid" : "INVALID");
	}

	BOOST_LOG_TRIVIAL(debug)
		<< "Inner iteration took " << iter << " steps";

	return iter;
}

double
Minimizer<NLopt_vector>::FunctionGradientCaller(
		const std::vector<double> &x,	/* x */
		std::vector<double> &g, 	/* grad */
		void *adata
		)
{
	// obtain minimizer object
	struct Minimizer<NLopt_vector> *minimizer =
			static_cast<Minimizer<NLopt_vector> *>(adata);

	// evaluate the function and store in internal values
	minimizer->tempoptimum = x;
	minimizer->optimum = minimizer->function_evaluator(x);
	minimizer->tempgradient = minimizer->gradient_evaluator(x);
	g = minimizer->tempgradient;
	++(minimizer->iter);

	// check stop conditions
	if (minimizer->iter >= minimizer->maxiterations)
		throw nlopt::forced_stop();
	if (minimizer->checkfunction(minimizer->tolerance)
			!= Minimization::gradient_continue)
		throw nlopt::forced_stop();

	// return
	return minimizer->optimum;
}
