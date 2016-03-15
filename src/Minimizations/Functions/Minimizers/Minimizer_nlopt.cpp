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
#include "Minimizations/Functions/Minimizers/MinimizerExceptions.hpp"

Minimizer<NLopt_vector>::Minimizer(
		const unsigned int _N
		) :
		FunctionMinimizer(_N),
#ifdef NLOPT_FOUND
		opt(nlopt::LD_MMA, _N),
#endif /* NLOPT_FOUND */
		N(_N),
		iter(0),
		constant_positivity(0.),
		tolerance(1e-4)
{
	// prepare objective function
#ifdef NLOPT_FOUND
	opt.set_min_objective(&FunctionGradientCaller,
			static_cast<void *>(this));
#endif /* NLOPT_FOUND */
}

enum FunctionMinimizer::GradientStatus
Minimizer<NLopt_vector>::checkGradient(
		const double _Tol) const
{
	// we do nothing as nlopt::opt::set_xtol_rel does all the work instead
	if (iter >= maxiterations)
		return FunctionMinimizer::gradient_success;
	else
		return FunctionMinimizer::gradient_continue;
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
		lb[*it] = 0.;
		// make sure _startvalue adheres bounds
		if (_startvalue[*it] < constant_positivity)
			_startvalue[*it] = constant_positivity;
	}
#ifdef NLOPT_FOUND
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
				<< (_checkfunction(_Tol) == FunctionMinimizer::gradient_success ?
						"valid" : "INVALID");
	}
	if (isnan(optimum) || isinf(optimum))
		throw MinimizerIllegalNumber_exception()
		<< MinimizerIllegalNumber_variablename("optimum");
	if (isnan(minf) || isinf(minf))
		throw MinimizerIllegalNumber_exception()
		<< MinimizerIllegalNumber_variablename("minf");

	BOOST_LOG_TRIVIAL(debug)
		<< "Inner iteration took " << iter << " steps";
#else
	BOOST_LOG_TRIVIAL(error)
		<< "NLopt support not compiled in.";
	assert(0);
#endif /* NLOPT_FOUND */

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
#ifdef NLOPT_FOUND
	if (minimizer->iter >= minimizer->maxiterations)
		throw nlopt::forced_stop();
	if (minimizer->checkfunction(minimizer->tolerance)
			!= FunctionMinimizer::gradient_continue)
		throw nlopt::forced_stop();
#endif /* NLOPT_FOUND */

	// return
	return minimizer->optimum;
}

