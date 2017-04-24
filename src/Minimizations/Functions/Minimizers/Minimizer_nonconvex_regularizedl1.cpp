/*
 * Minimizer_nonconvex_regularizedl1.cpp
 *
 *  Created on: Jan 2, 2017
 *      Author: heber
 */


#include "BassoConfig.h"

#include "Minimizer_nlopt.hpp"

//#include <fstream>
#include <sstream>

#include <boost/bind.hpp>
#include <boost/math/tools/minima.hpp>

#include "Log/Logging.hpp"
#include "Minimizations/Functions/Minimizers/MinimizerExceptions.hpp"

Minimizer<NonConvexRegL1>::Minimizer(
		const unsigned int _N,
		const std::vector< std::vector<double> > &_sections_per_direction
		) :
		FunctionMinimizer(_N),
#ifdef NLOPT_FOUND
		opt(nlopt::LD_MMA, _N),
#endif /* NLOPT_FOUND */
		N(_N),
		iter(0),
		constant_positivity(0.),
		tolerance(1e-4),
		sections_per_direction(_sections_per_direction)
{
	// prepare objective function
#ifdef NLOPT_FOUND
	opt.set_min_objective(&FunctionGradientCaller,
			static_cast<void *>(this));
#endif /* NLOPT_FOUND */
}

enum FunctionMinimizer::GradientStatus
Minimizer<NonConvexRegL1>::checkGradient(
		const double _Tol) const
{
	// we do nothing as nlopt::opt::set_xtol_rel does all the work instead
	if (iter >= maxiterations)
		return FunctionMinimizer::gradient_success;
	else
		return FunctionMinimizer::gradient_continue;
}

static bool checkIteratorVectorAtEnd(
		const Minimizer<NonConvexRegL1>::sections_per_direction_t &_sections_per_direction,
		const std::vector< Minimizer<NonConvexRegL1>::sections_t::const_iterator > &_iter_vector)
{
	assert(_iter_vector.size() == _sections_per_direction.size());
	bool status = true;
	for (size_t i=0;i<_iter_vector.size();++i)
		if (_iter_vector[i] != _sections_per_direction[i].end())
			status = false;
	return status;
}

static size_t advanceIteratorVector(
		const Minimizer<NonConvexRegL1>::sections_per_direction_t &_sections_per_direction,
		std::vector< Minimizer<NonConvexRegL1>::sections_t::const_iterator > &_iter_vector)
{
	assert(_iter_vector.size() == _sections_per_direction.size());
	size_t i=0;
	for (;i<_iter_vector.size();++i)
		if (_iter_vector[i] != _sections_per_direction[i].end()) {
			// skip equivalent interval boundaries
			const double oldvalue = *_iter_vector[i];
			while ((fabs(oldvalue - *_iter_vector[i]) < BASSOTOLERANCE)
					&& (_iter_vector[i] != _sections_per_direction[i].end()))
				++(_iter_vector[i]);
			// set all previous ones to begin again
			for (size_t j=0;j<i;++j)
				_iter_vector[j] = _sections_per_direction[j].begin();
			break;
		}
	assert( (i >= 0) && (i<_iter_vector.size()) );
	return i;
}

static double FunctionCall(
		void *adata,
		const double &_value)
{
	FunctionMinimizer::array_type x(1,_value);
	FunctionMinimizer::array_type g(1,0.);
	const double fval = Minimizer<NonConvexRegL1>::FunctionGradientCaller(x,g,adata);
//	LOG(debug, "F(" << x[0] << ") = " << fval);
	return fval;
}

const unsigned int
Minimizer<NonConvexRegL1>::minimize(
		const double _Tol,
		array_type &_startvalue,
		const check_function_t &_checkfunction
		)
{
	checkfunction = _checkfunction;
	tolerance = _Tol;

//	/// for debugging we keep a map of values to visualize the functional
//	typedef double functionvalue_t;
//	typedef std::map< array_type, functionvalue_t > functionvalues_t;
//	functionvalues_t functionvalues;

	// here, we store the lowest function value and the associated parameters
	array_type startvalue(_startvalue);
	array_type bestvalue(_startvalue);
	double lowestf = std::numeric_limits<double>::max();

	// to avoid recursion, we use a vector of iterators which we loop through
	std::vector< sections_t::const_iterator > iter_vector;
	for (sections_per_direction_t::const_iterator diriter = sections_per_direction.begin();
			diriter != sections_per_direction.end(); ++diriter)
		iter_vector.push_back(diriter->begin());

	iter = 0;
	// we have to go through the complete combinatorics of all possible combination
	// of intervals
	std::vector<double> lb(N, -HUGE_VAL);
	std::vector<double> ub(N, +HUGE_VAL);
	bool continue_flag = true;
	while (continue_flag) {
		// we need continue_flag because also the combination with all iterators
		// at end is a legal one
		continue_flag = !checkIteratorVectorAtEnd(sections_per_direction, iter_vector);
		// prepare bounds
		for (size_t i=0;i<N;++i) {
			if (iter_vector[i] == sections_per_direction[i].begin()) {
				lb[i] = -HUGE_VAL;
				ub[i] = *(iter_vector[i]);
			} else {
				sections_t::const_iterator previousiter = iter_vector[i]-1;
				lb[i] = *previousiter;
				if (iter_vector[i] == sections_per_direction[i].end())
					ub[i] = HUGE_VAL;
				else
					ub[i] = *(iter_vector[i]);
			}
			// enforce starting value to be inside interval
			if ((startvalue[i] < lb[i]) || (startvalue[i] > ub[i])) {
				if (lb[i] == -HUGE_VAL)
					_startvalue[i] = ub[i]-1.;
				else if (ub[i] == HUGE_VAL)
					_startvalue[i] = lb[i]+1.;
				else
					_startvalue[i] = (lb[i]+ub[i])/2.;
			} else
				_startvalue[i] = startvalue[i];
//			LOG(trace, "Intervals for dir #" << N << ": [" << lb[i] << ","
//					<< ub[i] << "], starting at " << _startvalue[i]);
		}

//		/// for debugging we sample the function in the interval [lb, ub]
//		array_type coordinates(N,0.);
//		array_type delta(N,0.);
//		array_type start(N,0.);
//		array_type end(N,0.);
//		const int samplepoints = 40;
//		for (int i=0;i<N;++i) {
//			start[i] = (lb[i] < -1000.) ? -1000. : lb[i];
//			end[i] = (ub[i] > 1000.) ? 1000. : ub[i];
//			delta[i] = (end[i] - start[i])/(double)samplepoints;
//			coordinates[i] = start[i];
//		}
//		while (coordinates < end) {
//			const double fval = function_evaluator(coordinates);
//			// we don't care when entry is already present
//			functionvalues.insert( std::make_pair(coordinates, fval));
//			int i=0;
//			for (;i<N;++i)
//				if (coordinates[i] != end[i]) {
//					coordinates[i] += delta[i];
//					break;
//				}
//			if (i != N)
//				for (int j=0;j<i;++j)
//					coordinates[j] = start[i];
//		}


		double minf;
		if(1) {
#ifdef NLOPT_FOUND
		nlopt::opt opt_local(nlopt::LD_MMA, N);
		opt_local.set_min_objective(&FunctionGradientCaller,
				static_cast<void *>(this));
		opt_local.set_lower_bounds(lb);
		opt_local.set_upper_bounds(ub);

		// set tolerance
		opt_local.set_xtol_rel(_Tol);

		nlopt::result result;
		try {
			result =
					opt_local.optimize(_startvalue, minf);
		} catch (nlopt::forced_stop &e) {
			LOG(debug,"NLopt stopped forcedly by checkfunction, "
					<< "checkfunction says: " << (_checkfunction(_Tol) == FunctionMinimizer::gradient_success ?
							"valid" : "INVALID"));
		}
		if (isnan(optimum) || isinf(optimum))
			throw MinimizerIllegalNumber_exception()
			<< MinimizerIllegalNumber_variablename("optimum");
		if (isnan(minf) || isinf(minf)) {
			const double fval = function_evaluator(_startvalue);
			std::stringstream output;
			std::copy(_startvalue.begin(), _startvalue.end(),
					std::ostream_iterator<double>(output, ","));
			LOG(debug, "minf is illegal, checking F(" << output.str() << ") = " << fval);
			if (isnan(fval) || isinf(fval))
				throw MinimizerIllegalNumber_exception()
				<< MinimizerIllegalNumber_variablename("minf");
			else
				minf = fval;
		}
		for (size_t i=0;i<N;++i) {
			if ((_startvalue[i] < lb[i]) && (_startvalue[i] > ub[i])) {
				std::stringstream name("x[");
				name << i << "]";
				throw MinimizerIllegalNumber_exception()
				<< MinimizerIllegalNumber_variablename(name.str());
			}
		}
#else
		LOG(error, "NLopt support not compiled in.");
		assert(0);
#endif /* NLOPT_FOUND */
		} else {
			if (lb[0] == -HUGE_VAL)
				lb[0] = ub[0]-1000.;
			if (ub[0] == HUGE_VAL)
				ub[0] = lb[0]+1000.;
			assert( N==1 );
			int bits = 32;
			boost::uintmax_t maxiter = 1000;
			std::pair<double, double> minpair =
					boost::math::tools::brent_find_minima(
							boost::bind(&FunctionCall,
									static_cast<void *>(this),
									_1),
							lb[0],
							ub[0],
							bits,
							maxiter);
			_startvalue[0] = minpair.first;
			minf = minpair.second;
		}

//		/// for debugging: add found minimum
//		functionvalues.insert( std::make_pair(_startvalue, minf));

		// check and pick best
//		LOG(trace, "Minimum value is " << minf);
		if ((minf < lowestf) && (_startvalue[0] > 0.)) {
			lowestf = minf;
			bestvalue = _startvalue;
			std::stringstream output;
			std::copy(bestvalue.begin(), bestvalue.end(),
					std::ostream_iterator<double>(output, ","));
			LOG(debug, "New best minimum value is " << lowestf << " with " << output.str());
		}

		// next combination
		if (continue_flag)
			advanceIteratorVector(sections_per_direction, iter_vector);
	}
	LOG(debug, "Inner iteration took " << iter << " steps");
	_startvalue = bestvalue;

//	/// for debugging: write function values to temporary file
//	{
//		std::ofstream file_output("test.csv");
//		for (int i=0;i<N;++i)
//			file_output << "x" << i+1 << "\t";
//		file_output << "value\n";
//		for (functionvalues_t::const_iterator iter = functionvalues.begin();
//				iter != functionvalues.end(); ++iter) {
//			std::copy(iter->first.begin(), iter->first.end(),
//					std::ostream_iterator<double>(file_output, "\t"));
//			file_output << iter->second << "\n";
//		}
//	}

	return iter;
}

double
Minimizer<NonConvexRegL1>::FunctionGradientCaller(
		const std::vector<double> &x,	/* x */
		std::vector<double> &g, 	/* grad */
		void *adata
		)
{
	// obtain minimizer object
	struct Minimizer<NonConvexRegL1> *minimizer =
			static_cast<Minimizer<NonConvexRegL1> *>(adata);

	// evaluate the function and store in internal values
	minimizer->tempoptimum = x;
	minimizer->optimum = minimizer->function_evaluator(x);
	minimizer->tempgradient = minimizer->gradient_evaluator(x);
	g = minimizer->tempgradient;
	++(minimizer->iter);

	// check stop conditions
//#ifdef NLOPT_FOUND
//	if (minimizer->iter >= minimizer->maxiterations)
//		throw nlopt::forced_stop();
//	if (minimizer->checkfunction(minimizer->tolerance)
//			!= FunctionMinimizer::gradient_continue)
//		throw nlopt::forced_stop();
//#endif /* NLOPT_FOUND */

	// return
	return minimizer->optimum;
}

