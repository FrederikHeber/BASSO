/*
 * FunctionalMinimizerFactory_impl.hpp
 *
 *  Created on: Dec 11, 2015
 *      Author: heber
 */

#ifndef MINIMIZATIONS_FUNCTIONS_MINIMIZERS_FUNCTIONALMINIMIZERFACTORY_IMPL_HPP_
#define MINIMIZATIONS_FUNCTIONS_MINIMIZERS_FUNCTIONALMINIMIZERFACTORY_IMPL_HPP_

#include "BassoConfig.h"

#include "FunctionalMinimizerFactory.hpp"
#include "Minimizations/Functions/Minimizers/FunctionalMinimizer.hpp"
#include "Minimizations/Functions/Minimizers/FunctionalMinimizer_exactLinesearch.hpp"
#include "Minimizations/Functions/Minimizers/FunctionalMinimizer_inexactLinesearch.hpp"
#include "Minimizations/Functions/Minimizers/Minimizer_gsl.hpp"
#include "Minimizations/Functions/Minimizers/Minimizer_nlopt.hpp"

template <class T>
typename FunctionalMinimizer<T>::ptr_t
FunctionalMinimizerFactory::create(
		const unsigned int _N,
		const MinimizationFunctional<T> &_functional,
		T &_value
		)
{
	FunctionMinimizer::ptr_t instance;
	switch (CurrentMinLib) {
	case gnuscientificlibrary:
		instance.reset(new Minimizer<gsl_vector>(_N));
		break;
	case nonlinearoptimization:
		instance.reset(new Minimizer<NLopt_vector>(_N));
		break;
	default:
		BOOST_LOG_TRIVIAL(error)
				<< "Unknown instance desired in FunctionalMinimizerFactory";
		assert(0);
		break;
	}
	instance->setMaxIterations(maxiterations);

	typename FunctionalMinimizer<T>::ptr_t fmin(
			new FunctionalMinimizer_exactLinesearch<T>(
					_functional,
					instance,
					_value));

	return fmin;
}

template <class T>
typename FunctionalMinimizer<T>::ptr_t
FunctionalMinimizerFactory::create(
		const unsigned int _N,
		const MinimizationFunctional<T> &_functional,
		T &_value,
		const double _constant_positivity,
		const std::vector<unsigned int> &_Wolfe_indexset,
		const double _constant_interpolation
		)
{
	FunctionMinimizer::ptr_t instance;
	switch (CurrentMinLib) {
	case gnuscientificlibrary:
		BOOST_LOG_TRIVIAL(error)
				<< "Wolfe line search is not implemented with gsl minimization.";
		assert(0);
		break;
	case nonlinearoptimization:
	{
		Minimizer<NLopt_vector> *minimizer_nlopt =
				new Minimizer<NLopt_vector>(_N);
		minimizer_nlopt->setConstantPositivity(_constant_positivity);
		minimizer_nlopt->setPositivityBoundIndices(_Wolfe_indexset);
		instance.reset(minimizer_nlopt);
		break;
	}
	default:
		BOOST_LOG_TRIVIAL(error)
				<< "Unknown instance desired in FunctionalMinimizerFactory";
		assert(0);
		break;
	}
	instance->setMaxIterations(maxiterations);

	FunctionalMinimizer_inexactLinesearch<T> *fmin_ptr =
			new FunctionalMinimizer_inexactLinesearch<T>(
					_functional,
					instance,
					_value);
	fmin_ptr->setInexactLinesearchParameters(
			_Wolfe_indexset,
			_constant_positivity,
			_constant_interpolation);
	typename FunctionalMinimizer<T>::ptr_t fmin(fmin_ptr);

	return fmin;
}


#endif /* MINIMIZATIONS_FUNCTIONS_MINIMIZERS_FUNCTIONALMINIMIZERFACTORY_IMPL_HPP_ */
