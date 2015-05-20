/*
 * HyperplaneProjection.hpp
 *
 *  Created on: Sep 4, 2014
 *      Author: heber
 */

#ifndef HYPERPLANEPROJECTION_HPP_
#define HYPERPLANEPROJECTION_HPP_

#include "BassoConfig.h"

#include <vector>

#include "Log/Logging.hpp"

#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Functions/BregmanProjectionFunctional.hpp"
#include "Minimizations/Functions/Minimizers/MinimizationFunctional.hpp"
#include "Minimizations/types.hpp"

/** Structure containing all parameters to call ..ProjectionFunctional's functions.
 *
 * This is required to use function minimization that only allows
 * to pass a void* pointer to pass on information to the function to be
 * minimized.
 *
 * It calculates the projection onto the intersection of hyperplanes.
 *
 * \sa BregmanProjectionFunctional, MetricProjectionFunctional
 *
 */
template <class F>
struct HyperplaneProjection :
		public MinimizationFunctional< std::vector<double> >
{
	typedef typename MinimizationFunctional< std::vector<double> >::array_type array_type;

	const F &functional;
	const SpaceElement_ptr_t &x;
	const std::vector<SpaceElement_ptr_t> &U;
	const std::vector<double> &alpha;

	/** Constructor to initialize refs.
	 *
	 */
	HyperplaneProjection(
		const F &_functional,
		const SpaceElement_ptr_t &_x,
		const std::vector<SpaceElement_ptr_t> &_U,
		const std::vector<double> &_alpha
		) :
			functional(_functional),
			x(_x),
			U(_U),
			alpha(_alpha)
	{}

	double function(
			const std::vector<double> &_value
			) const
	{
		const double returnvalue =
				functional(
						_value,
						x,
						U,
						alpha);
		BOOST_LOG_TRIVIAL(trace)
			<< "func() evaluates to " << returnvalue;
		return returnvalue;
	}

	const std::vector<double> gradient(
			const std::vector<double> &_value
			) const
		{
			const std::vector<double> grad =
					functional.gradient(
							_value,
							x,
							U,
							alpha);
			std::stringstream outputstream;
			outputstream << grad;
			BOOST_LOG_TRIVIAL(trace)
				<< "grad() evaluates to " << outputstream.str();
			return grad;
		}

	void convertInternalTypeToArrayType(
			const std::vector<double> &_t,
			array_type & _x
			) const
	{
		_x = _t;
	}

	void convertArrayTypeToInternalType(
			const array_type & _x,
			std::vector<double> &_t
			) const
	{
		_t = _x;
	}
};



#endif /* HYPERPLANEPROJECTION_HPP_ */
