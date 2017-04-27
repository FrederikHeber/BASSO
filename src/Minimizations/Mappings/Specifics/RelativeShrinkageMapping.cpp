/*
 * RelativeShrinkageMapping.cpp
 *
 *  Created on: Oct 13, 2014
 *      Author: heber
 */

#include "RelativeShrinkageMapping.hpp"

#include "Log/Logging.hpp"
#include "Math/Helpers.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Mappings/IllegalDualityMapping.hpp"
#include "Minimizations/Minimizers/MinimizationExceptions.hpp"
#include "Minimizations/Norms/L1Norm.hpp"
#include "Minimizations/Norms/Specifics/RelativeShrinkageCoefficient.hpp"

RelativeShrinkageMapping::RelativeShrinkageMapping(
		const NormedSpace_weakptr_t &_NormedSpaceRef,
		const double _lambda) :
	L1DualityMapping(_NormedSpaceRef, 2.),
	lambda(_lambda),
	count(0),
	timing(boost::chrono::nanoseconds(0))
{
	if (!(lambda > 0))
		throw MinimizationIllegalValue_exception()
		<< MinimizationIllegalValue_name("lambda");
}

template <class T>
struct applyShrinkrage : public std::binary_function<T,T,double>
{
	double operator()(const T _value, const double _coefficient) const
	{
		const T abs_value = fabs(_value);
		return (abs_value < _coefficient ?
				0. :
				(abs_value-_coefficient)*Helpers::sign(_value));
	}
};

const double RelativeShrinkageMapping::getRelativeShrinkage(const SpaceElement_ptr_t &_x) const
{
	return RelativeShrinkageCoefficient::get(_x, lambda);
}

void RelativeShrinkageMapping::operator()(
		const SpaceElement_ptr_t &_x,
		SpaceElement_ptr_t &_Jx) const
{
	// start timing
	const boost::chrono::high_resolution_clock::time_point timing_start =
			boost::chrono::high_resolution_clock::now();

	const double coefficient = getRelativeShrinkage(_x);
	RepresentationAdvocate::set(
			_Jx,
			(1./lambda)*RepresentationAdvocate::get(_x).unaryExpr(
					std::bind2nd(applyShrinkrage<double>(), coefficient))
	);

	// finish timing
	const boost::chrono::high_resolution_clock::time_point timing_end =
			boost::chrono::high_resolution_clock::now();
	timing += timing_end - timing_start;
	++count;
}

const Mapping_ptr_t RelativeShrinkageMapping::getAdjointMapping() const
{
	Mapping_ptr_t instance(
			new IllegalDualityMapping
	);
	return instance;
}

