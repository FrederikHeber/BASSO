/*
 * RegularizedL1DualityMapping.cpp
 *
 *  Created on: May 01, 2017
 *      Author: heber
 */

#include "BassoConfig.h"

#include "RegularizedL1DualityMapping.hpp"

#include <cmath>

#include "Minimizations/Elements/ElementCreator.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Elements/RepresentationAdvocate.hpp"
#include "Minimizations/Mappings/Specifics/RelativeShrinkageMapping.hpp"
#include "Minimizations/Norms/L1Norm.hpp"
#include "Math/Helpers.hpp"

/** General function to calculate the duality mapping.
 *
 * This is a single-valued selection for the duality mapping to the non-smooth
 * regularized l1 Banach space.
 *
 * \param _x vector
 * \param _Jx duality mapped \a _x
 */
void RegularizedL1DualityMapping::operator()(
		const SpaceElement_ptr_t &_x,
		SpaceElement_ptr_t &_Jx
		) const
{
	// start timing
	const boost::chrono::high_resolution_clock::time_point timing_start =
			boost::chrono::high_resolution_clock::now();

	// single-valued selection for l1 part
	// J=norm(x,1)^(q-1)*sign(x);
	assert( getSourceSpace().get() == _x->getSpace().get() );
	assert( getTargetSpace().get() == _Jx->getSpace().get() );
	const L1Norm &l1norm = static_cast<L1Norm &>(*_x->getSpace()->getNorm());

	const double factor = ::pow(l1norm.L1Norm::operator()(_x), (double)power-1.);
	assert( _Jx->getSpace()->getDimension() == _x->getSpace()->getDimension() );
	SpaceElement_ptr_t _Jx1 = getTargetSpace()->createElement();
	const Eigen::VectorXd &vector = RepresentationAdvocate::get(_x);
	RepresentationAdvocate::set(_Jx1,
			factor * vector.array().cwiseProduct(
					vector.array().abs().cwiseInverse()));
	SpaceElement_ptr_t _Jx2 = getTargetSpace()->createElement();
	RepresentationAdvocate::set(_Jx1, vector);
	_Jx = _Jx1 + lambda * _Jx2;

	// finish timing
	const boost::chrono::high_resolution_clock::time_point timing_end =
			boost::chrono::high_resolution_clock::now();
	timing += timing_end - timing_start;
	++count;
}

void RegularizedL1DualityMapping::getMinimumInfimum(
		const SpaceElement_ptr_t &_x,
		const SpaceElement_ptr_t &_y,
		SpaceElement_ptr_t &_Jx) const
{
	// start timing
	const boost::chrono::high_resolution_clock::time_point timing_start =
			boost::chrono::high_resolution_clock::now();

	// l1 is (component-wise) only not differentiable where we intersect with axis
	// hence, for all intersections check whether flipped sign is better
	assert( getSourceSpace().get() == _x->getSpace().get() );
	assert( getTargetSpace().get() == _Jx->getSpace().get() );
	const L1Norm &l1norm = static_cast<L1Norm &>(*_x->getSpace()->getNorm());

	const double factor = ::pow(l1norm.L1Norm::operator()(_x), (double)power-1.);
	assert( _Jx->getSpace()->getDimension() == _x->getSpace()->getDimension() );
	SpaceElement_ptr_t _Jx1 = getTargetSpace()->createElement();
	const Eigen::VectorXd &vector = RepresentationAdvocate::get(_x);
	RepresentationAdvocate::set(_Jx1,
			factor * vector.array().cwiseProduct(
					vector.array().abs().cwiseInverse()));
	for (unsigned int i=0;i<_x->getSpace()->getDimension();++i) {
		if (fabs((*_x)[i]) < BASSOTOLERANCE ) {
			if ((*_Jx1)[i] * (*_y)[i] > 0)
				(*_Jx1)[i] *= -1.;
		}
	}
	SpaceElement_ptr_t _Jx2 = getTargetSpace()->createElement();
	RepresentationAdvocate::set(_Jx1, vector);
	_Jx = _Jx1 + lambda * _Jx2;

	// finish timing
	const boost::chrono::high_resolution_clock::time_point timing_end =
			boost::chrono::high_resolution_clock::now();
	timing += timing_end - timing_start;
	++count;
}

const Mapping_ptr_t RegularizedL1DualityMapping::getAdjointMapping() const
{
	Mapping_ptr_t instance(
			new RelativeShrinkageMapping(
					getTargetSpace(), lambda)
	);
	return instance;
}
