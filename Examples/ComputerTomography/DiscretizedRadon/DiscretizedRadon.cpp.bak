/*
 * DiscretizedRadon.cpp
 *
 *  Created on: Jul 9, 2015
 *      Author: heber
 */


#include "BassoConfig.h"

#include "DiscretizedRadon.hpp"

#include "Minimizations/Elements/ElementCreator.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Elements/RepresentationAdvocate.hpp"

DiscretizedRadon::DiscretizedRadon(
		const NormedSpace_weakptr_t _SourceSpaceRef,
		const NormedSpace_weakptr_t _TargetSpaceRef,
		const unsigned int _num_pixel_x,
		const unsigned int _num_pixel_y,
		const unsigned int _num_angles,
		const unsigned int _num_offsets) :
	Mapping(_SourceSpaceRef, _TargetSpaceRef),
	radon_matrix(_num_pixel_x, _num_pixel_y, _num_angles, _num_offsets),
	MatrixVectorProductCounts(0)
{}

void DiscretizedRadon::operator()(
		const SpaceElement_ptr_t &_sourceelement,
		SpaceElement_ptr_t &_destelement
		) const
{
	assert( _sourceelement->getSpace() == getSourceSpace() );
	++MatrixVectorProductCounts;

	const boost::chrono::high_resolution_clock::time_point timing_start =
			boost::chrono::high_resolution_clock::now();
	Eigen::VectorXd tempvector;
	tempvector = radon_matrix.getMatrix() * RepresentationAdvocate::get(_sourceelement);
	const boost::chrono::high_resolution_clock::time_point timing_end =
			boost::chrono::high_resolution_clock::now();
	MatrixVectorProductTimings += timing_end - timing_start;

	RepresentationAdvocate::set(_destelement, tempvector);
}

const Mapping_ptr_t DiscretizedRadon::getAdjointMapping() const
{
	// this throws if AdjointLinearMapping is expired
	return Mapping_ptr_t(AdjointLinearMapping);
}

void DiscretizedRadon::setAdjointMapping(const Mapping_weakptr_t &_adjoint)
{
	assert( AdjointLinearMapping.expired() );
	const_cast<Mapping_weakptr_t &>(AdjointLinearMapping) = _adjoint;
}
