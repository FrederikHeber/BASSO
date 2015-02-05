/*
 * InverseProblem.cpp
 *
 *  Created on: Oct 27, 2014
 *      Author: heber
 */


#include "BassoConfig.h"

#include "InverseProblem.hpp"

#include <cassert>

#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Mappings/Mapping.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"

InverseProblem::InverseProblem(
		const Mapping_ptr_t &_A,
		const NormedSpace_ptr_t &_SourceSpace,
		const NormedSpace_ptr_t &_TargetSpace,
		const SpaceElement_ptr_t &_y
		) :
	SourceSpace(_SourceSpace),
	DualSourceSpace(_SourceSpace->getDualSpace()),
	TargetSpace(_TargetSpace),
	DualTargetSpace(_TargetSpace->getDualSpace()),
	A(_A),
	A_t(_A->getAdjointMapping()),
	x(_SourceSpace->createElement()),
	y(_y)
{
	assert( A->getTargetSpace() == y->getSpace() );
}

