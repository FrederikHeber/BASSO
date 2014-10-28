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
		const SpaceElement_ptr_t &_y
		) :
	A(_A),
	y(_y) /* ,
	x( A->getSourceSpace()->createElement() ) */
{
//	assert( A->getTargetSpace() == y->getSpace() );
}

