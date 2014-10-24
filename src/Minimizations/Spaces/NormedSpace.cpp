/*
 * NormedSpace.cpp
 *
 *  Created on: Oct 23, 2014
 *      Author: heber
 */


#include "BassoConfig.h"

#include "NormedSpace.hpp"

#include "Minimizations/Elements/SpaceElement.hpp"

SpaceElement_ptr_t NormedSpace::createElement() const
{
	SpaceElement_ptr_t newelement(new SpaceElement(getSpace()));
	newelement->setSelfRef(newelement);
	return newelement;
}
