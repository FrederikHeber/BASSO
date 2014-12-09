/*
 * NemirovskyDirection.cpp
 *
 *  Created on: Nov 17, 2014
 *      Author: heber
 */


#include "BassoConfig.h"

#include "NemirovskyDirection.hpp"

#include <cassert>

#include "Log/Logging.hpp"

#include "Minimizations/Elements/SpaceElement.hpp"

NemirovskyDirection::NemirovskyDirection(
		const NormedSpace_ptr_t &_SearchDirectionSpace_ptr) :
		Searchspace(_SearchDirectionSpace_ptr, 2)
{}

void NemirovskyDirection::update(
		const SpaceElement_ptr_t &_newdir,
		const double _alpha,
		const SpaceElement_ptr_t &_dual_iterate,
		const SpaceElement_ptr_t &_iterate
		)
{
	assert( _newdir->getSpace() == SearchDirectionSpace_ptr );
	assert( _dual_iterate->getSpace() == SearchDirectionSpace_ptr );

	// update search direction
	*(U[0]) = _newdir;
	alphas[0] = _alpha;

	// update Nemirovsky direction
	*(U[1]) = _dual_iterate;
	alphas[1] = _dual_iterate * _iterate;
}
