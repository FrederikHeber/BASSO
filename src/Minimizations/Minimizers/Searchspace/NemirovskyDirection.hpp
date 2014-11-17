/*
 * NemirovskyDirection.hpp
 *
 *  Created on: Nov 17, 2014
 *      Author: heber
 */

#ifndef MINIMIZATIONS_MINIMIZERS_SEARCHSPACE_NEMIROVSKYDIRECTION_HPP_
#define MINIMIZATIONS_MINIMIZERS_SEARCHSPACE_NEMIROVSKYDIRECTION_HPP_

#include <Minimizations/Minimizers/Searchspace/SearchspaceUpdater.hpp>
#include "BassoConfig.h"

#include "Searchspace.hpp"

class NemirovskyDirection : public Searchspace
{
	/** This function performs the actual update of the search space.
	 *
	 * @param _iterate current dual iterate
	 * @param _newdir current search direction
	 * @param _alpha current alpha to this \a _newdir
	 */
	void update(
			const SpaceElement_ptr_t &_iterate,
			const SpaceElement_ptr_t &_newdir,
			const double _alpha);

};



#endif /* MINIMIZATIONS_MINIMIZERS_SEARCHSPACE_NEMIROVSKYDIRECTION_HPP_ */
