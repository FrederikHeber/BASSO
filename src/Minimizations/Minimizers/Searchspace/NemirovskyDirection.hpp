/*
 * NemirovskyDirection.hpp
 *
 *  Created on: Nov 17, 2014
 *      Author: heber
 */

#ifndef MINIMIZATIONS_MINIMIZERS_SEARCHSPACE_NEMIROVSKYDIRECTION_HPP_
#define MINIMIZATIONS_MINIMIZERS_SEARCHSPACE_NEMIROVSKYDIRECTION_HPP_

#include "BassoConfig.h"

#include "Searchspace.hpp"

class NemirovskyDirection : public Searchspace
{
public:
	/** Constructor for class NemirovskyDirection.
	 *
	 * We have here two search directions, the current one
	 * and the current dual iterate.
	 *
	 * @param _SearchDirectionSpace_ptr search direction space (for checks)
	 * products in the subspace (i.e. search space)
	 */
	NemirovskyDirection(
			const NormedSpace_ptr_t &_SearchDirectionSpace_ptr);

	/** This function performs the actual update of the search space.
	 *
	 * @param _newdir current search direction
	 * @param _alpha current alpha to this \a _newdir
	 * @param _dual_iterate current dual iterate
	 * @param _iterate current iterate
	 */
	void update(
			const SpaceElement_ptr_t &_newdir,
			const double _alpha,
			const SpaceElement_ptr_t &_dual_iterate,
			const SpaceElement_ptr_t &_iterate
			);

	/** Returns 0, as is always the search direction.
	 *
	 * @return 0
	 */
	const unsigned int getIndex() const
	{ return 0; }

	/** Returns 0, as is always the search direction.
	 *
	 * @return vector(1,0)
	 */
	const std::vector<unsigned int> getLastIndices() const
	{ return std::vector<unsigned int>(1,0); }
};



#endif /* MINIMIZATIONS_MINIMIZERS_SEARCHSPACE_NEMIROVSKYDIRECTION_HPP_ */
