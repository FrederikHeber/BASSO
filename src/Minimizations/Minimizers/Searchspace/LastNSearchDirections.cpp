/*
 * LastNSearchDirections.cpp
 *
 *  Created on: Nov 17, 2014
 *      Author: heber
 */


#include <Minimizations/Minimizers/Searchspace/LastNSearchDirections.hpp>
#include "BassoConfig.h"

#include <algorithm>
#include <Eigen/Dense>
#include <list>

#include "Log/Logging.hpp"

#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Norms/Norm.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"

struct unique_number
{
	unique_number() : number(0)
	{}

	unsigned int operator()()
	{ return number++; }

private:
	unsigned int number;
};

// static entities
LastNSearchDirections::UpdateAlgorithmType
LastNSearchDirections::updateIndexType = LastNSearchDirections::RoundRobin;
bool LastNSearchDirections::enforceRandomMapping = false;

LastNSearchDirections::LastNSearchDirections(
		const NormedSpace_ptr_t &_SearchDirectionSpace_ptr,
		const unsigned int _N) :
	Searchspace(_SearchDirectionSpace_ptr,_N),
	index(0)
{}

void LastNSearchDirections::update(
		const SpaceElement_ptr_t &_newdir,
		const double _alpha,
		const SpaceElement_ptr_t &,
		const SpaceElement_ptr_t &
		)
{
	assert( _newdir->getSpace() == SearchDirectionSpace_ptr );

	switch (updateIndexType) {
		case RoundRobin:
			index = advanceIndex();
			break;
		case MostParallel:
			index = updateIndexToMostParallel(
					*_newdir->getSpace()->getNorm(),
					_newdir);
			break;
		case MostOrthogonal:
			index = updateIndexToMostOrthogonal(
					*_newdir->getSpace()->getNorm(),
					_newdir);
			break;
		default:
			BOOST_LOG_TRIVIAL(error)
				<< "Unknown updateIndex algorithm.";
			assert(0);
	}
	BOOST_LOG_TRIVIAL(debug)
		<< "Updated index is " << index;
	U.col(index) = _newdir->getVectorRepresentation();
	alphas(index) = _alpha;
}

unsigned int
LastNSearchDirections::advanceIndex() const
{
	return ((index + 1) % getDimension());
}

void
LastNSearchDirections::replenishIndexset(
		indexset_t &_indexset) const
{
	std::list<unsigned int> templist(getDimension(), 0);
	std::generate(
			templist.begin(),
			templist.end(),
			unique_number());
	_indexset.clear();
	_indexset.insert(templist.begin(), templist.end());
	assert( current_indexset.size() == getDimension());
}

unsigned int
LastNSearchDirections::updateIndexToMostParallel(
		const Norm &_Norm,
		const SpaceElement_ptr_t &_newdir) const
{
	// calculate the angles
	const angles_t angles = calculateBregmanAngles(
			_Norm,
			_newdir);

	angles_t::const_iterator indexiter = angles.begin();
	if (enforceRandomMapping) {
		// check whether we have to refresh current_indexset
		if (current_indexset.empty())
			replenishIndexset(current_indexset);

		// we look for largest element ourselves, constraint to
		// current_indexset
		double largest_angle = 0.;
		for (indexset_t::const_iterator iter = current_indexset.begin();
				iter != current_indexset.end(); ++iter) {
			if (angles[*iter] >= largest_angle) {
				largest_angle = angles[*iter];
				std::advance(
						indexiter,
						*iter - std::distance(angles.begin(), indexiter)
						);
			}
		}
	} else {
		unsigned int i=0;
		for (;i<getDimension();++i)
			if (U.col(i).isZero())
				break;
		// always fill a zero column if present (for the first N steps)
		indexiter =
				i < getDimension() ?
				angles.begin()+i :
				std::max_element(angles.begin(), angles.end());
	}

	// and return its index
	const unsigned int newindex =
			std::distance(
					const_cast<const angles_t &>(angles).begin(),
					indexiter);
	// also remove it from set
	if (enforceRandomMapping)
		current_indexset.erase(newindex);
	return newindex;
}

unsigned int
LastNSearchDirections::updateIndexToMostOrthogonal(
		const Norm &_Norm,
		const SpaceElement_ptr_t &_newdir) const
{
	// calculate the angles
	angles_t angles = calculateBregmanAngles(
			_Norm,
			_newdir);

	angles_t::const_iterator indexiter = angles.begin();
	if (enforceRandomMapping) {
		// check whether we have to refresh current_indexset
		if (current_indexset.empty())
			replenishIndexset(current_indexset);

		// we look for largest element ourselves, constraint to
		// current_indexset
		double smallest_angle = 1.;
		for (indexset_t::const_iterator iter = current_indexset.begin();
				iter != current_indexset.end(); ++iter) {
			if (angles[*iter] <= smallest_angle) {
				smallest_angle = angles[*iter];
				std::advance(
						indexiter,
						*iter - std::distance(
								const_cast<const angles_t &>(angles).begin(),
								indexiter)
						);
			}
		}
	} else {
		unsigned int i=0;
		for (;i<getDimension();++i)
			if (U.col(i).isZero())
				break;
		// find min angle (this automatically fills all zero columns)
		indexiter =
				i < getDimension() ?
				angles.begin()+i :
				std::min_element(angles.begin(), angles.end());
	}

	// and return its index
	const unsigned int newindex =
			std::distance(
						const_cast<const angles_t &>(angles).begin(),
						indexiter);
	// also remove it from set
	if (enforceRandomMapping)
		current_indexset.erase(newindex);
	return newindex;
}
