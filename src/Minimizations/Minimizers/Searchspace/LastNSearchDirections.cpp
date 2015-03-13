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
#include <functional>
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
		const unsigned int _N,
		const bool _orthogonal_directions) :
	Searchspace(_SearchDirectionSpace_ptr,_N),
	index(0),
	lastIndices(_N, 0),
	orthogonal_directions(_orthogonal_directions)
{}

void LastNSearchDirections::update(
		const SpaceElement_ptr_t &_newdir,
		const double _alpha,
		const SpaceElement_ptr_t &,
		const SpaceElement_ptr_t &
		)
{
	assert( _newdir->getSpace() == SearchDirectionSpace_ptr );

	// add one to all in lastIndices
	std::transform(
			lastIndices.begin(), lastIndices.end(),
			lastIndices.begin(), std::bind2nd(std::plus<unsigned int>(), 1));

	switch (updateIndexType) {
		case RoundRobin:
			index = advanceIndex();
			break;
		case MostParallel:
			index = updateIndexToMostParallel(_newdir);
			break;
		case MostOrthogonal:
			index = updateIndexToMostOrthogonal(_newdir);
			break;
		default:
			BOOST_LOG_TRIVIAL(error)
				<< "Unknown updateIndex algorithm.";
			assert(0);
	}
	BOOST_LOG_TRIVIAL(debug)
		<< "Updated index is " << index;
	lastIndices[index] = 0; // reset current search direction offset

	/// orthogonalize with respect to all old search directions
	if (orthogonal_directions) {
		SpaceElement_ptr_t newdir = _newdir->getSpace()->createElement();
		*newdir = _newdir;
		// we need to first orthogonalize w.r.t last search direction, then
		// last but one, ...
		std::vector<unsigned int> orderOfApplication(lastIndices.size(), 0);
		for (size_t l = 0;l<lastIndices.size(); ++l)
			orderOfApplication[ lastIndices[l] ] = l;
		const double prenorm = newdir->Norm();
		for (size_t l = 0;l<orderOfApplication.size(); ++l) {
			const double searchdir_norm = U[ orderOfApplication[l] ]->Norm();
			if (searchdir_norm < std::numeric_limits<double>::epsilon())
				continue;
			const std::pair<double, double> tmp =
					projector(
							U[ orderOfApplication[l] ],
							newdir,
							1e-8);
			const double projected_distance = tmp.second;
			const double searchdir_distance = pow(searchdir_norm,2);
			*newdir -=
					projected_distance/searchdir_distance * U[ orderOfApplication[l] ];
		}
		const double postnorm = newdir->Norm();
		BOOST_LOG_TRIVIAL(info)
			<< "Norm after Bregman projection changed from "
			<< prenorm << " to " << postnorm;

		*(U[index]) = newdir;
	} else {
		*(U[index]) = _newdir;
	}

	// TODO: Can we use the "old" alpha or do we require a projection
	// here as well?
	alphas[index] = _alpha;
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
		const SpaceElement_ptr_t &_newdir) const
{
	// calculate the angles
	const angles_t angles = calculateBregmanAngles(_newdir);

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
			if (U[i]->isZero())
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
		const SpaceElement_ptr_t &_newdir) const
{
	// calculate the angles
	angles_t angles = calculateBregmanAngles(_newdir);

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
			if (U[i]->isZero())
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
