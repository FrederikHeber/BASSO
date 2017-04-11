/*
 * Project: BASSO - BAnach Sequential Subspace Optimizer
 * Description: C++ library for solving optimization problems in Banach spaces
 * Copyright (C)  2014-2016 University of Saarland. All rights reserved.
 *
 *
 *   This file is part of the BASSO library.
 *
 *    BASSO is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 2 of the License, or
 *    (at your option) any later version.
 *
 *    BASSO is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with BASSO.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

/*
 * LastNSearchDirections.cpp
 *
 *  Created on: Nov 17, 2014
 *      Author: heber
 */


#include <Minimizations/Minimizers/Searchspace/LastNSearchDirections.hpp>
#include "BassoConfig.h"

#include <algorithm>
#include <functional>
#include <list>

#include "Log/Logging.hpp"

#include "Math/Helpers.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Functions/MetricProjectionFunctional.hpp"
#include "Minimizations/Functions/HyperplaneProjection.hpp"
#include "Minimizations/Functions/Minimizers/FunctionalMinimizerFactory.hpp"
#include "Minimizations/Functions/Minimizers/MinimizationFunctional.hpp"
#include "Minimizations/Functions/Minimizers/FunctionalMinimizerFactory.hpp"
#include "Minimizations/Functions/Minimizers/MinimizerExceptions.hpp"
#include "Minimizations/Mappings/DualityMapping.hpp"
#include "Minimizations/Mappings/PowerTypeDualityMapping.hpp"
#include "Minimizations/Norms/Norm.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"

// static entities
LastNSearchDirections::UpdateAlgorithmType
LastNSearchDirections::updateIndexType = LastNSearchDirections::RoundRobin;
bool LastNSearchDirections::enforceRandomMapping = false;

LastNSearchDirections::LastNSearchDirections(
		const NormedSpace_ptr_t &_SearchDirectionSpace_ptr,
		const unsigned int _N,
		const OrthogonalizationType _orthogonalization_type) :
	Searchspace(_SearchDirectionSpace_ptr,_N),
	index(0),
	lastIndices(_N, 0),
	orthogonalization_type(_orthogonalization_type),
	new_orthogonalized_dir(_SearchDirectionSpace_ptr->createElement())
{}

const unsigned int calculateMetricProjection(
		const SpaceElement_ptr_t& dual_solution,
		const std::vector<SpaceElement_ptr_t> &_searchspace,
		std::vector<double> & tmin
		)
{
	// some variables available in SequentialSubspaceMinimizer::calculateStepwidth()
	double TolY = 1e-6;
	const Norm &DualNormX = *dual_solution->getSpace()->getNorm();
	const Mapping & J_q = *dual_solution->getSpace()->getDualityMapping();
	const size_t MaxInnerIterations = 100;

	// consistency checks
	const size_t Ndirs = _searchspace.size();
	assert( Ndirs == tmin.size() );

	MetricProjectionFunctional metric(
			DualNormX, J_q, J_q.getPower(), _searchspace);
	const HyperplaneProjection<MetricProjectionFunctional> functional(
			metric, dual_solution);

	FunctionalMinimizer< std::vector<double> >::ptr_t functionminimizer =
			FunctionalMinimizerFactory::create< std::vector<double> >(
				Ndirs,
				functional,
				tmin);
	functionminimizer->setMaxIterations(MaxInnerIterations);

	unsigned int inner_iterations = (*functionminimizer)(
			Ndirs, TolY, tmin);

	return inner_iterations;
}

const unsigned int calculateBregmanProjection(
		const SpaceElement_ptr_t &_newdir,
		const SpaceElement_ptr_t &_olddir,
		std::vector<double> & tmin
		)
{
	const double searchdir_norm = _olddir->Norm();
	if (searchdir_norm < std::numeric_limits<double>::epsilon()) {
		tmin[0] = 0.;
		return 0;
	}
//			const std::pair<double, double> tmp =
//					projector(
//							_olddir,
//							_newdir,
//							1e-8);
//			const double projected_distance = tmp.second;
	const DualityMapping &J_q = static_cast<const DualityMapping &>(
			*_newdir->getSpace()->getDualityMapping());
	const double q = J_q.getPower();
	const double p = _newdir->getSpace()->getDualSpace()->getDualityMapping()->getPower();
	const double gamma_projected_distance =
			J_q(_newdir) * _olddir / searchdir_norm;
	const double searchdir_distance = searchdir_norm;
	tmin[0] = ::pow(
			gamma_projected_distance,
			p/q)/searchdir_distance;
//			const double projection_coefficient =
//					projected_distance/searchdir_distance;
//	LOG(info, "Projection coefficient is " << gamma_projected_distance << "/"
//		<< searchdir_distance << " = " << tmin[0]);
//	LOG(info, "Compare numerator to " << J_q(newdir) * _olddir / searchdir_norm);
	return 1;
}


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
			LOG(error, "Unknown updateIndex algorithm.");
			assert(0);
	}
	LOG(debug, "Updated index is " << index);
	lastIndices[index] = 0; // reset current search direction offset

	/// orthogonalize with respect to all old search directions
	if (orthogonalization_type != NoOrthogonalization) {
		// we need to first orthogonalize w.r.t last search direction, then
		// last but one, ...
		std::vector<unsigned int> orderOfApplication(lastIndices.size(), 0);
		for (size_t l = 0;l<lastIndices.size(); ++l)
			orderOfApplication[ lastIndices[l] ] = l;
		// put updated index (lastIndices[index]=0) at end of list as it is the
		// oldest search direction
		std::vector<unsigned int> indices_to_orthogonalize;
		indices_to_orthogonalize.insert(
				indices_to_orthogonalize.begin(),
				orderOfApplication.begin()+1,
				orderOfApplication.end());
		indices_to_orthogonalize.push_back(*orderOfApplication.begin());
		std::vector<double> tmin(indices_to_orthogonalize.size(), 0.);
		switch (orthogonalization_type) {
			case MetricOrthogonalization:
				calculateMetricProjection(
						_newdir,
						U,
						tmin
						);
				break;
			case BregmanOrthogonalization:
				for (std::vector<unsigned int>::const_iterator iter = indices_to_orthogonalize.begin();
						iter != indices_to_orthogonalize.end(); ++iter) {
					std::vector<double> tmin_tmp(1, 0.);
					calculateBregmanProjection(
							_newdir,
							U[ *iter ],
							tmin_tmp
							);
					tmin[ *iter ] = tmin_tmp[0];
				}
				break;
			default:
				LOG(error, "We cannot get here.");
				assert(0);
		}
		*new_orthogonalized_dir = _newdir;
		double new_orthogonalized_alpha = _alpha;
		for (std::vector<unsigned int>::const_iterator iter = indices_to_orthogonalize.begin();
				iter != indices_to_orthogonalize.end(); ++iter) {
			const double projection_coefficient = tmin[*iter];
			LOG(debug, "Projection coefficient of "
				<< (orthogonalization_type == MetricOrthogonalization ?
						"metric" : "bregman" ) << " orthogonalization is " << projection_coefficient);

			new_orthogonalized_dir->scaledAddition(-1.*projection_coefficient, U[ *iter ]);
			new_orthogonalized_alpha -= projection_coefficient * alphas[ *iter ];
		}
		*(U[index]) = new_orthogonalized_dir;
		alphas[index] = new_orthogonalized_alpha;
		LOG(debug, "Norm after projection changed from " << _newdir->Norm() << " to " << U[index]->Norm());
	} else {
		*(U[index]) = _newdir;
		alphas[index] = _alpha;
	}
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
			Helpers::unique_number());
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

std::string LastNSearchDirections::getOrthogonalizationTypeName(
		const OrthogonalizationType _type)
{
	switch (_type) {
	case NoOrthogonalization:
		return std::string("None");
		break;
	case MetricOrthogonalization:
		return std::string("metric");
		break;
	case BregmanOrthogonalization:
		return std::string("Bregman");
		break;
	default:
		return std::string("unknown");
		break;
	}
}
