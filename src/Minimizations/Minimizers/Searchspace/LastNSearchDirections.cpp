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

#include "Math/Helpers.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Functions/MetricProjectionFunctional.hpp"
#include "Minimizations/Functions/HyperplaneProjection.hpp"
#include "Minimizations/Functions/Minimizers/FunctionalMinimizer.hpp"
#include "Minimizations/Functions/Minimizers/MinimizationFunctional.hpp"
#include "Minimizations/Functions/Minimizers/Minimizer.hpp"
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
		const bool _orthogonal_directions) :
	Searchspace(_SearchDirectionSpace_ptr,_N),
	index(0),
	lastIndices(_N, 0),
	orthogonal_directions(_orthogonal_directions)
{}

const unsigned int calculateStepWidth(
		const SpaceElement_ptr_t& dual_solution,
		std::vector<double> & tmin,
		const std::vector<SpaceElement_ptr_t> &_searchspace,
		const std::vector<double> &_alphas
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

	// tmin=fminunc(@(t) MetricProjectionFunctional(t,Jx,u,alpha+d,DualNormX,DualPowerX,TolX),t0,BregmanOptions);
	MetricProjectionFunctional metric(DualNormX,
			dynamic_cast<const PowerTypeDualityMapping&>(J_q),
			J_q.getPower(), _searchspace);
	const HyperplaneProjection<MetricProjectionFunctional> functional(
			metric, dual_solution);

	// due to templation we need to instantiate both, as user
	// decides during runtime which we need
	Minimizer<gsl_vector> minimizer_gsl(Ndirs);
//	Minimizer<NLopt_vector> minimizer_nlopt(Ndirs);
//	if (MinLib == gnuscientificlibrary) {
		minimizer_gsl.setMaxIterations(MaxInnerIterations);
//	} else if (MinLib == nonlinearoptimization) {
//		minimizer_nlopt.setMaxIterations(MaxInnerIterations);
//	}

	unsigned int inner_iterations = 0;
	FunctionalMinimizer<std::vector<double>, gsl_vector> fmin_gsl(
			functional, minimizer_gsl, tmin);
//	FunctionalMinimizer<std::vector<double>, NLopt_vector> fmin_nlopt(
//			functional, minimizer_nlopt, tmin);
//	switch (MinLib) {
//	case gnuscientificlibrary:
		inner_iterations = fmin_gsl(Ndirs, TolY, tmin);
//		break;
//	case nonlinearoptimization:
//		inner_iterations = fmin_nlopt(Ndirs, TolY, tmin);
//		break;
//	default:
//		throw MinimizationIllegalValue_exception()
//				<< MinimizationIllegalValue_name("MinLib");
//		break;
//	}
	std::stringstream output_stepwidth;
	std::copy(tmin.begin(), tmin.end(),
			std::ostream_iterator<double>(output_stepwidth, " "));
	BOOST_LOG_TRIVIAL(debug)<< "tmin " << output_stepwidth.str()
	<< " found in " << inner_iterations << " iterations.";

	return inner_iterations;
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
		double alpha = _alpha;
		// we need to first orthogonalize w.r.t last search direction, then
		// last but one, ...
		std::vector<unsigned int> orderOfApplication(lastIndices.size(), 0);
		for (size_t l = 0;l<lastIndices.size(); ++l)
			orderOfApplication[ lastIndices[l] ] = l;
		// put updated index (lastIndices[index]=0) at end of list as it is the
		// oldes search direction
		std::vector<unsigned int> indices_to_orthogonalize;
		indices_to_orthogonalize.insert(
				indices_to_orthogonalize.begin(),
				orderOfApplication.begin()+1,
				orderOfApplication.end());
		indices_to_orthogonalize.push_back(*orderOfApplication.begin());
		for (std::vector<unsigned int>::const_iterator iter = indices_to_orthogonalize.begin();
				iter != indices_to_orthogonalize.end(); ++iter) {
			const std::vector<SpaceElement_ptr_t> searchspace(1, U[ *iter ]);
			const std::vector<double> alphas(1, 0.);
			std::vector<double> tmin(1, 0.);
			calculateStepWidth(
					newdir,
					tmin,
					searchspace,
					alphas
					);
			const double projection_coefficient = tmin[0];
			BOOST_LOG_TRIVIAL(info)
				<< "Projection coefficient is " << projection_coefficient;

//			const double searchdir_norm =  U[ *iter ]->Norm();
//			if (searchdir_norm < std::numeric_limits<double>::epsilon())
//				continue;
////			const std::pair<double, double> tmp =
////					projector(
////							U[ *iter ],
////							newdir,
////							1e-8);
////			const double projected_distance = tmp.second;
//			const DualityMapping &J_q = static_cast<const DualityMapping &>(
//					*newdir->getSpace()->getDualityMapping());
//			const double q = J_q.getPower();
//			const double p = newdir->getSpace()->getDualSpace()->getDualityMapping()->getPower();
//			const double gamma_projected_distance =
//					J_q(newdir) * U[ *iter ] / searchdir_norm;
//			const double searchdir_distance = searchdir_norm;
//			const double projection_coefficient = ::pow(
//					gamma_projected_distance,
//					p/q)/searchdir_distance;
////			const double projection_coefficient =
////					projected_distance/searchdir_distance;
//			BOOST_LOG_TRIVIAL(info)
//				<< "Projection coefficient is " << gamma_projected_distance << "/"
//				<< searchdir_distance << " = " << projection_coefficient;
////			BOOST_LOG_TRIVIAL(info)
////				<< "Compare numerator to "
////				<< J_q(newdir) * U[ *iter ] / searchdir_norm;
				*newdir -= projection_coefficient * U[ *iter ];
				alpha -= projection_coefficient * alphas[ *iter ];
		}
		const double prenorm = _newdir->Norm();
		const double postnorm = newdir->Norm();
		BOOST_LOG_TRIVIAL(info)
			<< "Norm after Bregman projection changed from "
			<< prenorm << " to " << postnorm;

		*(U[index]) = newdir;
		alphas[index] = alpha;
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
