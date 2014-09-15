/*
 * SequentialSubspaceMinimizer_IterationState.cpp
 *
 *  Created on: Sep 11, 2014
 *      Author: heber
 */


#include "BassoConfig.h"

#include "SequentialSubspaceMinimizer.hpp"

#include <algorithm>
#include <Eigen/Dense>
#include <iterator>
#include <limits>
#include <list>
#include <vector>

#include <boost/log/trivial.hpp>

#include "Minimizations/types.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
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

SequentialSubspaceMinimizer::IterationState::IterationState() :
	isInitialized(false),
	index(-1),
	updateIndex(&SequentialSubspaceMinimizer::IterationState::advanceIndex)
{}

void SequentialSubspaceMinimizer::IterationState::set(
		const SpaceElement_ptr_t &_x0,
		const SpaceElement_ptr_t &_residual,
		const double _residuum,
		const unsigned int _N
		)
{
	/// -# initialize return structure
	NumberOuterIterations = 0;
	// set iterate 'x' as start vector 'x0'
	m_solution = _x0->getSpace()->createElement();
	m_solution = _x0;
	// calculate starting residual and norm
	residuum = _residuum;
	m_residual = _residual->getSpace()->createElement();
	*m_residual = _residual;
	U = Eigen::MatrixXd::Zero(_x0->getSpace()->getDimension(),_N);
	alphas = Eigen::VectorXd::Zero(_N);
	index = 0;
	isInitialized = true;
}

const Eigen::MatrixXd &
SequentialSubspaceMinimizer::IterationState::getSearchSpace() const
{
	static const Eigen::MatrixXd zeroMatrix = Eigen::MatrixXd::Zero(0,0);
	if (isInitialized)
		return U;
	else
		return zeroMatrix;
}

const Eigen::VectorXd &
SequentialSubspaceMinimizer::IterationState::getAlphas() const
{
	static const Eigen::VectorXd zeroVector = Eigen::VectorXd::Zero(0);
	if (isInitialized)
		return alphas;
	else
		return zeroVector;
}

void
SequentialSubspaceMinimizer::IterationState::updateSearchSpace(
		const SpaceElement_ptr_t &_newdir,
		const double _alpha)
{
	index = updateIndex(this, _newdir);
	BOOST_LOG_TRIVIAL(debug)
		<< "Updated index is " << index;
	U.col(index) = _newdir->getVectorRepresentation();
	alphas(index) = _alpha;
}

unsigned int
SequentialSubspaceMinimizer::IterationState::advanceIndex(
		const SpaceElement_ptr_t &_newdir) const
{
	return ((index + 1) % getDimension());
}

const SequentialSubspaceMinimizer::IterationState::angles_t
SequentialSubspaceMinimizer::IterationState::calculateBregmanAngles(
		const Norm &_Norm,
		const VectorProjection &_projector,
		const SpaceElement_ptr_t &_newdir) const
{
	// calculate all angles
	const Eigen::MatrixXd &U = getSearchSpace();
	const unsigned int N = getDimension();
	typedef std::vector<double> angles_t;
	angles_t angles(N, 0.);

	for (unsigned int l=0;l<N;++l) {
		if (U.col(l).norm() < std::numeric_limits<double>::epsilon())
			continue;
		// first: minimum, second: minimizer (equals length in LpNorm here)
		const std::pair<double, double> tmp =
				_projector(
						U.col(l),
						_newdir->getVectorRepresentation(),
						1e-4);
		const double projected_distance = tmp.second;
		const double original_distance = _Norm(_newdir);
		if (fabs(original_distance) > std::numeric_limits<double>::epsilon()*1e2) {
			angles[l] = fabs(projected_distance / original_distance);
		} else {
			angles[l] = 0.;
		}

		BOOST_LOG_TRIVIAL(info)
			<< "Bregman Angles #" << l << " is " << angles[l];
	}

	return angles;
}

const SequentialSubspaceMinimizer::IterationState::angles_t
SequentialSubspaceMinimizer::IterationState::calculateAngles(
		const Norm &_Norm,
		const SpaceElement_ptr_t &_newdir) const
{
	// calculate all angles
	const Eigen::MatrixXd &U = getSearchSpace();
	const unsigned int N = getDimension();
	angles_t angles(N, 0.);

	for (unsigned int l=0;l<N;++l) {
		if (U.col(l).norm() < std::numeric_limits<double>::epsilon()) {
			continue;
		}
		// first: minimum, second: minimizer (equals length in LpNorm here)
		const Eigen::VectorXd Ucol_transposed = U.col(l).transpose();
		const double projected_distance = Ucol_transposed.dot(
				_newdir->getVectorRepresentation())
				/ _Norm(U.col(l));
		const double original_distance = _Norm(_newdir);
		if (fabs(original_distance) > std::numeric_limits<double>::epsilon()*1e2) {
			angles[l] = fabs(projected_distance / original_distance);
		} else {
			angles[l] = 0.;
		}
		BOOST_LOG_TRIVIAL(info)
			<< "Angles #" << l << " is " << angles[l];
	}

	return angles;
}

void
SequentialSubspaceMinimizer::IterationState::replenishIndexset(
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
SequentialSubspaceMinimizer::IterationState::updateIndexToMostParallel(
		const Norm &_Norm,
		const VectorProjection &_projector,
		const SpaceElement_ptr_t &_newdir) const
{
	// calculate the angles
	const angles_t angles = calculateBregmanAngles(
			_Norm,
			_projector,
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
		// always fill a zero column if present (for the first N steps)
		indexiter =
				(unsigned int)NumberOuterIterations < getDimension() ?
				angles.begin()+(unsigned int)NumberOuterIterations :
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
SequentialSubspaceMinimizer::IterationState::updateIndexToMostOrthogonal(
		const Norm &_Norm,
		const VectorProjection &_projector,
		const SpaceElement_ptr_t &_newdir) const
{
	// calculate the angles
	angles_t angles = calculateBregmanAngles(
			_Norm,
			_projector,
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
		// find min angle (this automatically fills all zero columns)
		indexiter =
				(unsigned int)NumberOuterIterations < getDimension() ?
				angles.begin()+(unsigned int)NumberOuterIterations :
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
