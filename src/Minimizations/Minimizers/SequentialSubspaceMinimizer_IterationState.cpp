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
#include <limits>
#include <vector>

#include <boost/log/trivial.hpp>

#include "Minimizations/types.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"

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
		if (fabs(projected_distance) > std::numeric_limits<double>::epsilon()*1e2) {
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
		if (fabs(projected_distance) > std::numeric_limits<double>::epsilon()*1e2) {
			angles[l] = fabs(projected_distance / original_distance);
		} else {
			angles[l] = 0.;
		}
		BOOST_LOG_TRIVIAL(info)
			<< "Angles #" << l << " is " << angles[l];
	}

	return angles;
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

	// always fill a zero column if present (for the first N steps)
	const angles_t::const_iterator indexiter =
			std::max_element(angles.begin(), angles.end());

	// and return its index
	return std::distance(
			const_cast<const angles_t &>(angles).begin(),
			indexiter);
}
