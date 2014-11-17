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
#include "Minimizations/Minimizers/Searchspace/LastNSearchDirections.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"

SequentialSubspaceMinimizer::IterationState::IterationState() :
	isInitialized(false)
{}

void SequentialSubspaceMinimizer::IterationState::set(
		const SpaceElement_ptr_t &_x0,
		const SpaceElement_ptr_t &_residual,
		const double _residuum,
		const unsigned int _N,
		const OperationCounter<
			Eigen::internal::scalar_product_traits<typename Eigen::internal::traits<Eigen::VectorXd>::Scalar, typename Eigen::internal::traits<Eigen::VectorXd>::Scalar>::ReturnType,
			const Eigen::MatrixBase<Eigen::VectorXd>&,
			const Eigen::MatrixBase<Eigen::VectorXd>&
			> &_ScalarVectorProduct_subspace
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
	// use last N search directions
	searchspace.reset(
			new LastNSearchDirections(
					_x0->getSpace()->getDualSpace(),
					_N,
					_ScalarVectorProduct_subspace));
	isInitialized = true;
}

void
SequentialSubspaceMinimizer::IterationState::reset()
{
	searchspace.reset();
	isInitialized = false;
}

unsigned int
SequentialSubspaceMinimizer::IterationState::getDimension() const
{
	if (isInitialized)
		return searchspace->getDimension();
	else
		return 0;
}

void
SequentialSubspaceMinimizer::IterationState::updateSearchSpace(
		const SpaceElement_ptr_t &_dualiterate,
		const SpaceElement_ptr_t &_newdir,
		const double _alpha)
{
	searchspace->update(_dualiterate, _newdir, _alpha);
}

const Eigen::MatrixXd &
SequentialSubspaceMinimizer::IterationState::getSearchSpace() const
{
	static const Eigen::MatrixXd zeroMatrix = Eigen::MatrixXd::Zero(0,0);
	if (isInitialized)
		return searchspace->getSearchSpace();
	else
		return zeroMatrix;
}

const Eigen::VectorXd &
SequentialSubspaceMinimizer::IterationState::getAlphas() const
{
	static const Eigen::VectorXd zeroVector = Eigen::VectorXd::Zero(0);
	if (isInitialized)
		return searchspace->getAlphas();
	else
		return zeroVector;
}

const SequentialSubspaceMinimizer::IterationState::angles_t
SequentialSubspaceMinimizer::IterationState::calculateAngles(
		const Norm &_Norm,
		const SpaceElement_ptr_t &_newdir) const
{
	return searchspace->calculateAngles(_Norm, _newdir);
}

const SequentialSubspaceMinimizer::IterationState::angles_t
SequentialSubspaceMinimizer::IterationState::calculateBregmanAngles(
		const Norm &_Norm,
		const SpaceElement_ptr_t &_newdir) const
{
	return searchspace->calculateBregmanAngles(
			_Norm, _newdir);
}
