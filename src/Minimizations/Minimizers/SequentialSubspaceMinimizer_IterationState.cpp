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

#include <boost/log/trivial.hpp>

SequentialSubspaceMinimizer::IterationState::IterationState() :
	isInitialized(false),
	index(-1),
	updateIndex(&SequentialSubspaceMinimizer::IterationState::advanceIndex)
{}

void SequentialSubspaceMinimizer::IterationState::set(
		const unsigned int _dimension,
		const unsigned int _N
		)
{
	U = Eigen::MatrixXd::Zero(_dimension,_N);
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
		const Eigen::VectorXd &_newdir,
		const Eigen::VectorXd &_iterate,
		const double _alpha)
{
	index = updateIndex(this, _newdir, _iterate);
	U.col(index) = _newdir;
	alphas(index) = _alpha;
}

unsigned int
SequentialSubspaceMinimizer::IterationState::advanceIndex(
		const Eigen::VectorXd &_newdir,
		const Eigen::VectorXd &_iterate) const
{
	return ((index + 1) % getDimension());
}
