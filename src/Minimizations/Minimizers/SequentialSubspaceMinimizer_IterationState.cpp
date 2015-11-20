/*
 * SequentialSubspaceMinimizer_IterationState.cpp
 *
 *  Created on: Sep 11, 2014
 *      Author: heber
 */


#include "BassoConfig.h"

#include "SequentialSubspaceMinimizer.hpp"

#include <algorithm>
#include <iterator>
#include <limits>
#include <list>
#include <vector>

#include <boost/log/trivial.hpp>

#include "Minimizations/types.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Minimizers/Searchspace/SearchspaceFactory.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"

SequentialSubspaceMinimizer::IterationState::IterationState() :
	isInitialized(false)
{
	status = notbegun;
}

void SequentialSubspaceMinimizer::IterationState::set(
		const SpaceElement_ptr_t &_x0,
		const SpaceElement_ptr_t &_dual_x0,
		const SpaceElement_ptr_t &_residual,
		const double _residuum,
		const unsigned int _N,
		const LastNSearchDirections::OrthogonalizationType _orthogonalization_type
		)
{
	/// -# initialize return structure
	NumberOuterIterations = 0;
	// set iterate 'x' as start vector 'x0'
	m_solution = _x0->getSpace()->createElement();
	m_solution = _x0;
	m_dual_solution = _dual_x0->getSpace()->createElement();
	m_dual_solution = _dual_x0;
	// calculate starting residual and norm
	residuum = _residuum;
	m_residual = _residual->getSpace()->createElement();
	*m_residual = _residual;
	// use pre-specified type
	searchspace = SearchspaceFactory::create(
					_x0->getSpace()->getDualSpace(),
					_N,
					_orthogonalization_type);
	isInitialized = true;
	status = starting;
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
		const SpaceElement_ptr_t &_newdir,
		const double _alpha,
		const SpaceElement_ptr_t &_dual_iterate,
		const SpaceElement_ptr_t &_iterate)
{
	searchspace->update(_newdir, _alpha, _dual_iterate, _iterate);
}

const std::vector<SpaceElement_ptr_t> &
SequentialSubspaceMinimizer::IterationState::getSearchSpace() const
{
	static const std::vector<SpaceElement_ptr_t> zeroSearchspace;
	if (isInitialized)
		return searchspace->getSearchSpace();
	else
		return zeroSearchspace;
}

const std::vector<double> &
SequentialSubspaceMinimizer::IterationState::getAlphas() const
{
	static const std::vector<double> zeroVector;
	if (isInitialized)
		return searchspace->getAlphas();
	else
		return zeroVector;
}

const SequentialSubspaceMinimizer::IterationState::angles_t
SequentialSubspaceMinimizer::IterationState::calculateAngles(
		const SpaceElement_ptr_t &_newdir) const
{
	return searchspace->calculateAngles(_newdir);
}

const SequentialSubspaceMinimizer::IterationState::angles_t
SequentialSubspaceMinimizer::IterationState::calculateBregmanAngles(
		const SpaceElement_ptr_t &_newdir) const
{
	return searchspace->calculateBregmanAngles(_newdir);
}
