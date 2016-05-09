/*
 * LastNSearchDirections.hpp
 *
 *  Created on: Nov 17, 2014
 *      Author: heber
 */

#ifndef MINIMIZATIONS_MINIMIZERS_SEARCHSPACE_LASTNSEARCHDIRECTIONS_HPP_
#define MINIMIZATIONS_MINIMIZERS_SEARCHSPACE_LASTNSEARCHDIRECTIONS_HPP_

#include "BassoConfig.h"

#include "Searchspace.hpp"

#include <set>

#include "Minimizations/types.hpp"

class SequentialSubspaceMinimizer;

/** Uses the last N search directions.
 *
 */
class LastNSearchDirections : public Searchspace
{
public:
	//!> enumeration of all available orthogonalization types
	enum OrthogonalizationType {
		NoOrthogonalization=0,
		MetricOrthogonalization=1,
		BregmanOrthogonalization=2,
		MAX_Orthogonalization
	};

	/** Constructor of class LastNSearchDirections.
	 *
	 * Sets update function to simple advanceIndex().
	 *
	 * @param _SearchDirectionSpace_ptr search direction space (for checks)
	 * @param _N number of search directions in the subspace (i.e. search space)
	 * @param _orthogonalization_type orthogonalize new search direction?
	 */
	LastNSearchDirections(
			const NormedSpace_ptr_t &_SearchDirectionSpace_ptr,
			const unsigned int _N,
			const OrthogonalizationType _orthogonalization_type);

	/** This function performs the actual update of the search space.
	 *
	 * @param _newdir current search direction
	 * @param _alpha current alpha to this \a _newdir
	 * @param _dual_iterate not used
	 * @param _iterate not used
	 */
	void update(
			const SpaceElement_ptr_t &_newdir,
			const double _alpha,
			const SpaceElement_ptr_t &,
			const SpaceElement_ptr_t &
			);

	/** Const ref getter for \a index.
	 *
	 * @return current index
	 */
	const unsigned int getIndex() const {
		return index;
	}

	/** Const ref getter for \a lastIndices.
	 *
	 * @return array of offset w.r.t current search direction
	 */
	const std::vector<unsigned int> getLastIndices() const {
		return lastIndices;
	}

	//!> enumeration of available update index methods
	enum UpdateAlgorithmType {
		RoundRobin=0,
		MostParallel=1,
		MostOrthogonal=2,
		MAX_UpdateAlgorithmType
	};

	/** Returns a name for the given OrthogonalizationType.
	 *
	 * @param type chosen type
	 * @return string containing name of type
	 */
	static std::string getOrthogonalizationTypeName(
			const OrthogonalizationType type);

private:
	//!> grant SequentialSubspaceMinimizer access to updateIndex and enforceRandomMapping
	friend class SequentialSubspaceMinimizer;

	static UpdateAlgorithmType updateIndexType;

	//!> bool whether to enforce update index to be a random mapping or not
	static bool enforceRandomMapping;

	/** Setter for whether the update algorithm is forced to be a random
	 * mapping.
	 *
	 * @param _enforceRandomMapping true - enforce update algorithm to be
	 *		  random mapping, false - no enforcement
	 */
	void setEnforceRandomMapping(const bool _enforceRandomMapping)
	{  enforceRandomMapping = _enforceRandomMapping; }

	/** Function that simply advances the index in a round-robin fashion.
	 *
	 * @return index+1 mod N
	 */
	unsigned int advanceIndex() const;

	/** Function that compares new search direction with present ones
	 * in the state and selects the one that is most parallel,
	 * i.e. the one with maximum angle.
	 *
	 * @param _newdir new direction to compare to present ones
	 * @return index whose search direction is most parallel to \a newdir
	 */
	unsigned int
	updateIndexToMostParallel(
			const SpaceElement_ptr_t &_newdir) const;

	/** Function that compares new search direction with present ones
	 * in the state and selects the one that is most orthogonal,
	 * i.e. the one with minimum angle.
	 *
	 * @param _newdir new direction to compare to present ones
	 * @return index whose search direction is most orthogonal to \a newdir
	 */
	unsigned int
	updateIndexToMostOrthogonal(
			const SpaceElement_ptr_t &_newdir) const;

	//!> typedef for a set of indices
	typedef std::set<unsigned int> indexset_t;

	/** Helper function for enforcing RandomMapping.
	 *
	 * This recreates the full index set.
	 *
	 * @param _indexset indexset to replenish
	 */
	void
	replenishIndexset(
			indexset_t &_indexset) const;

protected:
	//!> index of the last updated search direction
	unsigned int index;

	/** Array of offsets to indicate sequence of search direction
	 * updates in SearchSpace, i.e. current search direction has offset
	 * of 0.
	 */
	std::vector<unsigned int> lastIndices;

	//!> current set of indices remaining for enforcing random mapping
	mutable indexset_t current_indexset;

	//!> whether and how to orthogonalize new search directions w.r.t old ones
	const OrthogonalizationType orthogonalization_type;

	//!> internal updated new direction for orthogonalization
	SpaceElement_ptr_t new_orthogonalized_dir;
};



#endif /* MINIMIZATIONS_MINIMIZERS_LASTNSEARCHDIRECTIONS_HPP_ */
