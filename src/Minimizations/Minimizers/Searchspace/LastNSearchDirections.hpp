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
	/** Constructor of class LastNSearchDirections.
	 *
	 * Sets update function to simple advanceIndex().
	 *
	 * @param _SearchDirectionSpace_ptr search direction space (for checks)
	 * @param _N number of search directions
	 * @param _ScalarVectorProduct_subspace counts for vector-vector
	 * products in the subspace (i.e. search space)
	 */
	LastNSearchDirections(
			const NormedSpace_ptr_t &_SearchDirectionSpace_ptr,
			const unsigned int _N,
			const OperationCounter<
				Eigen::internal::scalar_product_traits<typename Eigen::internal::traits<Eigen::VectorXd>::Scalar, typename Eigen::internal::traits<Eigen::VectorXd>::Scalar>::ReturnType,
				const Eigen::MatrixBase<Eigen::VectorXd>&,
				const Eigen::MatrixBase<Eigen::VectorXd>&
				> &_ScalarVectorProduct_subspace);

	/** This function performs the actual update of the search space.
	 *
	 * @param _iterate current dual iterate
	 * @param _newdir current search direction
	 * @param _alpha current alpha to this \a _newdir
	 */
	void update(
			const SpaceElement_ptr_t &_iterate,
			const SpaceElement_ptr_t &_newdir,
			const double _alpha);

	/** Const ref getter for \a index.
	 *
	 * @return current index
	 */
	const unsigned int getIndex() const {
		return index;
	}

	//!> enumeration of available update index methods
	enum UpdateAlgorithmType {
		RoundRobin=0,
		MostParallel=1,
		MostOrthogonal=2,
		MAX_UpdateAlgorithmType
	};

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
	 * @param _Norm norm object to calculate norms
	 * @param _newdir new direction to compare to present ones
	 * @return index whose search direction is most parallel to \a newdir
	 */
	unsigned int
	updateIndexToMostParallel(
			const Norm &_Norm,
			const SpaceElement_ptr_t &_newdir) const;

	/** Function that compares new search direction with present ones
	 * in the state and selects the one that is most orthogonal,
	 * i.e. the one with minimum angle.
	 *
	 * @param _Norm norm object to calculate norms
	 * @param _newdir new direction to compare to present ones
	 * @return index whose search direction is most orthogonal to \a newdir
	 */
	unsigned int
	updateIndexToMostOrthogonal(
			const Norm &_Norm,
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

	//!> current set of indices remaining for enforcing random mapping
	mutable indexset_t current_indexset;
};



#endif /* MINIMIZATIONS_MINIMIZERS_LASTNSEARCHDIRECTIONS_HPP_ */
