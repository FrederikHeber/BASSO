/*
 * Searchspace.hpp
 *
 *  Created on: Nov 17, 2014
 *      Author: heber
 */

#ifndef MINIMIZATIONS_MINIMIZERS_SEARCHSPACE_SEARCHSPACE_HPP_
#define MINIMIZATIONS_MINIMIZERS_SEARCHSPACE_SEARCHSPACE_HPP_

#include "BassoConfig.h"

#include <Eigen/Dense>
#include <vector>

#include "MatrixIO/OperationCounter.hpp"

#include "Minimizations/types.hpp"
#include "Minimizations/Functions/VectorProjection.hpp"

class Norm;

/** This class contains the search space and interface to access
 * its representations.
 */
class Searchspace
{
public:
	/** Constructor for class Searchspace, initializes matrix and
	 * vector representations.
	 *
	 * @param _SearchDirectionSpace_ptr search direction space (for checks)
	 * @param _N number of search directions
	 * @param _ScalarVectorProduct_subspace counts for vector-vector
	 * products in the subspace (i.e. search space)
	 */
	Searchspace(
			const NormedSpace_ptr_t &_SearchDirectionSpace_ptr,
			const unsigned int _N,
			const OperationCounter<
				Eigen::internal::scalar_product_traits<typename Eigen::internal::traits<Eigen::VectorXd>::Scalar, typename Eigen::internal::traits<Eigen::VectorXd>::Scalar>::ReturnType,
				const Eigen::MatrixBase<Eigen::VectorXd>&,
				const Eigen::MatrixBase<Eigen::VectorXd>&
				> &_ScalarVectorProduct_subspace);

	virtual ~Searchspace() {}

	/** Const ref getter for \a U.
	 *
	 * @return const ref to U
	 */
	const Eigen::MatrixXd & getSearchSpace() const
	{ return U; }

	/** Const ref getter for \a alphas.
	 *
	 * @return const ref to alphas
	 */
	const Eigen::VectorXd & getAlphas() const
	{ return alphas; }

	/** Getter for the dimension of the search directions in \a U.
	 *
	 * @return size of search direction vectors,
	 * 		   0 - if state not initialized
	 */
	virtual unsigned int getDimension() const;

	/** Returns the index of the current search direction in the search space.
	 *
	 * @return index of the current search direction
	 */
	virtual const unsigned int getIndex() const = 0;

	/** This function performs the update of the search space.
	 *
	 * @param _iterate current dual iterate
	 * @param _newdir current search direction
	 * @param _alpha current alpha to this \a _newdir
	 */
	virtual void update(
			const SpaceElement_ptr_t &_iterate,
			const SpaceElement_ptr_t &_newdir,
			const double _alpha) = 0;

	//!> typedef for a vector of angles
	typedef std::vector<double> angles_t;

	/** Helper function to calculate the angles between each search
	 * direction in ::U and the given _newdir.
	 *
	 * @param _Norm norm object to calculate norms
	 * @param _newdir new direction to compare to present ones
	 * @return vector of doubles, the angles
	 */
	const angles_t
	calculateAngles(
			const Norm &_Norm,
			const SpaceElement_ptr_t &_newdir) const;

	/** Helper function to calculate the angles between each search
	 * direction in ::U and the given _newdir using Bregman projections
	 * and distance.
	 *
	 * @param _Norm norm object to calculate norms
	 * @param _newdir new direction to compare to present ones
	 * @return vector of doubles, the angles
	 */
	const angles_t
	calculateBregmanAngles(
			const Norm &_Norm,
			const SpaceElement_ptr_t &_newdir) const;

protected:
	//!> reference to Space of search directions for checking
	NormedSpace_ptr_t SearchDirectionSpace_ptr;

	//!> Vector Projection instance for calculating angles in Banach space
	VectorProjection projector;

	//!> subspace matrix with search directions as column vectors
	Eigen::MatrixXd U;
	//!> offset of hyperplanes of search directions for projection
	Eigen::VectorXd alphas;
};



#endif /* MINIMIZATIONS_MINIMIZERS_SEARCHSPACE_SEARCHSPACE_HPP_ */
