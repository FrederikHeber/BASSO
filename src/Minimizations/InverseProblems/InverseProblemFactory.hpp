/*
 * InverseProblemFactory.hpp
 *
 *  Created on: Oct 27, 2014
 *      Author: heber
 */

#ifndef INVERSEPROBLEMFACTORY_HPP_
#define INVERSEPROBLEMFACTORY_HPP_

#include "BassoConfig.h"

#include <Eigen/Dense>

#include "Minimizations/types.hpp"

/** This factory creates InverseProblem instances, that are the abstract
 * interface with everything the Minimizers need, from given parameters and
 * structures, i.e. eventually we always want to solve \f$ Ax=y \f$ where
 * \a A is a matrix, \a y is the right-hand side, and we are looking for
 * the right \a x.
 */
struct InverseProblemFactory
{
	/** Creates an inverse problem in Lp spaces.
	 *
	 * @param _p p value of the space X
	 * @param _powerx power for the duality mapping X -> dualX
	 * @param _r p value of the space Y
	 * @param _powery power for the duality mapping Y -> dualY
	 * @param _matrix matrix as the finite-dimensional representation of
	 * 		  the linear mapping X -> Y
	 * @param _rhs finite-dimensional representation of right-hand side.
	 * @return
	 */
	static InverseProblem_ptr_t createLpInstance(
			const double _p,
			const double _powerx,
			const double _r,
			const double _powery,
			const Eigen::MatrixXd &_matrix,
			const Eigen::VectorXd &_rhs);

	/** Creates an inverse problem in regalurized L1 space.
	 *
	 * @param _lambda regularization parameter
	 * @param _r p value of the space Y
	 * @param _powery power for the duality mapping Y -> dualY
	 * @param _matrix matrix as the finite-dimensional representation of
	 * 		  the linear mapping X -> Y
	 * @param _rhs finite-dimensional representation of right-hand side.
	 * @return
	 */
	static InverseProblem_ptr_t createRegularizedL1Instance(
			const double _lambda,
			const double _r,
			const double _powery,
			const Eigen::MatrixXd &_matrix,
			const Eigen::VectorXd &_rhs);
};



#endif /* INVERSEPROBLEMFACTORY_HPP_ */
