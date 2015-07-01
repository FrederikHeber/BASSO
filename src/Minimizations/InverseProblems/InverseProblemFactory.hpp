/*
 * InverseProblemFactory.hpp
 *
 *  Created on: Oct 27, 2014
 *      Author: heber
 */

#ifndef INVERSEPROBLEMFACTORY_HPP_
#define INVERSEPROBLEMFACTORY_HPP_

#include "BassoConfig.h"

#include <string>
#include <vector>

#include <boost/any.hpp>

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
	//!> typedef for the vector of arbitrary arguments.
	typedef std::vector<boost::any> args_t;

	/** Creates an inverse problem in Lp spaces.
	 *
	 * @param  _type_SpaceX type of space for X
	 * @param _args_SpaceX arguments required to create X
	 * @param _type_SpaceY type of space for Y
	 * @param _args_SpaceY arguments required to create Y
	 * @param _matrix matrix as the finite-dimensional representation of
	 * 		  the linear mapping X -> Y
	 * @param _rhs finite-dimensional representation of right-hand side.
	 * @return created instance
	 */
	static InverseProblem_ptr_t create(
			const std::string &_type_SpaceX,
			const args_t &_args_SpaceX,
			const std::string &_type_SpaceY,
			const args_t &_args_SpaceY,
			const Eigen::MatrixXd &_matrix,
			const Eigen::VectorXd &_rhs);
};



#endif /* INVERSEPROBLEMFACTORY_HPP_ */
