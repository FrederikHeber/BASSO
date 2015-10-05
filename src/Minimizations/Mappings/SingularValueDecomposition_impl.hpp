/*
 * SingularValueDecomposition_impl.hpp
 *
 *  Created on: Oct 5, 2015
 *      Author: heber
 */

#ifndef MINIMIZATIONS_MAPPINGS_SINGULARVALUEDECOMPOSITION_IMPL_HPP_
#define MINIMIZATIONS_MAPPINGS_SINGULARVALUEDECOMPOSITION_IMPL_HPP_

#include "BassoConfig.h"

#include <boost/shared_ptr.hpp>
#include <Eigen/Dense>

/** PIMPL pattern forSingularValueDecomposition.
 *
 */
struct SingularValueDecomposition_impl
{
	//!> wrap this entity in shared_ptr
	typedef boost::shared_ptr<SingularValueDecomposition_impl> ptr_t;

	SingularValueDecomposition_impl(const Eigen::MatrixXd &_matrix);

	//!> internal svd structure from Eigen
	const Eigen::JacobiSVD<Eigen::MatrixXd> svd;
};

#endif /* MINIMIZATIONS_MAPPINGS_SINGULARVALUEDECOMPOSITION_IMPL_HPP_ */

typedef SingularValueDecomposition_impl::ptr_t SingularValueDecomposition_impl_ptr_t;
