/*
 * SingularValueDecomposition_impl.cpp
 *
 *  Created on: Oct 5, 2015
 *      Author: heber
 */


#include "BassoConfig.h"

#include "SingularValueDecomposition_impl.hpp"

SingularValueDecomposition_impl::SingularValueDecomposition_impl(
		const Eigen::MatrixXd &_matrix) :
	svd(_matrix.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV))
{}

