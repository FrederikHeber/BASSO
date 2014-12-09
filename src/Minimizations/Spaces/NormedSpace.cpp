/*
 * NormedSpace.cpp
 *
 *  Created on: Oct 23, 2014
 *      Author: heber
 */


#include "BassoConfig.h"

#include "NormedSpace.hpp"

#include <boost/bind.hpp>
#include <Eigen/Dense>

#include "Minimizations/Elements/SpaceElement.hpp"

NormedSpace::NormedSpace(
		const unsigned int _dimension) :
	dimension(_dimension),
	scalar_vector_fctor(
			boost::bind(
					static_cast<Eigen::internal::scalar_product_traits<typename Eigen::internal::traits<Eigen::VectorXd>::Scalar, typename Eigen::internal::traits<Eigen::VectorXd>::Scalar>::ReturnType
						(Eigen::MatrixBase<Eigen::VectorXd>::*)(const Eigen::MatrixBase<Eigen::VectorXd>&) const>(
								&Eigen::MatrixBase<Eigen::VectorXd>::dot),
								_1, _2
			)
	),
	ScalarVectorProduct(scalar_vector_fctor)
{}

NormedSpace::NormedSpace(
		const unsigned int _dimension,
		const Norm_ptr_t &_norm,
		const Mapping_ptr_t &_dualitymapping
		) :
	norm(_norm),
	dualitymapping(_dualitymapping),
	dimension(_dimension),
	scalar_vector_fctor(
			boost::bind(
					static_cast<Eigen::internal::scalar_product_traits<typename Eigen::internal::traits<Eigen::VectorXd>::Scalar, typename Eigen::internal::traits<Eigen::VectorXd>::Scalar>::ReturnType
						(Eigen::MatrixBase<Eigen::VectorXd>::*)(const Eigen::MatrixBase<Eigen::VectorXd>&) const>(
								&Eigen::MatrixBase<Eigen::VectorXd>::dot),
								_1, _2
			)
	),
	ScalarVectorProduct(scalar_vector_fctor)
{}

SpaceElement_ptr_t NormedSpace::createElement() const
{
	SpaceElement_ptr_t newelement(new SpaceElement(getSpace()));
	newelement->setSelfRef(newelement);
	return newelement;
}
