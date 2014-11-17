/*
 * Searchspace.cpp
 *
 *  Created on: Nov 17, 2014
 *      Author: heber
 */


#include "BassoConfig.h"

#include "Searchspace.hpp"

#include <cmath>
#include <limits>

#include "Log/Logging.hpp"

#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Functions/VectorProjection.hpp"
#include "Minimizations/Mappings/PowerTypeDualityMapping.hpp"
#include "Minimizations/Norms/Norm.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"

Searchspace::Searchspace(
		const NormedSpace_ptr_t &_SearchDirectionSpace_ptr,
		const unsigned int _N,
		const OperationCounter<
				Eigen::internal::scalar_product_traits<typename Eigen::internal::traits<Eigen::VectorXd>::Scalar, typename Eigen::internal::traits<Eigen::VectorXd>::Scalar>::ReturnType,
				const Eigen::MatrixBase<Eigen::VectorXd>&,
				const Eigen::MatrixBase<Eigen::VectorXd>&
				> &_ScalarVectorProduct_subspace) :
	SearchDirectionSpace_ptr(_SearchDirectionSpace_ptr),
	projector(*_SearchDirectionSpace_ptr->getNorm(),
			dynamic_cast<const PowerTypeDualityMapping &>(
					*_SearchDirectionSpace_ptr->getDualityMapping()
					),
					_SearchDirectionSpace_ptr->getDualityMapping()->getPower(),
			_ScalarVectorProduct_subspace)
{
	U = Eigen::MatrixXd::Zero(
			SearchDirectionSpace_ptr->getDimension(),
			_N);
	alphas = Eigen::VectorXd::Zero(_N);
}

unsigned int
Searchspace::getDimension() const
{
	return getSearchSpace().outerSize();
}

const Searchspace::angles_t
Searchspace::calculateBregmanAngles(
		const Norm &_Norm,
		const SpaceElement_ptr_t &_newdir) const
{
	// calculate all angles
	const Eigen::MatrixXd &U = getSearchSpace();
	const unsigned int N = getDimension();
	angles_t angles(N, 0.);

	for (unsigned int l=0;l<N;++l) {
		if (U.col(l).norm() < std::numeric_limits<double>::epsilon())
			continue;
		// first: minimum, second: minimizer (equals length in LpNorm here)
		const std::pair<double, double> tmp =
				projector(
						U.col(l),
						_newdir->getVectorRepresentation(),
						1e-4);
		const double projected_distance = tmp.second;
		const double original_distance = _Norm(_newdir);
		if (fabs(original_distance) > std::numeric_limits<double>::epsilon()*1e2) {
			angles[l] = fabs(projected_distance / original_distance);
		} else {
			angles[l] = 0.;
		}

		BOOST_LOG_TRIVIAL(debug)
			<< "Bregman Angles #" << l << " is " << angles[l];
	}

	return angles;
}

const Searchspace::angles_t
Searchspace::calculateAngles(
		const Norm &_Norm,
		const SpaceElement_ptr_t &_newdir) const
{
	// calculate all angles
	const Eigen::MatrixXd &U = getSearchSpace();
	const unsigned int N = getDimension();
	angles_t angles(N, 0.);

	for (unsigned int l=0;l<N;++l) {
		if (U.col(l).norm() < std::numeric_limits<double>::epsilon()) {
			continue;
		}
		// first: minimum, second: minimizer (equals length in LpNorm here)
		const Eigen::VectorXd Ucol_transposed = U.col(l).transpose();
		const double projected_distance = Ucol_transposed.dot(
				_newdir->getVectorRepresentation())
				/ _Norm(U.col(l));
		const double original_distance = _Norm(_newdir);
		if (fabs(original_distance) > std::numeric_limits<double>::epsilon()*1e2) {
			angles[l] = fabs(projected_distance / original_distance);
		} else {
			angles[l] = 0.;
		}
		BOOST_LOG_TRIVIAL(debug)
			<< "Angles #" << l << " is " << angles[l];
	}

	return angles;
}
