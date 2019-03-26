/*
 * Project: BASSO - BAnach Sequential Subspace Optimizer
 * Description: C++ library for solving optimization problems in Banach spaces
 * Copyright (C)  2014-2016 University of Saarland. All rights reserved.
 * Copyright (C)  2016-2019 Frederik Heber. All rights reserved.
 *
 *
 *   This file is part of the BASSO library.
 *
 *    BASSO is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 2 of the License, or
 *    (at your option) any later version.
 *
 *    BASSO is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with BASSO.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

/*
 * TruncatedSVD.cpp
 *
 *  Created on: Mar 21, 2016
 *      Author: heber
 */


#include "BassoConfig.h"

#include <iostream>

#include <list>

#include <Eigen/SparseCore>
#include <unsupported/Eigen/ArpackSupport>

#include <boost/assign.hpp>

#include "TruncatedSVD.hpp"

#include "Log/Logging.hpp"

static TruncatedSVD::TypeToString_t initTypeToString()
{
	TruncatedSVD::TypeToString_t returnmap =
			boost::assign::list_of< TruncatedSVD::TypeToString_t::value_type >
	( TruncatedSVD::LargestMagnitude, "LM")
	( TruncatedSVD::SmallestMagnitude, "SM")
	( TruncatedSVD::LargestAlgebraic, "LA")
	( TruncatedSVD::SmallestAlgebraic, "SA")
	;
	return returnmap;
}

TruncatedSVD::TruncatedSVD(const MatrixType_t &_matrix) :
	TypeToString(initTypeToString()),
	matrix(_matrix)
{}

TruncatedSVD::~TruncatedSVD()
{}

inline
double sign(const double &_value)
{
	// zero is positive
	if (fabs(_value) > 0.)
		return (_value/fabs(_value));
	else
		return +1.;
}

template <class MatrixType, class ColumnType>
std::pair<double,double>
calculateLeftAndRightSign(
		const MatrixType &_matrix,
		const ColumnType &_U,
		const ColumnType &_V)
{
	std::pair<double,double> signs(std::make_pair(0.,0.));

	// left signs
	for (int j=0;j<_matrix.cols();++j) {
		const double tmp = _U.dot(_matrix.col(j));
		signs.first += sign(tmp)*tmp*tmp;
	}
	LOG(trace, "sign, left " << signs.first);
	// right signs
	for (int j=0;j<_matrix.rows();++j) {
		const double tmp = _V.dot(_matrix.row(j));
		signs.second += sign(tmp)*tmp*tmp;
	}
	LOG(trace, "sign, right " << signs.second);

	return signs;
}

/** Corrects the signs of the singular vectors.
 *
 * In an Singular Value Decomposition there is generally a ambiguity with
 * respect to sign because of the multiplicative setting of its three
 * components. While the singular values are non-negative, the signs of
 * U and V are arbitrary as long as USV resembles the original matrix.
 *
 * In a Technical report [Bro, Acer, Kolda '07] suggest to assign the
 * sign of each singular vector in such way that it points in the same
 * directions of the majority of data vectors it represents.
 *
 */
void TruncatedSVD::correctSigns(const bool _correlated)
{
	std::vector<std::pair<double,double> > signs;
	for (int k=0;k<S.rows();++k) {
		Eigen::MatrixXd residual(matrix);
		if (_correlated) {
			// this is only needed if data is correlated: this should
			// not be the case for a general low-rank factorization
			for (int m=0;m<S.rows();++m) {
				if (m==k)
					continue;
				residual -= U.col(m)*S(m)*V.col(m).transpose();
			}
			LOG(info, "Residual " << residual);
		}
		signs.push_back(
				calculateLeftAndRightSign(residual, U.col(k), V.col(k)));
	}

	for (int k=0;k<S.rows();++k) {
		if ((signs[k].first * signs[k].second) < 0.) {
			if (signs[k].first < signs[k].second)
				signs[k].first *= -1.;
			else
				signs[k].second *= -1.;
		}
		LOG(info, "signs left[" << k << "] " << signs[k].first << ", signs right[" << k << "] " << signs[k].second);

		U.col(k) *= sign(signs[k].first);
		V.col(k) *= sign(signs[k].second);
	}
}

std::string spellOutComputationInfo(const Eigen::ComputationInfo &_info)
{
	switch (_info) {
	case Eigen::Success:
		return std::string("Success");
		break;
	case Eigen::NumericalIssue:
		return std::string("NumericalIssue");
		break;
	case Eigen::NoConvergence:
		return std::string("NoConvergence");
		break;
	case Eigen::InvalidInput:
		return std::string("InvalidInput");
		break;
	default:
		return std::string("unknown");
		break;
	}
}


void TruncatedSVD::operator()(
		unsigned int k,
		enum EigenvalueType_t _type)
{
	const std::string desired_type = TypeToString[_type];

//	{
//		Eigen::ArpackGeneralizedSelfAdjointEigenSolver< MatrixType_t > solver(
//				matrix*matrix.transpose(),
//				k,
//				desired_type,
//				Eigen::ComputeEigenvectors,
//				BASSOTOLERANCE
//				);
//		LOG(info, "State o Convergence: " << spellOutComputationInfo(solver.info()));
//		LOG(info, "Number of iterations: " << solver.getNbrIterations());
//		LOG(info, "Number of converged eigenvalues: " << solver.getNbrConvergedEigenValues());
//
//		S = solver.eigenvalues();
//		for (int i=0;i<S.rows();++i)
//			S.coeffRef(i,0) = sqrt(S.coeffRef(i,0));
//		U = solver.eigenvectors();
//		solver.info();
//	}
	{
		Eigen::ArpackGeneralizedSelfAdjointEigenSolver< MatrixType_t > solver(
				matrix.transpose()*matrix,
				k,
				desired_type,
				Eigen::ComputeEigenvectors,
				BASSOTOLERANCE
				);
		LOG(debug, "State of Convergence: " << spellOutComputationInfo(solver.info()));
		LOG(debug, "Number of iterations: " << solver.getNbrIterations());
		LOG(debug, "Number of converged eigenvalues: " << solver.getNbrConvergedEigenValues());

		S = solver.eigenvalues();
		for (int i=0;i<S.rows();++i)
			S.coeffRef(i,0) = sqrt(S.coeffRef(i,0));
//		assert( (S - solver.eigenvalues()).norm() < sqrt(BASSOTOLERANCE));
		V = solver.eigenvectors();
	}

	// obtain U by 1/sigma_i * A * v_i
	std::list<int> zerovalues;
	U = Eigen::MatrixXd(matrix.rows(),S.rows());
	for (int i=0;i<S.rows();++i) {
		if (fabs(S.coeffRef(i,0)) > BASSOTOLERANCE) {
			U.col(i) = 1./S.coeffRef(i,0) * (matrix * V.col(i));
		} else {
			zerovalues.push_back(i);
		}
	}
	assert( zerovalues.empty() );
//	// all zero singular must be made orthogonal to all other vectors
//	assert( zerovalues.size() <= S.rows() );
//	for (; !zerovalues.empty();) {
//		const int index = zerovalues.back();
//		// create orthogonal vector
//		if (zerovalues.size()-1 == S.rows()) {
//			// there's just one other vector
//			// search for the only index of an orthogonal vector
//			int otherindex = -1;
//			for (int i=0;i<S.rows();++i) {
//				const std::list<int>::const_iterator iter =
//						std::find(zerovalues.begin(),
//								zerovalues.end(), i);
//				if (iter == zerovalues.end()) {
//					otherindex = i;
//					break;
//				}
//			}
//			assert( otherindex != -1 );
//			// create a linear independent vector
//			U.col(index).setZero();
//			if (fabs(U(0, otherindex) < BASSOTOLERANCE))
//				U(0,index) = 1.;
//			else if (fabs(U(1, otherindex) < BASSOTOLERANCE))
//				U(1,index) = 1.;
//			else {
//				U(0,index) = -U(1, otherindex);
//				U(1,index) = U(0, otherindex);
//			}
//			// and do Gram-Schmidt
//			const double expfactor = U.col(index).dot(U.col(otherindex));
//			assert( fabs(expfactor - 1.) > BASSOTOLERANCE );
//			U.col(index) = U.col(index) - U.col(index).dot(U.col(otherindex))*U.col(index);
//			U.col(index) *= 1./U.col(index).norm();
//		} else {
//			U.col(index).setZero();
//			for (int i=0;i<S.rows();++i) {
//				const std::list<int>::const_iterator iter =
//						std::find(zerovalues.begin(),
//								zerovalues.end(), i);
//				if (iter != zerovalues.end())
//					continue;
//				if (U.col(index).isZero())
//					U.col(index) = U.col(i);
//				else
//					U.col(index) = U.col(index).cross(U.col(i));
//			}
//		}
//	}

	assert( matrix.rows() == U.rows() );
	assert( matrix.cols() == V.rows() );
}
