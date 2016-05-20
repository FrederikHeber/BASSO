/*
 * TruncatedSVD.cpp
 *
 *  Created on: Mar 21, 2016
 *      Author: heber
 */


#include "BassoConfig.h"

#include <iostream>

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


void TruncatedSVD::operator()(
		unsigned int k,
		enum EigenvalueType_t _type)
{
	const std::string desired_type = TypeToString[_type];

	{
		Eigen::ArpackGeneralizedSelfAdjointEigenSolver< MatrixType_t > solver(
				matrix*matrix.transpose(),
				k,
				desired_type,
				Eigen::ComputeEigenvectors,
				BASSOTOLERANCE
				);

		S = solver.eigenvalues();
		for (int i=0;i<S.rows();++i)
			S.coeffRef(i,0) = sqrt(S.coeffRef(i,0));
		U = solver.eigenvectors();
	}
	{
		Eigen::ArpackGeneralizedSelfAdjointEigenSolver< MatrixType_t > solver(
				matrix.transpose()*matrix,
				k,
				desired_type,
				Eigen::ComputeEigenvectors,
				BASSOTOLERANCE
				);

//		assert( (S - solver.eigenvalues()).norm() < sqrt(BASSOTOLERANCE));
		V = solver.eigenvectors();
	}
	assert( matrix.rows() == U.rows() );
	assert( matrix.cols() == V.rows() );

	/** Basically, we do here something similar to [Bro, Acer, Kolda '07],
	 * \sa TruncatedSVD::correctSigns().
	 *
	 * However, in our case the problem is even larger: The signs are not
	 * even set in such a way as to fulfill the above condition of resembling
	 * the original data matrix. This is because both are obtained via a
	 * truncated eigenvalue calculation separately. Hence, singular vectors
	 * in U and V know nothing about each other, especially not with respect
	 * to their sign, as the overall sign of an eigenvector is arbitrary.
	 *
	 * Hence, we use the ansatz of the Technical Report of comparing the
	 * singular vector with the majority of the data vectors and flipping
	 * the sign if the majority ends up pointing in the other direction.
	 * Afterwards, the correctSigns() may still be applied.
	 *
	 */
	for (int k=0;k<S.rows();++k) {
		const Eigen::MatrixXd fullmatrix(matrix);
		const std::pair<double, double> signs(
				calculateLeftAndRightSign(fullmatrix, U.col(k), V.col(k)));

		U.col(k) *= sign(signs.first);
		V.col(k) *= sign(signs.second);
	}
}
