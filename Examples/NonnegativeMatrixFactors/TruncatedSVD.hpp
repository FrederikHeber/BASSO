/*
 * TruncatedSVD.hpp
 *
 *  Created on: Mar 21, 2016
 *      Author: heber
 */

#ifndef TRUNCATEDSVD_HPP_
#define TRUNCATEDSVD_HPP_

#include "BassoConfig.h"

#include <map>

/** This class wraps calls to Eigen's unsupported ArpackSupport
 * module for computing a truncated SVD.
 */
class TruncatedSVD
{
public:
	typedef Eigen::SparseMatrix<double> MatrixType_t;

	TruncatedSVD(const MatrixType_t &_matrix);
	~TruncatedSVD();

	//!> enumeration of the k eigenvalues of which type: largest, smallest, ...
	enum EigenvalueType_t {
		LargestMagnitude,
		SmallestMagnitude,
		LargestAlgebraic,
		SmallestAlgebraic,
		MAX_EigenvalueType
	};

	typedef std::map<EigenvalueType_t, std::string> TypeToString_t;
	//!> translation map from enumerated type to string to pass to ARPACK
	TypeToString_t TypeToString;

	void operator()(
			unsigned int k,
			enum EigenvalueType_t _type);

	const Eigen::MatrixXd& getSingularvalues() const
	{ return S; }

	const Eigen::MatrixXd& getLeftSingularvectors() const
	{ return U; }

	const Eigen::MatrixXd& getRightSingularvectors() const
	{ return V; }

	void correctSigns(const bool _correlated = false);

private:
	// ref to matrix to perform truncated svd on
	const MatrixType_t &matrix;
	//!> contains computed eigenvalues after operator()
	Eigen::MatrixXd S;
	//!> contains computed left eigenvectors after operator()
	Eigen::MatrixXd U;
	//!> contains computed right eigenvectors after operator()
	Eigen::MatrixXd V;
};



#endif /* TRUNCATEDSVD_HPP_ */
