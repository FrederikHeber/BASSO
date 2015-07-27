/*
 * MatrixIO.hpp
 *
 *  Created on: Apr 1, 2014
 *      Author: heber
 */

#ifndef MATRIXIO_HPP_
#define MATRIXIO_HPP_

#include "BassoConfig.h"

#include <iostream>
#include <sstream>
#include <string>

#include <boost/filesystem.hpp>

#include <Eigen/Dense>
#include <Eigen/SparseCore>

#include "MatrixIO/MatrixIOExceptions.hpp"

/** This namespace combines functions to read and write matrices and
 * vectors from the Eigen library
 */
namespace MatrixIO {

/** Output contents of \a m into stream \ost in Matlab-like format.
 *
 * \param ost output stream
 * \param m matrix to write in Matlab-like format
 * \return stream for concatenation
 */
inline std::ostream & operator<<(
		std::ostream &ost,
		const Eigen::MatrixXd &m)
{
	if (ost.good()) {
		const Eigen::MatrixXd::Index rows = m.innerSize();
		const Eigen::MatrixXd::Index cols = m.outerSize();
		ost << rows << " " << cols << std::endl;
		Eigen::MatrixXd::Index row = 0;
		for (; (row < rows) && ost.good(); ++row ) {
			Eigen::MatrixXd::Index col = 0;
			for (; (col < cols) && ost.good(); ++col )
				ost << m(row, col) << " ";
			ost << std::endl;
			if (col != cols)
				throw MatrixIOStreamEnded_exception();
		}
		if (row != rows)
			throw MatrixIOStreamEnded_exception();
	}
	return ost;
}

/** Output contents of \a m into stream \ost in Matlab-like format.
 *
 * \param ost output stream
 * \param m matrix to write in Matlab-like format
 * \return stream for concatenation
 */
template <typename type,
	int options>
inline std::ostream & operator<<(
		std::ostream &ost,
		const Eigen::SparseMatrix<type, options> &m)
{
	typedef typename Eigen::SparseMatrix<type, options>::Index index_t;
	if (ost.good()) {
		const index_t rows = m.innerSize();
		const index_t cols = m.outerSize();
		ost << rows << " " << cols << std::endl;
		index_t row = 0;
		for (; (row < rows) && ost.good(); ++row ) {
			index_t col = 0;
			for (; (col < cols) && ost.good(); ++col )
				ost << m.coeff(row, col) << " ";
			ost << std::endl;
			if (col != cols)
				throw MatrixIOStreamEnded_exception();
		}
		if (row != rows)
			throw MatrixIOStreamEnded_exception();
	}
	return ost;
}

/** Output contents of \a v into stream \ost in Matlab-like format.
 *
 * \param ost output stream
 * \param v vector to write in Matlab-like format
 * \return stream for concatenation
 */
inline std::ostream & operator<<(std::ostream &ost, const Eigen::VectorXd &v)
{
	if (ost.good()) {
		const Eigen::MatrixXd::Index rows = v.innerSize();
		ost << rows << " " << 1 << std::endl;
		Eigen::MatrixXd::Index row = 0;
		for (; (row < rows) && ost.good(); ++row )
			ost << v(row) << std::endl;
		if (row != rows)
			throw MatrixIOStreamEnded_exception();
	}
	return ost;
}


/** Parse a Matlab-like formatted stream \a ist into the matrix \a m.
 *
 * @param ist input stream
 * @param m matrix to parse
 * @return stream for concatenation
 */
inline std::istream & operator>>(
		std::istream &ist,
		Eigen::MatrixXd &m)
{
	if (ist.good()) {
		// parse first line with dimensions
		Eigen::MatrixXd::Index rows = 0;
		Eigen::MatrixXd::Index cols = 0;
		{
			std::string contents;
			getline(ist, contents);
			std::stringstream sstr(contents);
			sstr >> rows;
			sstr >> cols;
		}
		m.resize(rows, cols);

		// parse following lines
		Eigen::MatrixXd::Index i=0;
		for (; (i < rows) && ist.good(); ++i) {
			std::string contents;
			getline(ist, contents);
			std::stringstream sstr(contents);
			Eigen::MatrixXd::Index j=0;
			for (; (j < cols) && sstr.good(); ++j)
				sstr >> m.coeffRef(i,j);
			if (j != cols)
				throw MatrixIOStreamEnded_exception();
		}
		if (i != rows)
			throw MatrixIOStreamEnded_exception();
	}
	return ist;
}

//!> typedef for a pair of dimensions, representing row and column counts
typedef std::pair<Eigen::MatrixXd::Index,Eigen::MatrixXd::Index> dims_t;

/** Determines the dimensions for a vector or matrix object contained
 * in the stream \a ist.
 *
 * \note we do not change the position inside the stream but the stream
 * position must be right at the start of the object
 *
 * @param ist stream to parse dimensions of object from
 * @return pair with rows and columns of oject
 */
inline dims_t getDimensions(std::istream &ist)
{
	std::pair<Eigen::MatrixXd::Index,Eigen::MatrixXd::Index> result;

	size_t position = ist.tellg();
	if (ist.good()) {
		// parse first line with dimensions
		result.first = 0;
		result.second = 0;
		{
			std::string contents;
			getline(ist, contents);
			std::stringstream sstr(contents);
			sstr >> result.first;
			sstr >> result.second;
		}
	}
	ist.seekg(position);

	return result;
}

/** Checks whether the string represents the path of a present file.
 *
 * @param _name name of file
 * @return true file present, false - else
 */
inline bool isPresentFile(const std::string &_name)
{
	boost::filesystem::path matrix_file(_name);
	return boost::filesystem::exists(matrix_file);
}

/** Parse a Matlab-like formatted stream \a ist into the matrix \a m.
 *
 * @param ist input stream
 * @param m matrix to parse
 * @return stream for concatenation
 */
template <
	typename type,
	int options>
inline std::istream & operator>>(
		std::istream &ist,
		Eigen::SparseMatrix<type, options> &m)
{
	typedef typename Eigen::SparseMatrix<type, options>::Index index_t;
	if (ist.good()) {
		// parse first line with dimensions
		index_t rows = 0;
		index_t cols = 0;
		{
			std::string contents;
			getline(ist, contents);
			std::stringstream sstr(contents);
			sstr >> rows;
			sstr >> cols;
		}
		m.resize(rows, cols);

		// parse following lines
		index_t i=0;
		for (; (i < rows) && ist.good(); ++i) {
			std::string contents;
			getline(ist, contents);
			std::stringstream sstr(contents);
			index_t j=0;
			for (; (j < cols) && sstr.good(); ++j)
				sstr >> m.coeffRef(i,j);
			if (j != cols)
				throw MatrixIOStreamEnded_exception();
		}
		if (i != rows)
			throw MatrixIOStreamEnded_exception();
	}
	return ist;
}

inline std::istream & operator>>(std::istream &ist, Eigen::VectorXd &v)
{
	Eigen::MatrixXd m;
	ist >> m;
	v.resize(m.outerSize(), 1);
	v = m.leftCols(1);
	return ist;
}

/** Parse a Matlab-like formatted stream \a ist into the matrix \a m.
 *
 * \param ist input stream with Matlab-like format
 * \param m matrix to set
 * \return stream for concatenation
 */
//std::ostream & operator<<(std::ostream &ost, const Eigen::MatrixXd &m)
//{
//	return ost;
//}

};

#endif /* MATRIXIO_HPP_ */
