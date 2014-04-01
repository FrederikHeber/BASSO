/*
 * MatrixIO.hpp
 *
 *  Created on: Apr 1, 2014
 *      Author: heber
 */

#ifndef MATRIXIO_HPP_
#define MATRIXIO_HPP_

#include <Eigen/Dense>
#include <iostream>
#include <sstream>
#include <string>


/** This namespace combines functions to read and write matrices and
 * vectors from the Eigen library
 */
namespace MatrixIO {

/** Output contents of \a m into stream \ost in Matlab-like format.
 *
 * \param ist output stream
 * \param m matrix to write in Matlab-like format
 * \return stream for concatenation
 */
std::istream & operator>>(std::istream &ist, Eigen::MatrixXd &m)
{
	if (ist.good()) {
		// parse first line with dimensions
		unsigned int rows = 0;
		unsigned int cols = 0;
		{
			std::string contents;
			getline(ist, contents);
			std::stringstream sstr(contents);
			sstr >> rows;
			sstr >> cols;
		}
		m.resize(rows, cols);

		// parse following lines
		for (unsigned int i=0; (i < rows) && ist.good(); ++i) {
			std::string contents;
			getline(ist, contents);
			std::stringstream sstr(contents);
			for (unsigned int j=0; (j < cols) && sstr.good(); ++j)
				sstr >> m(i,j);
		}
	}
	return ist;
}

std::istream & operator>>(std::istream &ist, Eigen::VectorXd &v)
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
