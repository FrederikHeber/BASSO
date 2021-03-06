/*
 * Eigen_matrix_serizalization.hpp
 *
 *  Created on: Jun 29, 2015
 *      Author: heber
 */

#ifndef MINIMIZATIONS_ELEMENTS_EIGEN_MATRIX_SERIALIZATION_HPP_
#define MINIMIZATIONS_ELEMENTS_EIGEN_MATRIX_SERIALIZATION_HPP_

#include "BassoConfig.h"

#include <Eigen/Dense>

#include <boost/serialization/array.hpp>
#include <boost/serialization/split_free.hpp>

/** Answer from
 * http://stackoverflow.com/questions/12851126/serializing-eigens-matrix-using-boost-serialization
 */
namespace boost
{
	namespace serialization {
		template<class Archive, typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
		inline void load(
			Archive & ar,
			Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> & t,
			const unsigned int file_version
		)
		{
			size_t rows;
			size_t cols;
			ar >> rows;
			ar >> cols;
			if ((rows != -1) && (cols != -1)) {
				if ((rows != t.rows()) || (cols != t.cols()))
					t.resize( rows, cols );
				ar >> boost::serialization::make_array(t.data(), t.size());
			}
		}

		template<class Archive, typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
		inline void save(
			Archive & ar,
			const Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> & t,
			const unsigned int file_version
		)
		{
			size_t rows = t.rows();
			size_t cols = t.cols();
			ar << rows;
			ar << cols;
			if ((rows != -1) && (cols != -1))
				ar << boost::serialization::make_array(t.data(), t.size());
		}

		template<class Archive, typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
		inline void serialize(
		    Archive & ar,
			Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> & t,
		    const unsigned int file_version
		){
		    split_free(ar, t, file_version);
		}

	} /* namespace serialization */
} /* namespace boost */

#endif /* MINIMIZATIONS_ELEMENTS_EIGEN_MATRIX_SERIALIZATION_HPP_ */
