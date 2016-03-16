/*
 * detail.hpp
 *
 *  Created on: Jun 26, 2015
 *      Author: heber
 */

#ifndef MATRIXFACTORIZERBASE_HELPERS_DETAIL_HPP_
#define MATRIXFACTORIZERBASE_HELPERS_DETAIL_HPP_

#include "BassoConfig.h"

#include <Eigen/Dense>

#include "MatrixIO/MatrixIO.hpp"
#include "Minimizations/Elements/RepresentationAdvocate.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"

class MatrixFactorizerOptions;

namespace detail
{
	/** Calculates the residual of the Matrix Factorization problem.
	 *
	 * @param _data given matrix to factorize
	 * @param _spectral_matrix first factor
	 * @param _pixel_matrix second factor
	 * @return residual in l2 norm
	 */
	inline double calculateResidual(
			const Eigen::MatrixXd &_data,
			const Eigen::MatrixXd &_spectral_matrix,
			const Eigen::MatrixXd &_pixel_matrix
			)
	{
		const Eigen::MatrixXd difference_matrix =
				_data - _spectral_matrix * _pixel_matrix;
		return difference_matrix.norm();
	}

	/** Parses the options from the command-line.
	 *
	 * @param _argc argument count
	 * @param _argv argument array
	 * @param _opts on return parsed options
	 * @return returncode (0 - success)
	 */
	int parseOptions(
			int _argc,
			char ** _argv,
			MatrixFactorizerOptions &_opts);

	/** Parses the matrix to factorize from a file.
	 *
	 * @param _filename filename of matrix file
	 * @param _data on return parsed matrix
	 * @return returncode (0 - success)
	 */
	template <class T>
	int parseDataFile(
			const std::string &_filename,
			T &_data)
	{
		using namespace MatrixIO;

		if (!MatrixIO::parse(_filename, "data matrix", _data))
			return 255;
		return 0;
	}

	/** Parses the matrix factors from a file as starting point
	 *
	 * @param _filename filename of matrix file
	 * @param _data on return parsed matrix
	 * @param _name name of factor
	 * @return returncode (0 - success)
	 */
	int parseFactorFile(
			const std::string &_filename,
			Eigen::MatrixXd &_data,
			const std::string &_name);

	/** Constructs random starting matrices.
	 *
	 * @param _matrix matrix to construct
	 */
	void constructRandomMatrix(
			Eigen::MatrixXd &_matrix
			);

	/** Constructs zero starting matrices.
	 *
	 * @param _matrix matrix to construct
	 */
	void constructZeroMatrix(
			Eigen::MatrixXd &_matrix
			);

#ifdef MPI_FOUND

	//!> Unique tags for the MPI communication
	enum MPITags {
		InitialId,
		ColumnWise,
		MAX_MPITags
	};

#endif

} /* namespace detail */



#endif /* MATRIXFACTORIZERBASE_HELPERS_DETAIL_HPP_ */
