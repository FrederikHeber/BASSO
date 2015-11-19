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
	int parseDataFile(
			const std::string &_filename,
			Eigen::MatrixXd &_data);

	/** Constructs suitable starting matrices of either factor for the
	 * matrix factorization problem.
	 *
	 * @param _spectral_matrix first factor
	 * @param _pixel_matrix second factor
	 * @param _nonnegative whether to construct non-negative starting matrices
	 */
	void constructStartingMatrices(
			Eigen::MatrixXd &_spectral_matrix,
			Eigen::MatrixXd &_pixel_matrix,
			const bool _nonnegative
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
