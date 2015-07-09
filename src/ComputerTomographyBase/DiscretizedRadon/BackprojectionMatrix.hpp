/*
 * BackprojectionMatrix.hpp
 *
 *  Created on: Jul 8, 2015
 *      Author: heber
 */

#ifndef COMPUTERTOMOGRAPHYBASE_DISCRETIZEDRADON_BACKPROJECTIONMATRIX_HPP_
#define COMPUTERTOMOGRAPHYBASE_DISCRETIZEDRADON_BACKPROJECTIONMATRIX_HPP_

#include "BassoConfig.h"

#include <boost/filesystem/path.hpp>

#include "ComputerTomographyBase/DiscretizedRadon/point_t.hpp"

/** This class implements the back projection which is the adjoint operator
 * to the Radon transform.
 */
class BackprojectionMatrix
{
public:
	//!> constant for dimension of space
	enum { dim=2 };

	/** Constructor for class BackprojectionMatrix.
	 *
	 * Dimension of the vector of the right-hand side of the problem is
	 * \a _num_angles times \a _num_offsets. The dimension of the solution
	 * vector is \a _num_pixel_x \times \a _num_pixel_y.
	 *
	 * @param _num_pixel_x number of pixel in x direction of the image to reconstruct
	 * @param _num_pixel_y number of pixel in y direction of the image to reconstruct
	 * @param _num_beams number of beams (angles, offsets) in the measured data
	 * @param _num_offsets number of beams (angles, offsets) in the measured data
	 */
	BackprojectionMatrix(
			const unsigned int _num_pixel_x,
			const unsigned int _num_pixel_y,
			const unsigned int _num_angles,
			const unsigned int _num_offsets);

	/** Destructor for class DiscretizedRadonMatrix.
	 *
	 */
	~BackprojectionMatrix() {}

	/** Loads the matrix from a given \a _file.
	 *
	 * @param _file filename to load from
	 * @return 255 - failure, 0 - good
	 */
	int load(const boost::filesystem::path &_file);

	/** Getter for the constructed matrix.
	 *
	 * @return const ref to matrix object
	 */
	const Eigen::MatrixXd& getMatrix() const
	{ return matrix; }

private:

private:
	//!> number of pixel in x direction of the image to reconstruct
	const unsigned int num_pixel_x;
	//!> number of pixel in y direction of the image to reconstruct
	const unsigned int num_pixel_y;
	//!> number of beams (angles, offsets) in the measured data
	const unsigned int num_angles;
	//!> number of beams (angles, offsets) in the measured data
	const unsigned int num_offsets;

	//!> matrix containing discretized backprojection
	Eigen::MatrixXd matrix;
};



#endif /* COMPUTERTOMOGRAPHYBASE_DISCRETIZEDRADON_BACKPROJECTIONMATRIX_HPP_ */
