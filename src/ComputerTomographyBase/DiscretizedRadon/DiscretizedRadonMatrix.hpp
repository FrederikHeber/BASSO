/*
 * DiscretizedRadonMatrix.hpp
 *
 *  Created on: Jul 6, 2015
 *      Author: heber
 */

#ifndef COMPUTERTOMOGRAPHYBASE_DISCRETIZEDRADONMATRIX_HPP_
#define COMPUTERTOMOGRAPHYBASE_DISCRETIZEDRADONMATRIX_HPP_

#include "BassoConfig.h"

#include <iosfwd>
#include <set>
#include <vector>

#include <boost/array.hpp>
#include <boost/filesystem/path.hpp>

#include <Eigen/SparseCore>

#include "ComputerTomographyBase/DiscretizedRadon/point_t.hpp"

class DiscretizedRadonMatrixUnitTest;

/** This class reconstructs a matrix resembling the discretized
 * Radon transformation of some given data.
 *
 * We assume the domain is [-1,1]^dim.
 *
 */
class DiscretizedRadonMatrix
{
	//!> grant unit test access to private functions
	friend class DiscretizedRadonMatrixUnitTest;
public:
	/** Constructor for class DiscretizedRadonMatrix.
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
	DiscretizedRadonMatrix(
			const unsigned int _num_pixel_x,
			const unsigned int _num_pixel_y,
			const unsigned int _num_angles,
			const unsigned int _num_offsets);

	/** Destructor for class DiscretizedRadonMatrix.
	 *
	 */
	~DiscretizedRadonMatrix() {}

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
	const Eigen::SparseMatrix<double, Eigen::RowMajor>& getMatrix() const
	{ return matrix; }

	//!> typedef for a vector of intersection points
	typedef std::set<point_t> intersections_t;

	struct pixel_t : public boost::array<unsigned int, point_t::dim>
	{
	};

	typedef std::vector< pixel_t > pixels_t;

	void setMatrix(const Eigen::SparseMatrix<double, Eigen::RowMajor>& _matrix)
	{ matrix = _matrix; }

private:
	/** Calculates intersection points between a line and all pixel
	 * boundaries, i.e. the discretized grid lines.
	 *
	 * @param _phi angle of line with respect to x axis
	 * @param _s offset of line with respect to x axis
	 * @return sorted list of intersection points
	 */
	intersections_t calculatePixelBoundaryIntersections(
			const double _phi,
			const double _s) const;

	/** Removes adjacent points that are closer than BASSOTOLERANCE to
	 * each other.
	 *
	 * @param _intersections intersections to scan for identical points
	 */
	void removeIdenticalAdjacentPoints(
			intersections_t &_intersections) const;

	/** Removes illegal points that translate to pixels outside the matrix's
	 * array sizes.
	 *
	 * @param _intersections intersections to scan for illegal points
	 */
	void removeIllegalPoints(
			intersections_t &_intersections) const;

	/** Flips pixels that are on the upper or right boundary.
	 *
	 * @param _angle angle to check what beam we have
	 * @param _offset offset to check what beam we have
	 * @param _pixels pixels to modify
	 */
	void correctBoundaryPixels(
			const unsigned int _angle,
			const unsigned int _offset,
			pixels_t &_pixels) const;

	/** Calculates for each intersection point the corresponding pixel
	 * indices.
	 *
	 * @param _intersections
	 * @return
	 */
	pixels_t PointsToPixels(
			const intersections_t &_intersections) const;

	/** Calculates the length of the adjacent intersection associated
	 * with a specific pixel and adds the length onto the \a matrix entry.
	 *
	 * @param _angle angle index
	 * @param _offset offset index
	 * @param _pixels set of pixels associated with adjacent intersections
	 * @param _intersections set of sorted intersections
	 */
	void calculateLengths(
			unsigned int _angle,
			unsigned int _offset,
			const pixels_t &_pixels,
			const intersections_t &_intersections);

private:
	//!> number of pixel in x direction of the image to reconstruct
	const unsigned int num_pixel_x;
	//!> number of pixel in y direction of the image to reconstruct
	const unsigned int num_pixel_y;
	//!> number of beams (angles, offsets) in the measured data
	const unsigned int num_angles;
	//!> number of beams (angles, offsets) in the measured data
	const unsigned int num_offsets;

	//!> matrix containing discretized radon transform
	Eigen::SparseMatrix<double, Eigen::RowMajor> matrix;

	bool debugflag;
};

/** Output operator for class pixel_t.
 *
 * @param _out output stream
 * @param _pixel pixel to print
 * @return stream for concatenation
 */
std::ostream& operator<<(
		std::ostream &ost,
		const DiscretizedRadonMatrix::pixel_t &_pixel);

#endif /* COMPUTERTOMOGRAPHYBASE_DISCRETIZEDRADONMATRIX_HPP_ */
