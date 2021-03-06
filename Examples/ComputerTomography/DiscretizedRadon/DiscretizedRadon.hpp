/*
 * DiscretizedRadon.hpp
 *
 *  Created on: Jul 8, 2015
 *      Author: heber
 */

#ifndef COMPUTERTOMOGRAPHYBASE_DISCRETIZEDRADON_DISCRETIZEDRADON_HPP_
#define COMPUTERTOMOGRAPHYBASE_DISCRETIZEDRADON_DISCRETIZEDRADON_HPP_

#include "BassoConfig.h"

#include "ComputerTomography/DiscretizedRadon/DiscretizedRadonMatrix.hpp"

#include "Minimizations/Mappings/Mapping.hpp"

/** This implements the RadonTransform operator.
 *
 * The problem is that we need a specific adjoint operator, namely the
 * backprojection.
 */
class DiscretizedRadon : public Mapping
{
public:
	/** Constructor for class DiscretizedRadonTransform.
	 *
	 * Dimension of the vector of the right-hand side of the problem is
	 * \a _num_angles times \a _num_offsets. The dimension of the solution
	 * vector is \a _num_pixel_x \times \a _num_pixel_y.
	 *
	 * @param _SourceSpaceRef source space
	 * @param _TargetSpaceRef target space
	 * @param _num_pixel_x number of pixel in x direction of the image to reconstruct
	 * @param _num_pixel_y number of pixel in y direction of the image to reconstruct
	 * @param _num_beams number of beams (angles, offsets) in the measured data
	 * @param _num_offsets number of beams (angles, offsets) in the measured data
	 */
	DiscretizedRadon(
			const NormedSpace_weakptr_t _SourceSpaceRef,
			const NormedSpace_weakptr_t _TargetSpaceRef,
			const unsigned int _num_pixel_x,
			const unsigned int _num_pixel_y,
			const unsigned int _num_angles,
			const unsigned int _num_offsets);

	//!> expose overloaded operator method from base class
	using Mapping::operator();

	void operator()(
				const SpaceElement_ptr_t &_sourceelement,
				SpaceElement_ptr_t &_destelement
				) const;

	const Mapping_ptr_t getAdjointMapping() const;

	/** Returns the number of times the operator() was called.
	 *
	 * @return number of calls
	 */
	const unsigned int getCount() const
	{ return MatrixVectorProductCounts; }

	/** Returns the total runtime the program spent so far on
	 * its operator().
	 *
	 * \warning Due to Eigen's lazy evaluation we cannot really
	 * measure the time required on matrix-vector products.
	 *
	 * @return runtime summed over all calls
	 */
	const boost::chrono::nanoseconds getTiming() const
	{ return boost::chrono::nanoseconds(0); }

	/** Internal setter for the adjoint mapping
	 *
	 * When we create the adjoint mapping for the first time, the
	 * adjoint's adjoint mapping is the original one and hence must
	 * not be created.
	 *
	 * @param _adjoint adjoint instance to set
	 */
	void setAdjointMapping(const Mapping_weakptr_t &_adjoint);

	/** Getter for the internal matrix.
	 *
	 * @return const ref to matrix
	 */
	DiscretizedRadonMatrix& get()
	{ return radon_matrix; }

private:
	//!> internal pointer to adjoint
	const Mapping_weakptr_t AdjointLinearMapping;

	//!> internal radon matrix
	DiscretizedRadonMatrix radon_matrix;

	//!> typedef for counting mvps
	typedef unsigned int Count_t;

	//!> Counts the number of calculated matrix-vector products
	mutable Count_t MatrixVectorProductCounts;

	//!> counts the time used for matrix-vector products
	mutable boost::chrono::nanoseconds MatrixVectorProductTimings;
};



#endif /* COMPUTERTOMOGRAPHYBASE_DISCRETIZEDRADON_DISCRETIZEDRADON_HPP_ */
