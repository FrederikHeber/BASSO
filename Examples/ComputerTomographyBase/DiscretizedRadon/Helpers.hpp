/*
 * Helpers.hpp
 *
 *  Created on: Jul 22, 2015
 *      Author: heber
 */

#ifndef COMPUTERTOMOGRAPHYBASE_DISCRETIZEDRADON_HELPERS_HPP_
#define COMPUTERTOMOGRAPHYBASE_DISCRETIZEDRADON_HELPERS_HPP_

#include "BassoConfig.h"

#include <vector>

#include "ComputerTomographyBase/DiscretizedRadon/point_t.hpp"

namespace detail {

	/** Returns the center of the pixel indicated by the indices \a _i
	 * and \a _j, this is times 1 over h_x/h_y.
	 *
	 * @param _i index for x
	 * @param _j index for y
	 * @param _h_x length of a pixel in x direction
	 * @param _h_y length of a pixel in y direction
	 * @return pixel center coordinates
	 */
	point_t getPixelCenter(
			const unsigned int _i,
			const unsigned int _j,
			const double _h_x,
			const double _h_y
			);

	/** Function to generate a vector filled with omegas for a set of
	 * increasing angles.
	 *
	 * @param _num_angles number of angles
	 * @return vector with omegas for increasing angles [0:pi] in \a
	 * 		\a _num_angles steps
	 */
	std::vector< point_t > calculateOmegas(const unsigned int _num_angles);

	/** Functor structure for calculating omega vector from angle.
	 *
	 */
	struct calculateOmega
	{
		/** Function to calculate [cos(_phi), sin(_phi)].
		 *
		 * @param _phi angle
		 * @return [cos(_phi), sin(_phi)]
		 */
		point_t operator()(const double _phi) const;
	};

	/** Simple Functor to generate a vector of increasing angles.
	 *
	 */
	struct AngleIncrementer
	{
		/** Cstor for AngleIncrementer
		 *
		 * @param _delta delta width of the angle increment
		 */
		AngleIncrementer(const double _delta);

		/** Returns per call an angle always increased by \delta.
		 *
		 * @return increasing angle
		 */
		double operator()();

	private:
		//!> last given angle
		double angle;
		//!> step width for angle increment
		double delta;
	};

}; /* namespace detail */


#endif /* COMPUTERTOMOGRAPHYBASE_DISCRETIZEDRADON_HELPERS_HPP_ */
