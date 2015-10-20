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

#include "ComputerTomography/DiscretizedRadon/point_t.hpp"

namespace detail {

	/** Precalculates the pixel centers for on direction.
	 *
	 * @param _h length of a pixel in the respective direction
	 * @param _num_pixel number of pixel in the respective direction
	 * @return pixel center component for respective direction
	 */
	std::vector<double> calculatePixelCenter(
			const double _h,
			const unsigned int _num_pixel);

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
