/*
 * RelativeShrinkageCoefficient.cpp
 *
 *  Created on: Apr 25, 2017
 *      Author: heber
 */

#include <functional>
#include <map>

#include "RelativeShrinkageCoefficient.hpp"

#include "Log/Logging.hpp"
#include "Math/Helpers.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"

const double
RelativeShrinkageCoefficient::get(
		const SpaceElement_ptr_t &_x,
		const double _lambda
		)
{
	/// we use algorithm 1 given in [Schoepfer 2012]
	if (_x->isZero(BASSOTOLERANCE))
		return 0.;

	size_t active_set_size = 0;
	double l1norm_sum = 0.;
	double coefficient = 0.;
	typedef std::vector<size_t> active_set_t;
	active_set_t active_set;

	// create a map of largest coefficients
	typedef std::multimap<double, size_t> largest_coeffs_map_t;
	largest_coeffs_map_t largest_coeffs_map;
	for (size_t i=0;i<_x->getSpace()->getDimension();++i)
		largest_coeffs_map.insert( std::make_pair(fabs((*_x)[i]), i) );
	// loop over element's components and check against current coefficient
	do {
		// extract largest coefficient from map
		const largest_coeffs_map_t::reverse_iterator iter = largest_coeffs_map.rbegin();
		const std::pair<double,size_t> maxCoeff = *iter;
		largest_coeffs_map.erase(--(iter.base()));
		// compare with coefficient
//		LOG(info, "Is current max " << maxCoeff.first
//			<< " greater than coefficient " << coefficient << "?");
		if (coefficient > maxCoeff.first + BASSOTOLERANCE)
			break;
		// add to l1norm sum and to active set
//		LOG(info, "Adding " << maxCoeff.second << " to active_set.";
		l1norm_sum += maxCoeff.first;
		++active_set_size;
		active_set.push_back(maxCoeff.second);
		// and recalculate coefficient
		coefficient = l1norm_sum/((double)active_set_size + _lambda);
	} while (!largest_coeffs_map.empty());

//	// sensibility checks
//	// all component in set exceed coefficient
//	std::sort(active_set.begin(), active_set.end());
////	{
////		std::stringstream output;
////		std::copy(
////				active_set.begin(), active_set.end(),
////				std::ostream_iterator<int>(output, ";"));
////		LOG(info, "Final active_set is " << output.str();
////	}
//	for (active_set_t::const_iterator iter = active_set.begin();
//			iter != active_set.end(); ++iter) {
////		LOG(info, "Is |_x[" << *iter << "]| = " <<  fabs((*_x)[*iter]) << " >= " << coefficient << "?");
//		assert( fabs((*_x)[*iter]) + BASSOTOLERANCE >= coefficient );
//	}
//	// all other components don't
//	active_set_t complement_set(_x->getSpace()->getDimension(), 0);
//	std::generate(complement_set.begin(), complement_set.end(), Helpers::unique_number());
//	active_set_t::iterator eraseiter =
//			std::set_difference(
//				complement_set.begin(), complement_set.end(),
//				active_set.begin(), active_set.end(),
//				complement_set.begin());
//	for (active_set_t::const_iterator iter = complement_set.begin();
//			iter != eraseiter; ++iter) {
////		LOG(info, "Is |_x[" << *iter << "]| = " <<  fabs((*_x)[*iter]) << " < " << coefficient << "?");
//		assert( fabs((*_x)[*iter]) < coefficient + BASSOTOLERANCE);
//	}

	return coefficient;
}
