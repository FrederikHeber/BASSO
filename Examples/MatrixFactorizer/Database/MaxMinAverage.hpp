/*
 * MaxMinAverage.hpp
 *
 *  Created on: Nov 24, 2015
 *      Author: heber
 */

#ifndef MAXMINAVERAGE_HPP_
#define MAXMINAVERAGE_HPP_

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

/** This structure extracts minimum, maximum, and  average with variance out of
 * a given vector of values.
 */
template <class T>
struct MaxMinAverage
{
	//!> typedef for the vector of values
	typedef std::vector<T> values_t;

	MaxMinAverage(const values_t &_valuevector) :
		minimum(std::numeric_limits<T>::max()),
		maximum(std::numeric_limits<T>::min()),
		average(0),
		variance(0)
	{
		for (typename values_t::const_iterator valueiter = _valuevector.begin();
				valueiter != _valuevector.end(); ++valueiter) {
			maximum = std::max(maximum, *valueiter);
			minimum = std::min(minimum, *valueiter);
			average += *valueiter;
			variance += (*valueiter)*(*valueiter);
		}
		average = average/(double)_valuevector.size();
		variance = variance/(double)_valuevector.size();
		variance = sqrt(variance-average*average);
	}

	T minimum;
	T maximum;
	T average;
	T variance;
};


#endif /* MAXMINAVERAGE_HPP_ */
