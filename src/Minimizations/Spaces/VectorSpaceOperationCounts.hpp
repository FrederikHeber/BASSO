/*
 * VectorSpaceOperationCounts.hpp
 *
 *  Created on: Dec 12, 2014
 *      Author: heber
 */

#ifndef MINIMIZATIONS_SPACES_VECTORSPACEOPERATIONCOUNTS_HPP_
#define MINIMIZATIONS_SPACES_VECTORSPACEOPERATIONCOUNTS_HPP_

#include "BassoConfig.h"

#include <boost/chrono.hpp>
#include <utility>

/** This simple struct contains counts and timings for all operations
 * in a vector space.
 */
struct VectorSpaceOperationCounts
{
	//!> typedef for the count type
	typedef int Count_t;
	//!> typedef for the time-keeping type
	typedef boost::chrono::nanoseconds Timing_t;
	//!> typedef for the combined count and time-keeping type
	typedef std::pair<Count_t, Timing_t> CountTiming_t;

	VectorSpaceOperationCounts() :
		ElementCreation(std::make_pair(0, boost::chrono::nanoseconds(0))),
		VectorAddition(std::make_pair(0, boost::chrono::nanoseconds(0))),
		VectorAssignment(std::make_pair(0, boost::chrono::nanoseconds(0))),
		VectorComparison(std::make_pair(0, boost::chrono::nanoseconds(0))),
		VectorModification(std::make_pair(0, boost::chrono::nanoseconds(0))),
		VectorMultiplication(std::make_pair(0, boost::chrono::nanoseconds(0))),
		ScalarVectorMultiplication(std::make_pair(0, boost::chrono::nanoseconds(0)))
	{}

	//!> time and count for the creation of vectors
	CountTiming_t ElementCreation;
	//!> time and count for modifying a vector componentwise: O(n)
	CountTiming_t VectorModification;
	//!> time and count for comparing a vector to another or to constant: O(n)
	CountTiming_t VectorComparison;
	//!> time and count for assigning a vector: O(n)
	CountTiming_t VectorAssignment;
	//!> time and count for adding/subtracting of vectors: O(n)
	CountTiming_t VectorAddition;
	//!> time and count for the vector vector multiplication: O(n)
	CountTiming_t VectorMultiplication;
	//!> time and count for the vector multiplication with a scalar: O(n)
	CountTiming_t ScalarVectorMultiplication;

	/** This is a convenience struct that measures the time that
	 * passed by during its existence.
	 *
	 * In order to measure the time spent in a function, we only
	 * require to instantiate such a TimeKeeper object prior to
	 * when the time measuring should start (this allows excluding
	 * certain initial operations). The time measurement stops when
	 * the function leaves the scope and the TimeKeeper instance is
	 * destroyed.
	 *
	 * \code
	 * const double examplefunction(
	 *     const SpaceElement_ptr_t &_a,
	 *     const SpaceElement_ptr_t &_b)
	 * {
	 *   TimeKeeper(VectorMultiplication);
	 *   return _a * _b;
	 * }
	 * \endcode
	 *
	 */
	struct TimeKeeper {
		/** Constructor of TimeKeeper.
		 *
		 * Increments count and initializes time keeping
		 *
		 * @param _optype reference to one of the counttimes above
		 */
		TimeKeeper(CountTiming_t &_optype) :
			optype(_optype),
			timing_start(boost::chrono::high_resolution_clock::now())
		{
			// increment count
			++optype.first;
		}

		/** Destructor of TimeKeeper.
		 *
		 * Finalizes time keeping.
		 */
		~TimeKeeper()
		{
			const boost::chrono::high_resolution_clock::time_point timing_end =
					boost::chrono::high_resolution_clock::now();
			optype.second += timing_end - timing_start;
		}

	private:
		//!> stores the optype to write results to
		CountTiming_t &optype;
		//!> contains the starting time
		const boost::chrono::high_resolution_clock::time_point timing_start;
	};
};



#endif /* MINIMIZATIONS_SPACES_VECTORSPACEOPERATIONCOUNTS_HPP_ */
