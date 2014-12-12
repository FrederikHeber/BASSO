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
#include <boost/fusion/sequence.hpp>
#include <boost/mpl/at.hpp>
#include <utility>
#include <iosfwd>

#include "Minimizations/Spaces/OperationCountMap.hpp"

#ifdef USE_TIMINGS
#define TIMEKEEPER(type) \
	VectorSpaceOperationCounts::TimeKeeper(type)
#else
#define TIMEKEEPER(type)
#endif

namespace VectorSpaceOperations {

	/** Allows access to a certain timing object.
	 *
	 * @return ref to timing object
	 */
	template<class _type>
	CountTiming_t& getCountTiming(
			OperationCountMap_t &_opcounts
			)
	{
		return boost::fusion::at_key<_type>(_opcounts);
	}

	/** Allows access to a certain timing object.
	 *
	 * @return ref to timing object
	 */
	template<class _type>
	const CountTiming_t& getCountTiming(
			const OperationCountMap_t &_opcounts
			)
	{
		return boost::fusion::at_key<_type>(_opcounts);
	}

	template<class _type>
	struct TypeToString {
		const std::string operator()() const
		{ return ""; }
	};

	template <>
	struct TypeToString<VectorSpaceOperations::ElementCreation> {
		const std::string operator()() const
				{ return "ElementCreation"; }
	};
	template <>
	struct TypeToString<VectorSpaceOperations::ScalarVectorMultiplication> {
		const std::string operator()() const
		{ return "ScalarVectorMultiplication"; }
	};
	template <>
	struct TypeToString<VectorSpaceOperations::VectorAddition> {
		const std::string operator()() const
		{ return "VectorAddition"; }
	};
	template <>
	struct TypeToString<VectorSpaceOperations::VectorAssignment> {
		const std::string operator()() const
		{ return "VectorAssignment"; }
	};
	template <>
	struct TypeToString<VectorSpaceOperations::VectorComparison> {
		const std::string operator()() const
		{ return "VectorComparison"; }
	};
	template <>
	struct TypeToString<VectorSpaceOperations::VectorModification> {
		const std::string operator()() const
		{ return "VectorModification"; }
	};
	template <>
	struct TypeToString<VectorSpaceOperations::VectorMultiplication> {
		const std::string operator()() const
		{ return "VectorMultiplication"; }
	};
	template <>
	struct TypeToString<VectorSpaceOperations::VectorNorm> {
		const std::string operator()() const
		{ return "VectorNorm"; }
	};
};

/** This simple struct contains counts and timings for all operations
 * in a vector space.
 */
struct VectorSpaceOperationCounts
{
	OperationCountMap_t instance;

	/** Returns the total number of O(n) operations.
	 *
	 * \note This excludes ElementCreation and this is not necessarily
	 * O(n) but simply an allocation of O(n) memory. Also, VectorNorm
	 * is not included as this is here with Banach Spaces interesting
	 * on its own.
	 *
	 * \sa getTotalElementCreationCounts() and getTotalVectorNormCounts()
	 *
	 * @return total number of operations called
	 */
	VectorSpaceOperations::Count_t getTotalLinearCounts() const
	{ return boost::fusion::at_key<VectorSpaceOperations::ScalarVectorMultiplication>(instance).first
			+boost::fusion::at_key<VectorSpaceOperations::VectorAddition>(instance).first
			+boost::fusion::at_key<VectorSpaceOperations::VectorAssignment>(instance).first
			+boost::fusion::at_key<VectorSpaceOperations::VectorComparison>(instance).first
			+boost::fusion::at_key<VectorSpaceOperations::VectorModification>(instance).first
			+boost::fusion::at_key<VectorSpaceOperations::VectorMultiplication>(instance).first
			+boost::fusion::at_key<VectorSpaceOperations::VectorNorm>(instance).first;
	}

	VectorSpaceOperations::Count_t getTotalConstantCounts() const
	{ return boost::fusion::at_key<VectorSpaceOperations::ElementCreation>(instance).first; }

	/** Returns the total time spent in O(n) operations.
	 *
	 * \note This excludes ElementCreation and this is not necessarily
	 * O(n) but simply an allocation of O(n) memory. Also, VectorNorm
	 * is not included as this is here with Banach Spaces interesting
	 * on its own.
	 *
	 * \warning the times may not actually be meaningful due to lazy
	 * evaluation of \b Eigen. I.e. an operation may occur outside
	 * the scope of the time measured function at a later point where
	 * the result of the operation is actually needed.
	 *
	 * \sa getTotalElementCreationTimings() and getTotalVectorNormTimings()
	 *
	 * @return total time spent in calls of operation
	 */
	VectorSpaceOperations::Timing_t getTotalLinearTimings() const
	{ return boost::fusion::at_key<VectorSpaceOperations::ScalarVectorMultiplication>(instance).second
			+boost::fusion::at_key<VectorSpaceOperations::VectorAddition>(instance).second
			+boost::fusion::at_key<VectorSpaceOperations::VectorAssignment>(instance).second
			+boost::fusion::at_key<VectorSpaceOperations::VectorComparison>(instance).second
			+boost::fusion::at_key<VectorSpaceOperations::VectorModification>(instance).second
			+boost::fusion::at_key<VectorSpaceOperations::VectorMultiplication>(instance).second
			+boost::fusion::at_key<VectorSpaceOperations::VectorNorm>(instance).second;
	}

	VectorSpaceOperations::Timing_t getTotalConstantTimings() const
	{ return boost::fusion::at_key<VectorSpaceOperations::ElementCreation>(instance).second; }

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
		TimeKeeper(VectorSpaceOperations::CountTiming_t &_optype) :
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
		VectorSpaceOperations::CountTiming_t &optype;
		//!> contains the starting time
		const boost::chrono::high_resolution_clock::time_point timing_start;
	};
};



#endif /* MINIMIZATIONS_SPACES_VECTORSPACEOPERATIONCOUNTS_HPP_ */
