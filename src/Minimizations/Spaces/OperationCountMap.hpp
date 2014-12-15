/*
 * OperationCountMap.hpp
 *
 *  Created on: Dec 12, 2014
 *      Author: heber
 */

#ifndef MINIMIZATIONS_SPACES_OPERATIONCOUNTMAP_HPP_
#define MINIMIZATIONS_SPACES_OPERATIONCOUNTMAP_HPP_

#include "BassoConfig.h"

#include <boost/fusion/container.hpp>
#include <boost/mpl/vector.hpp>

namespace VectorSpaceOperations {
	//!> typedef for the count type
	typedef int Count_t;
	//!> typedef for the time-keeping type
	typedef boost::chrono::nanoseconds Timing_t;
	//!> typedef for the combined count and time-keeping type
	typedef std::pair<Count_t, Timing_t> CountTiming_t;

	struct ElementCreation {};
	struct ScalarVectorMultiplication {};
	struct VectorAddition {};
	struct VectorAssignment {};
	struct VectorComparison {};
	struct VectorModification {};
	struct VectorMultiplication {};
	struct VectorNorm {};
};

/** A boost::fusion::container allows to discern between multiple instances
 * of the same type.
 *
 */
typedef boost::fusion::map<
		boost::fusion::pair<VectorSpaceOperations::ElementCreation, VectorSpaceOperations::CountTiming_t>,
		boost::fusion::pair<VectorSpaceOperations::ScalarVectorMultiplication, VectorSpaceOperations::CountTiming_t>,
		boost::fusion::pair<VectorSpaceOperations::VectorAddition, VectorSpaceOperations::CountTiming_t>,
		boost::fusion::pair<VectorSpaceOperations::VectorAssignment, VectorSpaceOperations::CountTiming_t>,
		boost::fusion::pair<VectorSpaceOperations::VectorComparison, VectorSpaceOperations::CountTiming_t>,
		boost::fusion::pair<VectorSpaceOperations::VectorModification, VectorSpaceOperations::CountTiming_t>,
		boost::fusion::pair<VectorSpaceOperations::VectorMultiplication, VectorSpaceOperations::CountTiming_t>,
		boost::fusion::pair<VectorSpaceOperations::VectorNorm, VectorSpaceOperations::CountTiming_t>
> OperationCountMap_t;

/** In this boost::mpl::list we enumerate all posssible operation counts.
 *
 * \todo use boost::fusion::copy
 */
typedef boost::mpl::vector<
		VectorSpaceOperations::ElementCreation,
		VectorSpaceOperations::ScalarVectorMultiplication,
		VectorSpaceOperations::VectorAddition,
		VectorSpaceOperations::VectorAssignment,
		VectorSpaceOperations::VectorComparison,
		VectorSpaceOperations::VectorModification,
		VectorSpaceOperations::VectorMultiplication,
		VectorSpaceOperations::VectorNorm
> OperationCountVector_t;

/** In this boost::mpl::list we enumerate all posssible operation counts
 * that have constant scaling costs
 *
 */
typedef boost::mpl::vector<
		VectorSpaceOperations::ElementCreation
> ConstantOperationCountVector_t;

/** In this boost::mpl::list we enumerate all posssible operation counts
 * that have linear scaling costs
 *
 */
typedef boost::mpl::vector<
		VectorSpaceOperations::ScalarVectorMultiplication,
		VectorSpaceOperations::VectorAddition,
		VectorSpaceOperations::VectorAssignment,
		VectorSpaceOperations::VectorComparison,
		VectorSpaceOperations::VectorModification,
		VectorSpaceOperations::VectorMultiplication,
		VectorSpaceOperations::VectorNorm
> LinearOperationCountVector_t;


#endif /* MINIMIZATIONS_SPACES_OPERATIONCOUNTMAP_HPP_ */
