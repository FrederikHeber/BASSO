/*
 * RelativeShrinkageMapping.cpp
 *
 *  Created on: Oct 13, 2014
 *      Author: heber
 */

#include "RelativeShrinkageMapping.hpp"

#include "Log/Logging.hpp"
#include "Math/Helpers.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Mappings/IllegalDualityMapping.hpp"
#include "Minimizations/Minimizers/MinimizationExceptions.hpp"
#include "Minimizations/Norms/L1Norm.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"

RelativeShrinkageMapping::RelativeShrinkageMapping(
		const NormedSpace_weakptr_t &_NormedSpaceRef) :
	L1DualityMapping(_NormedSpaceRef, 2.),
	lambda(0.1),
	count(0),
	timing(boost::chrono::nanoseconds(0))
{}

RelativeShrinkageMapping::RelativeShrinkageMapping(
		const NormedSpace_weakptr_t &_NormedSpaceRef,
		const double _lambda) :
	L1DualityMapping(_NormedSpaceRef, 2.),
	lambda(_lambda),
	count(0),
	timing(boost::chrono::nanoseconds(0))
{
	if (!(lambda > 0))
		throw MinimizationIllegalValue_exception()
		<< MinimizationIllegalValue_name("lambda");
}

const SpaceElement_ptr_t
RelativeShrinkageMapping::operator()(
		const SpaceElement_ptr_t &_x) const
{
	// start timing
	const boost::chrono::high_resolution_clock::time_point timing_start =
			boost::chrono::high_resolution_clock::now();

	SpaceElement_ptr_t result = getTargetSpace()->createElement();
	for (unsigned int i=0;i<result->getSpace()->getDimension();++i)
		(*result)[i] = fabs((*_x)[i]) < lambda ?
				0. :
				(fabs((*_x)[i])-lambda)*Helpers::sign((*_x)[i]);
	*result *= 1./lambda;

	// finish timing
	const boost::chrono::high_resolution_clock::time_point timing_end =
			boost::chrono::high_resolution_clock::now();
	timing += timing_end - timing_start;
	++count;

	return result;
}

const double
RelativeShrinkageMapping::getRelativeShrinkage(
		const SpaceElement_ptr_t &_x
		) const
{
	/// we use algorithm 1 given in [Schoepfer 2012]
	if (_x->isZero(BASSOTOLERANCE))
		return 0.;

	SpaceElement_ptr_t tempAbsvector = _x->getSpace()->createElement();
	*tempAbsvector = _x->getAbsVector();
	size_t active_set_size = 0;
	double l1norm_sum = 0.;
	double coefficient = 0.;
	std::vector<int> active_set;

	// loop over element's components and check against current coefficient
	do {
		// extract max component
		const std::pair<double, int> maxCoeff =
				tempAbsvector->getMaxCoefficientAndIndex();
		// compare with coefficient
//		BOOST_LOG_TRIVIAL(info)
//			<< "Is current max " << maxCoeff.first
//			<< " greater than coefficient " << coefficient << "?";
		if (coefficient > maxCoeff.first + BASSOTOLERANCE)
			break;
		// add to l1norm sum and to active set
//		BOOST_LOG_TRIVIAL(info)
//			<< "Adding " << maxCoeff.second << " to active_set.";
		l1norm_sum += maxCoeff.first;
		++active_set_size;
		active_set.push_back(maxCoeff.second);
		// remove this component from max component search by zeroing
		(*tempAbsvector)[maxCoeff.second] = 0.;
		// and recalculate coefficient
		coefficient = l1norm_sum/((double)active_set_size + lambda);
	} while (active_set_size < tempAbsvector->getSpace()->getDimension());


	// sensibility checks
	// all component in set exceed coefficient
	std::sort(active_set.begin(), active_set.end());
//	{
//		std::stringstream output;
//		std::copy(
//				active_set.begin(), active_set.end(),
//				std::ostream_iterator<int>(output, ";"));
//		BOOST_LOG_TRIVIAL(info)
//				<< "Final active_set is " << output.str();
//	}
	for (std::vector<int>::const_iterator iter = active_set.begin();
			iter != active_set.end(); ++iter) {
//		BOOST_LOG_TRIVIAL(info)
//				<< "Is |_x[" << *iter << "]| = "
//				<<  fabs((*_x)[*iter]) << " >= " << coefficient << "?";
		assert( fabs((*_x)[*iter]) + BASSOTOLERANCE >= coefficient );
	}
	// all other components don't
	std::vector<int> complement_set(_x->getSpace()->getDimension(), 0);
	std::generate(complement_set.begin(), complement_set.end(), Helpers::unique_number());
	std::vector<int>::iterator eraseiter =
			std::set_difference(
				complement_set.begin(), complement_set.end(),
				active_set.begin(), active_set.end(),
				complement_set.begin());
	for (std::vector<int>::const_iterator iter = complement_set.begin();
			iter != eraseiter; ++iter) {
//		BOOST_LOG_TRIVIAL(info)
//				<< "Is |_x[" << *iter << "]| = "
//				<<  fabs((*_x)[*iter]) << " < " << coefficient << "?";
		assert( fabs((*_x)[*iter]) < coefficient + BASSOTOLERANCE);
	}

	return coefficient;
}

const Mapping_ptr_t RelativeShrinkageMapping::getAdjointMapping() const
{
	Mapping_ptr_t instance(
			new IllegalDualityMapping
	);
	return instance;
}

