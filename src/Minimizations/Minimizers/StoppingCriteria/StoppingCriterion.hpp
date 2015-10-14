/*
 * StoppingCriterion.hpp
 *
 *  Created on: Oct 14, 2015
 *      Author: heber
 */

#ifndef MINIMIZATIONS_MINIMIZERS_STOPPINGCRITERION_STOPPINGCRITERIA_HPP_
#define MINIMIZATIONS_MINIMIZERS_STOPPINGCRITERION_STOPPINGCRITERIA_HPP_

#include "BassoConfig.h"

#include <boost/chrono.hpp>
#include <boost/shared_ptr.hpp>

#include "Minimizations/Minimizers/StoppingCriteria/StoppingCriterion_impl.hpp"

/** Interface definition for a stopping criterion.
 *
 */
struct StoppingCriterion
{
	//!> typedef for shared ptr containing StoppingCriterion
	typedef boost::shared_ptr<StoppingCriterion_impl> ptr_t;

	/** Functor to check the respective stopping criterion.
	 *
	 * We have a broad interface, requesting all possible values
	 * where each of the stopping criteria needs only a few.
	 *
	 * @param _time current spent time on iteratging
	 * @param _current_outeriterations current number of iteration steps
	 * @param _residuum residuum at current iteration
	 * @param _ynorm norm of right-hand side (y)
	 * @return true - stop, false - continue
	 */
	bool operator()(
			const boost::chrono::duration<double> &_time,
			const int _current_outeriterations,
			const double _residuum,
			const double _ynorm) const;

private:
	//!> internal predicate evaluating the stopping criterion or a combination
	const ptr_t internal_impl;
};

// boolean operators to combine stopping criteria
StoppingCriterion::ptr_t operator&&(
		const StoppingCriterion::ptr_t &_a,
		const StoppingCriterion::ptr_t &_b);

StoppingCriterion::ptr_t operator||(
		const StoppingCriterion::ptr_t &_a,
		const StoppingCriterion::ptr_t &_b);

StoppingCriterion::ptr_t operator!(
		const StoppingCriterion::ptr_t &_a);

#endif /* MINIMIZATIONS_MINIMIZERS_STOPPINGCRITERION_STOPPINGCRITERIA_HPP_ */
