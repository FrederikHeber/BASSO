/*
 * StoppingCriterionStoppingCriterion.hpp.hpp
 *
 *  Created on: Oct 14, 2015
 *      Author: heber
 */

#ifndef MINIMIZATIONS_MINIMIZERS_STOPPINGCRITERIA_STOPPINGCRITERION_IMPL_HPP_
#define MINIMIZATIONS_MINIMIZERS_STOPPINGCRITERIA_STOPPINGCRITERION_IMPL_HPP_

#include "BassoConfig.h"

#include <boost/chrono.hpp>

/** Interface definition for a stopping criterion.
 *
 * This defines only the function itself. In StoppingCriterion we
 * wrap these "implementations" inside a hidden shared_ptr to allow
 * for combining these predicates.
 *
 */
struct StoppingCriterion_impl
{
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
	virtual bool operator()(
			const boost::chrono::duration<double> &_time,
			const int _current_outeriterations,
			const double _residuum,
			const double _ynorm) const = 0;

	/** Returns the specific name of this stopping condition."
	 *
	 * @return name of condition
	 */
	virtual const std::string & getName() const = 0;

	/** Gives a statement on who currently says we should stop.
	 *
	 * @param _time current spent time on iteratging
	 * @param _current_outeriterations current number of iteration steps
	 * @param _residuum residuum at current iteration
	 * @param _ynorm norm of right-hand side (y)
	 * @return name of stopping criterion if true, else - empty string
	 */
	virtual std::string whoIsTrue(
			const boost::chrono::duration<double> &_time,
			const int _current_outeriterations,
			const double _residuum,
			const double _ynorm
			) const
	{
		if (operator()(_time, _current_outeriterations, _residuum, _ynorm))
			return getName();
		else
			return std::string();
	}
};

#endif /* MINIMIZATIONS_MINIMIZERS_STOPPINGCRITERIA_STOPPINGCRITERION_IMPL_HPP_ */
