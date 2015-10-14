/*
 * StoppingCriterion_OR.hpp
 *
 *  Created on: Oct 14, 2015
 *      Author: heber
 */

#ifndef MINIMIZATIONS_MINIMIZERS_STOPPINGCRITERIA_STOPPINGCRITERION_OR_HPP_
#define MINIMIZATIONS_MINIMIZERS_STOPPINGCRITERIA_STOPPINGCRITERION_OR_HPP_

#include "BassoConfig.h"

#include "Minimizations/Minimizers/StoppingCriteria/StoppingCriterion.hpp"

struct StoppingCriterion_OR : public StoppingCriterion_impl
{
	StoppingCriterion_OR(
			const StoppingCriterion::ptr_t &_left,
			const StoppingCriterion::ptr_t &_right) :
				left(_left),
				right(_right)
	{}

	bool operator()(
			const boost::chrono::duration<double> &_time,
			const int _current_outeriterations,
			const double _residuum,
			const double _ynorm) const
	{
		return (left->operator()(_time, _current_outeriterations, _residuum, _ynorm)
			|| right->operator()(_time, _current_outeriterations, _residuum, _ynorm));
	}

	const StoppingCriterion::ptr_t left;
	const StoppingCriterion::ptr_t right;
};



#endif /* MINIMIZATIONS_MINIMIZERS_STOPPINGCRITERIA_STOPPINGCRITERION_OR_HPP_ */
