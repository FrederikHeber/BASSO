/*
 * StoppingCriterion_AND.hpp
 *
 *  Created on: Oct 14, 2015
 *      Author: heber
 */

#ifndef MINIMIZATIONS_MINIMIZERS_STOPPINGCRITERIA_STOPPINGCRITERION_NOT_HPP_
#define MINIMIZATIONS_MINIMIZERS_STOPPINGCRITERIA_STOPPINGCRITERION_NOT_HPP_

#include "BassoConfig.h"

#include "Minimizations/Minimizers/StoppingCriteria/StoppingCriterion.hpp"

struct StoppingCriterion_NOT : public StoppingCriterion_impl
{
	StoppingCriterion_NOT(
			const StoppingCriterion::ptr_t &_impl) :
				impl(_impl)
	{}

	bool operator()(
			const boost::chrono::duration<double> &_time,
			const int _current_outeriterations,
			const double _residuum,
			const double _ynorm) const
	{
		return !impl->operator()(_time, _current_outeriterations, _residuum, _ynorm);
	}

	const StoppingCriterion::ptr_t impl;
};



#endif /* MINIMIZATIONS_MINIMIZERS_STOPPINGCRITERIA_STOPPINGCRITERION_NOT_HPP_ */
