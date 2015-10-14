/*
 * StoppingCriterion.cpp
 *
 *  Created on: Oct 14, 2015
 *      Author: heber
 */


#include "BassoConfig.h"

#include "StoppingCriterion.hpp"

#include "Minimizations/Minimizers/StoppingCriteria/StoppingCriterion_AND.hpp"
#include "Minimizations/Minimizers/StoppingCriteria/StoppingCriterion_NOT.hpp"
#include "Minimizations/Minimizers/StoppingCriteria/StoppingCriterion_OR.hpp"

bool StoppingCriterion::operator()(
		const boost::chrono::duration<double> &_time,
		const int _current_outeriterations,
		const double _residuum,
		const double _ynorm) const
{
	return internal_impl->operator()(_time, _current_outeriterations, _residuum, _ynorm);
}

StoppingCriterion::ptr_t operator&&(
		const StoppingCriterion::ptr_t &_a,
		const StoppingCriterion::ptr_t &_b)
{
	StoppingCriterion::ptr_t criterion(new StoppingCriterion_AND(_a, _b));
	return criterion;
}

StoppingCriterion::ptr_t operator||(
		const StoppingCriterion::ptr_t &_a,
		const StoppingCriterion::ptr_t &_b)
{
	StoppingCriterion::ptr_t criterion(new StoppingCriterion_OR(_a, _b));
	return criterion;
}


StoppingCriterion::ptr_t operator!(
		const StoppingCriterion::ptr_t &_a)
{
	StoppingCriterion::ptr_t criterion(new StoppingCriterion_NOT(_a));
	return criterion;
}


