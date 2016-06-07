/*
 * CheckIterationCount.hpp
 *
 *  Created on: Oct 14, 2015
 *      Author: heber
 */

#ifndef MINIMIZATIONS_MINIMIZERS_STOPPINGCRITERIA_CHECKITERATIONCOUNT_HPP_
#define MINIMIZATIONS_MINIMIZERS_STOPPINGCRITERIA_CHECKITERATIONCOUNT_HPP_

#include "BassoConfig.h"


#include "Minimizations/Minimizers/StoppingCriteria/StoppingArguments.hpp"
#include "Minimizations/Minimizers/StoppingCriteria/StoppingCriterion_impl.hpp"

struct CheckIterationCount : public StoppingCriterion_impl
{
	CheckIterationCount(const StoppingArguments &_args) :
		args(_args)
	{}

	bool operator()(
			const boost::chrono::duration<double> &_time,
			const int _current_outeriterations,
			const double _residuum,
			const double _ynorm) const
	{
		return _current_outeriterations >= args.getMaxIterations();
	}

	const std::string& getName() const
	{
		static const std::string name("IterationCount");
		return name;
	}

	const StoppingArguments args;
};




#endif /* MINIMIZATIONS_MINIMIZERS_STOPPINGCRITERIA_CHECKITERATIONCOUNT_HPP_ */
