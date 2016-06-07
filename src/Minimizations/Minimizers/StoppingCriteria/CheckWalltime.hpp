/*
 * CheckWalltime.hpp
 *
 *  Created on: Oct 14, 2015
 *      Author: heber
 */

#ifndef MINIMIZATIONS_MINIMIZERS_STOPPINGCRITERIA_CHECKWALLTIME_HPP_
#define MINIMIZATIONS_MINIMIZERS_STOPPINGCRITERIA_CHECKWALLTIME_HPP_

#include "BassoConfig.h"

#include <boost/chrono.hpp>

#include "Minimizations/Minimizers/StoppingCriteria/StoppingArguments.hpp"
#include "Minimizations/Minimizers/StoppingCriteria/StoppingCriterion_impl.hpp"

struct CheckWalltime : public StoppingCriterion_impl
{
	CheckWalltime(const StoppingArguments &_args) :
		args(_args)
	{}

	bool operator()(
			const boost::chrono::duration<double> &_time,
			const int _current_outeriterations,
			const double _residuum,
			const double _ynorm) const
	{
		return _time >= args.getMaxWalltime();
	}

	const std::string& getName() const
	{
		static const std::string name("WallTime");
		return name;
	}

	const StoppingArguments args;
};


#endif /* MINIMIZATIONS_MINIMIZERS_STOPPINGCRITERIA_CHECKWALLTIME_HPP_ */
