/*
 * CheckDivergentResiduum.hpp
 *
 *  Created on: Oct 14, 2015
 *      Author: heber
 */

#ifndef MINIMIZATIONS_MINIMIZERS_STOPPINGCRITERIA_CHECKDIVERGENTRESIDUUM_HPP_
#define MINIMIZATIONS_MINIMIZERS_STOPPINGCRITERIA_CHECKDIVERGENTRESIDUUM_HPP_

#include "BassoConfig.h"

#include <limits>

#include "Minimizations/Minimizers/StoppingCriteria/StoppingArguments.hpp"
#include "Minimizations/Minimizers/StoppingCriteria/StoppingCriterion_impl.hpp"

struct CheckDivergentResiduum : public StoppingCriterion_impl
{
	CheckDivergentResiduum(const StoppingArguments &_args) :
		args(_args),
		oldresiduum(std::numeric_limits<double>::max()),
		violation_counts(0),
		MAX_COUNTS(3)
	{}

	virtual ~CheckDivergentResiduum() {}

	bool operator()(
			const boost::chrono::duration<double> &_time,
			const int _current_outeriterations,
			const double _residuum,
			const double _ynorm) const
	{
		if (oldresiduum >= _residuum) {
			oldresiduum = _residuum;
			return false;
		} else {
			oldresiduum = _residuum;
			++violation_counts;
			if (violation_counts > MAX_COUNTS)
				return true;
			else
				return false;
		}
	}

	const std::string& getName() const
	{
		static const std::string name("DivergentResiduum");
		return name;
	}

	const StoppingArguments args;

private:
	//!> stores old residuum from last iteration
	mutable double oldresiduum;
	//!> counts how often this check was violated
	mutable unsigned int violation_counts;
	//!> how often check may be violated before we trigger stop
	const unsigned int MAX_COUNTS;
};


#endif /* MINIMIZATIONS_MINIMIZERS_STOPPINGCRITERIA_CHECKDIVERGENTRESIDUUM_HPP_ */
