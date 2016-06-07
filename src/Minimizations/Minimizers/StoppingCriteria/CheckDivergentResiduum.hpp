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
		oldresiduum(std::numeric_limits<double>::max())
	{}

	bool operator()(
			const boost::chrono::duration<double> &_time,
			const int _current_outeriterations,
			const double _residuum,
			const double _ynorm) const
	{
		if (oldresiduum >= _residuum) {
			oldresiduum = _residuum;
			return false;
		} else
			return true;
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
};


#endif /* MINIMIZATIONS_MINIMIZERS_STOPPINGCRITERIA_CHECKDIVERGENTRESIDUUM_HPP_ */
