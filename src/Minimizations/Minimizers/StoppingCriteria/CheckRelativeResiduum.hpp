/*
 * CheckRelativeResiduum.hpp
 *
 *  Created on: Oct 14, 2015
 *      Author: heber
 */

#ifndef MINIMIZATIONS_MINIMIZERS_STOPPINGCRITERIA_CHECKRELATIVERESIDUUM_HPP_
#define MINIMIZATIONS_MINIMIZERS_STOPPINGCRITERIA_CHECKRELATIVERESIDUUM_HPP_

#include "BassoConfig.h"


#include "Minimizations/Minimizers/StoppingCriteria/StoppingArguments.hpp"
#include "Minimizations/Minimizers/StoppingCriteria/StoppingCriterion_impl.hpp"

struct CheckRelativeResiduum : public StoppingCriterion_impl
{
	CheckRelativeResiduum(const StoppingArguments &_args) :
		args(_args)
	{}

	virtual ~CheckRelativeResiduum() {}

	bool operator()(
			const boost::chrono::duration<double> &_time,
			const int _current_outeriterations,
			const double _residuum,
			const double _ynorm) const
	{
		if (fabs(_ynorm) > BASSOTOLERANCE)
			return _residuum/_ynorm <= args.getTolerance()*args.getDiscrepancyParameter();
		else
			return false;
	}

	const std::string& getName() const
	{
		static const std::string name("RelativeResiduum");
		return name;
	}

	const StoppingArguments args;
};


#endif /* MINIMIZATIONS_MINIMIZERS_STOPPINGCRITERIA_CHECKRELATIVERESIDUUM_HPP_ */
