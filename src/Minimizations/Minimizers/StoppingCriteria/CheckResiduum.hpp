/*
 * CheckResiduum.hpp
 *
 *  Created on: Oct 14, 2015
 *      Author: heber
 */

#ifndef MINIMIZATIONS_MINIMIZERS_STOPPINGCRITERIA_CHECKRESIDUUM_HPP_
#define MINIMIZATIONS_MINIMIZERS_STOPPINGCRITERIA_CHECKRESIDUUM_HPP_

#include "BassoConfig.h"


#include "Minimizations/Minimizers/StoppingCriteria/StoppingArguments.hpp"
#include "Minimizations/Minimizers/StoppingCriteria/StoppingCriterion_impl.hpp"

struct CheckResiduum : public StoppingCriterion_impl
{
	CheckResiduum(const StoppingArguments &_args) :
		args(_args)
	{}

	bool operator()(
			const boost::chrono::duration<double> &_time,
			const int _current_outeriterations,
			const double _residuum,
			const double _ynorm) const
	{
		return _residuum <= args.getTolerance()*args.getDiscrepancyParameter();
	}

	const StoppingArguments args;
};




#endif /* MINIMIZATIONS_MINIMIZERS_STOPPINGCRITERIA_CHECKRESIDUUM_HPP_ */
