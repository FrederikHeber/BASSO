/*
 * CheckRelativeChangeResiduum.hpp
 *
 *  Created on: Nov 05, 2015
 *      Author: heber
 */

#ifndef MINIMIZATIONS_MINIMIZERS_STOPPINGCRITERIA_CHECKRELATIVECHANGERESIDUUM_HPP_
#define MINIMIZATIONS_MINIMIZERS_STOPPINGCRITERIA_CHECKRELATIVECHANGERESIDUUM_HPP_

#include "BassoConfig.h"


#include "Minimizations/Minimizers/StoppingCriteria/StoppingArguments.hpp"
#include "Minimizations/Minimizers/StoppingCriteria/StoppingCriterion_impl.hpp"

struct CheckRelativeChangeResiduum : public StoppingCriterion_impl
{
	CheckRelativeChangeResiduum(const StoppingArguments &_args) :
		args(_args),
		oldresidual(0.)
	{}

	virtual ~CheckRelativeChangeResiduum() {}

	bool operator()(
			const boost::chrono::duration<double> &,
			const int,
			const double _residuum,
			const double) const
	{
		bool result = true;
		if (fabs(_residuum) > BASSOTOLERANCE)
			result = (fabs(oldresidual -_residuum)/_residuum)
					<= args.getTolerance()*args.getDiscrepancyParameter();
		oldresidual = _residuum;
		return result;
	}

	const std::string& getName() const
	{
		static const std::string name("RelativeChangeResiduum");
		return name;
	}

	const StoppingArguments args;
private:

	//!> contains residual from last iteration
	mutable double oldresidual;
};


#endif /* MINIMIZATIONS_MINIMIZERS_STOPPINGCRITERIA_CHECKRELATIVECHANGERESIDUUM_HPP_ */
