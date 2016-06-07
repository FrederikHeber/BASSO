/*
 * StoppingCriterion_AND.hpp
 *
 *  Created on: Oct 14, 2015
 *      Author: heber
 */

#ifndef MINIMIZATIONS_MINIMIZERS_STOPPINGCRITERIA_STOPPINGCRITERION_NOT_HPP_
#define MINIMIZATIONS_MINIMIZERS_STOPPINGCRITERIA_STOPPINGCRITERION_NOT_HPP_

#include "BassoConfig.h"

#include <sstream>

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

	const std::string& getName() const
	{
		static const std::string name("NOT");
		return name;
	}

	/** Gives a statement on who currently says we should stop.
	 *
	 * TODO: whoIsTrue() is not fully working for StoppingCriterion_NOT but
	 * it also not easy to implement because we need to invert the operator
	 * in every instance of a possible impl hierarchy.
	 *
	 * @param _time
	 * @param _current_outeriterations
	 * @param _residuum
	 * @param _ynorm
	 * @return
	 */
	std::string whoIsTrue(const boost::chrono::duration<double> &_time,
			const int _current_outeriterations,
			const double _residuum,
			const double _ynorm) const
	{
		std::stringstream output;
		if (impl->operator()(_time, _current_outeriterations, _residuum, _ynorm))
			output << "NOT( ";
		else
			output << "( ";
		output << impl->getName();
		output << ")";
		return output.str();
	}

	const StoppingCriterion::ptr_t impl;
};



#endif /* MINIMIZATIONS_MINIMIZERS_STOPPINGCRITERIA_STOPPINGCRITERION_NOT_HPP_ */
