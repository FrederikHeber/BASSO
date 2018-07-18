/*
 * StoppingCriterion_AND.hpp
 *
 *  Created on: Oct 14, 2015
 *      Author: heber
 */

#ifndef MINIMIZATIONS_MINIMIZERS_STOPPINGCRITERIA_STOPPINGCRITERION_AND_HPP_
#define MINIMIZATIONS_MINIMIZERS_STOPPINGCRITERIA_STOPPINGCRITERION_AND_HPP_

#include "BassoConfig.h"

#include <sstream>

#include "Minimizations/Minimizers/StoppingCriteria/StoppingCriterion.hpp"

struct StoppingCriterion_AND : public StoppingCriterion_impl
{
	StoppingCriterion_AND(
			const StoppingCriterion::ptr_t &_left,
			const StoppingCriterion::ptr_t &_right) :
				left(_left),
				right(_right)
	{}

	virtual ~StoppingCriterion_AND() {}

	bool operator()(
			const boost::chrono::duration<double> &_time,
			const int _current_outeriterations,
			const double _residuum,
			const double _ynorm) const
	{
		return (left->operator()(_time, _current_outeriterations, _residuum, _ynorm)
			&& right->operator()(_time, _current_outeriterations, _residuum, _ynorm));
	}
	const std::string& getName() const
	{
		static const std::string name("AND");
		return name;
	}

	std::string whoIsTrue(const boost::chrono::duration<double> &_time,
			const int _current_outeriterations,
			const double _residuum,
			const double _ynorm) const
	{
		const bool left_true = left->operator()(_time, _current_outeriterations, _residuum, _ynorm);
		const bool right_true = right->operator()(_time, _current_outeriterations, _residuum, _ynorm);
		if (left_true && right_true) {
			std::stringstream output;
			output << left->whoIsTrue(_time, _current_outeriterations, _residuum, _ynorm);
			output << " && ";
			output << right->whoIsTrue(_time, _current_outeriterations, _residuum, _ynorm);
			return output.str();
		} else {
			return std::string();
		}
	}

	const StoppingCriterion::ptr_t left;
	const StoppingCriterion::ptr_t right;
};



#endif /* MINIMIZATIONS_MINIMIZERS_STOPPINGCRITERIA_STOPPINGCRITERION_AND_HPP_ */
