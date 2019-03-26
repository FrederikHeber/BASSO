/*
 * Project: BASSO - BAnach Sequential Subspace Optimizer
 * Description: C++ library for solving optimization problems in Banach spaces
 * Copyright (C)  2014-2016 University of Saarland. All rights reserved.
 * Copyright (C)  2016-2019 Frederik Heber. All rights reserved.
 *
 *
 *   This file is part of the BASSO library.
 *
 *    BASSO is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 2 of the License, or
 *    (at your option) any later version.
 *
 *    BASSO is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with BASSO.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

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

const std::string& StoppingCriterion::getName() const
{
	return internal_impl->getName();
}

std::string StoppingCriterion::whoIsTrue(
		const boost::chrono::duration<double> &_time,
		const int _current_outeriterations,
		const double _residuum,
		const double _ynorm) const
{
	return internal_impl->whoIsTrue(_time, _current_outeriterations, _residuum, _ynorm);
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


