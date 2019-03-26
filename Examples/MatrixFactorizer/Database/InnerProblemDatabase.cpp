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
 * InnerProblemDatabase.cpp
 *
 *  Created on: Nov 23, 2015
 *      Author: heber
 */


#include "BassoConfig.h"

#include "InnerProblemDatabase.hpp"

#include <boost/assign.hpp>

#include "Database/Table_mock.hpp"

using namespace boost::assign;

InnerProblemDatabase::InnerProblemDatabase() :
	overall_table_accumulator(*this)
{}

InnerProblemDatabase::InnerProblemDatabase(
		const keys_t &_overall_keys
		) :
	overall_table_accumulator(
			*this,
			"data_overall",
			_overall_keys.begin(),
			_overall_keys.end())
{
	truetables += "data_overall";
}

Table& InnerProblemDatabase::addTable(const std::string &_name)
{
	if (truetables.count(_name) != 0) {
		// generate full table
		return TableDirectoryDatabase::addTable(_name);
	} else {
		// generate mockup
		Table::ptr insert_table(new Table_mock(_name));
		return directory.insertTable( _name, insert_table );
	}
}

void InnerProblemDatabase::clear()
{
	TableDirectoryDatabase::clear();
}

void InnerProblemDatabase::finish()
{
	// extract present data
	overall_table_accumulator.extractData();
}

void InnerProblemDatabase::insertAccumulatedValues(
		Table &_table,
		const std::string &_suffix
		) const
{
	overall_table_accumulator.insertValues(
			_table, _suffix);
}

void InnerProblemDatabase::insertValues(
		const std::vector<AccumulatedValues>::const_iterator &_values_first,
		const std::vector<AccumulatedValues>::const_iterator &_values_last
		)
{
	for (std::vector<AccumulatedValues>::const_iterator iter = _values_first;
			iter != _values_last; ++iter)
		overall_table_accumulator.addAccumulatedValues(*iter);
}
