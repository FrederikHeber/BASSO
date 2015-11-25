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
		Table &_table
		) const
{
	overall_table_accumulator.insertValues(
			_table);
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
