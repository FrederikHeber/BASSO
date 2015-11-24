/*
 * InnerProblemDatabase.cpp
 *
 *  Created on: Nov 23, 2015
 *      Author: heber
 */


#include "BassoConfig.h"

#include "InnerProblemDatabase.hpp"

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
{}

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
