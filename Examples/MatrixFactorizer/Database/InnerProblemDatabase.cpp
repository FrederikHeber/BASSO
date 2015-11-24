/*
 * InnerProblemDatabase.cpp
 *
 *  Created on: Nov 23, 2015
 *      Author: heber
 */


#include "BassoConfig.h"

#include "InnerProblemDatabase.hpp"

#include "Log/Logging.hpp"

InnerProblemDatabase::InnerProblemDatabase()
{}

InnerProblemDatabase::~InnerProblemDatabase()
{
	BOOST_LOG_TRIVIAL(info)
			<< "The directory contains " << size() << "tables.";
	for (TableDirectory::TableIterator_t iter = directory.begin();
			iter != directory.end(); ++iter)
		BOOST_LOG_TRIVIAL(info)
				<< "Table " << iter->first << " has " << iter->second->size() << " tuples.";
}
