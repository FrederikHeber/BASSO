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
 * Database_TableDirectory.cpp
 *
 *  Created on: Nov 23, 2015
 *      Author: heber
 */


#include "BassoConfig.h"

#include <sstream>

#include "Log/Logging.hpp"

#include "TableDirectory.hpp"

Table& TableDirectory::addTable(const std::string &_name)
{
	// check whether table of such a name is not already present
	tables_t::iterator iter = tables.find(_name);
	if (iter == tables.end()) {
		Table::ptr insert_table(new Table(_name));
		tables.insert( std::make_pair(_name, insert_table) );
		iter = tables.find(_name);
	}
	return *iter->second;
}

bool TableDirectory::removeTable(const std::string &_name)
{
	tables_t::iterator iter = tables.find(_name);
	const bool status = iter != tables.end();
	if (status)
		tables.erase(iter);
	return status;
}

void TableDirectory::clear()
{
	for(tables_t::iterator iter = tables.begin();
		iter != tables.end(); ++iter) {
		LOG(debug, "Clearing table " << iter->second->getName());
		iter->second->clear();
	}
}


bool TableDirectory::clearTable(const std::string &_name)
{
	tables_t::iterator iter = tables.find(_name);
	const bool status = iter != tables.end();
	if (status)
		iter->second->clear();
	return status;
}

Table& TableDirectory::getTable(const std::string &_name)
{
	tables_t::iterator iter = tables.find(_name);
	assert( iter != tables.end() );
	return *iter->second;
}

const Table& TableDirectory::getTableConst(const std::string &_name) const
{
	tables_t::const_iterator iter = tables.find(_name);
	assert( iter != tables.end() );
	return *iter->second;
}

Table& TableDirectory::insertTable(
		const std::string &_name,
		const Table::ptr &_ref)
{
	tables.insert( std::make_pair(_name, _ref) );
	tables_t::iterator iter = tables.find(_name);
	assert( iter != tables.end() );
	return *iter->second;
}
