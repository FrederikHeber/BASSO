/*
 * Project: BASSO - BAnach Sequential Subspace Optimizer
 * Description: C++ library for solving optimization problems in Banach spaces
 * Copyright (C)  2014-2016 University of Saarland. All rights reserved.
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
 * DatabaseUnitTest.cpp
 *
 *  Created on: Jun 12, 2014
 *      Author: heber
 */

#include "DatabaseUnitTest.hpp"

#include "Database/TableDirectoryDatabase_mock.hpp"
#include "Database/Table.hpp"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( DatabaseUnitTest );


void DatabaseUnitTest::setUp()
{
}


void DatabaseUnitTest::tearDown()
{
}

void DatabaseUnitTest::simpleTest()
{
    TableDirectoryDatabase_mock mydata;
    Table &mytable = mydata.addTable("mytable");
    Table::Tuple_t& tuple = mytable.getTuple();
    tuple.insert( std::make_pair("p", 1), Table::Parameter );
    mytable.addTuple(tuple);
    CPPUNIT_ASSERT_EQUAL( (size_t)1, mydata.size() );
}
