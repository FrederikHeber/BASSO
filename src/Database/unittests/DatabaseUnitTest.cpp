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
