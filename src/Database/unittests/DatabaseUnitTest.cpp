/*
 * DatabaseUnitTest.cpp
 *
 *  Created on: Jun 12, 2014
 *      Author: heber
 */

#include "DatabaseUnitTest.hpp"

#include "Database/Database.hpp"


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
    Database mydata;
    Database::Tuple_t tuple;
    tuple.insert( std::make_pair("p", 1) );
    mydata.addTuple(tuple);
    CPPUNIT_ASSERT_EQUAL( (size_t)1, mydata.size() );
}
