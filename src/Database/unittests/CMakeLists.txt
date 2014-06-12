### TESTS ###

########################## DualityMappingUnitTest #######################
add_executable (DatabaseUnitTest 
	${PROJECT_SOURCE_DIR}/src/unittests/TestRunner.cpp 
	DatabaseUnitTest.cpp) 
target_link_libraries (DatabaseUnitTest 
	${EXTRA_LIBS}
	${POCO_SQL_LIBRARIES}
	"${CPPUNIT_LIBRARY}")

add_test (DatabaseUnitTest 
	${PROJECT_BINARY_DIR}/src/Database/unittests/DatabaseUnitTest)

