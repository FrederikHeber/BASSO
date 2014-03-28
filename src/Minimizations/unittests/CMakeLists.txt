### TESTS ###

########################## DualityMappingUnitTest #######################
add_executable (DualityMappingUnitTest 
	${PROJECT_SOURCE_DIR}/src/unittests/TestRunner.cpp 
	DualityMappingUnitTest.cpp) 
target_link_libraries (DualityMappingUnitTest 
	${EXTRA_LIBS}
	"${CPPUNIT_LIBRARY}")

add_test (DualityMappingUnitTest 
	${PROJECT_BINARY_DIR}/src/Minimizations/unittests/DualityMappingUnitTest)

########################## BregmanFunctionalUnitTest #######################
add_executable (BregmanFunctionalUnitTest 
	${PROJECT_SOURCE_DIR}/src/unittests/TestRunner.cpp 
	BregmanFunctionalUnitTest.cpp) 
target_link_libraries (BregmanFunctionalUnitTest 
	${EXTRA_LIBS}
	"${CPPUNIT_LIBRARY}")

# does the application run
add_test (BregmanFunctionalUnitTest 
	${PROJECT_BINARY_DIR}/src/Minimizations/unittests/BregmanFunctionalUnitTest)

# mark as to fail if required
#set_tests_properties (BregmanFunctionalUnitTest
#	PROPERTIES WILL_FAIL true)
