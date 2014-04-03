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

add_test (BregmanFunctionalUnitTest 
	${PROJECT_BINARY_DIR}/src/Minimizations/unittests/BregmanFunctionalUnitTest)

########################## LpNormUnitTest #######################
add_executable (LpNormUnitTest 
	${PROJECT_SOURCE_DIR}/src/unittests/TestRunner.cpp 
	LpNormUnitTest.cpp) 
target_link_libraries (LpNormUnitTest 
	${EXTRA_LIBS}
	"${CPPUNIT_LIBRARY}")

add_test (LpNormUnitTest 
	${PROJECT_BINARY_DIR}/src/Minimizations/unittests/LpNormUnitTest)

########################## SmoothnessModulusUnitTest #######################
add_executable (SmoothnessModulusUnitTest 
	${PROJECT_SOURCE_DIR}/src/unittests/TestRunner.cpp 
	SmoothnessModulusUnitTest.cpp) 
target_link_libraries (SmoothnessModulusUnitTest 
	${EXTRA_LIBS}
	"${CPPUNIT_LIBRARY}")

add_test (SmoothnessModulusUnitTest 
	${PROJECT_BINARY_DIR}/src/Minimizations/unittests/SmoothnessModulusUnitTest)


# mark as to fail if required
#set_tests_properties (BregmanFunctionalUnitTest
#	PROPERTIES WILL_FAIL true)
