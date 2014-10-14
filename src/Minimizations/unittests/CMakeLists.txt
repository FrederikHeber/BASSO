### TESTS ###

########################## BregmanDistanceUnitTest #######################
add_executable (BregmanDistanceUnitTest 
	${PROJECT_SOURCE_DIR}/src/unittests/TestRunner.cpp 
	BregmanDistanceUnitTest.cpp) 
target_link_libraries (BregmanDistanceUnitTest 
	${EXTRA_LIBS}
	${Boost_LIBRARIES}
	${CMAKE_THREAD_LIBS_INIT}
	"${CPPUNIT_LIBRARY}"
	${LIBRT_LIBRARIES}) 

add_test (BregmanDistanceUnitTest 
	${PROJECT_BINARY_DIR}/src/Minimizations/unittests/BregmanDistanceUnitTest)

########################## BregmanProjectionFunctionalUnitTest #######################
add_executable (BregmanProjectionFunctionalUnitTest 
	${PROJECT_SOURCE_DIR}/src/unittests/TestRunner.cpp 
	BregmanProjectionFunctionalUnitTest.cpp) 
target_link_libraries (BregmanProjectionFunctionalUnitTest 
	${EXTRA_LIBS}
	${Boost_LIBRARIES}
	"${CPPUNIT_LIBRARY}"
	${LIBRT_LIBRARIES})

add_test (BregmanProjectionFunctionalUnitTest 
	${PROJECT_BINARY_DIR}/src/Minimizations/unittests/BregmanProjectionFunctionalUnitTest)

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
