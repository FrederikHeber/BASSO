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
	${PROJECT_BINARY_DIR}/src/Minimizations/Functions/unittests/BregmanDistanceUnitTest)

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
	${PROJECT_BINARY_DIR}/src/Minimizations/Functions/unittests/BregmanProjectionFunctionalUnitTest)

########################## MetricProjectionFunctionalUnitTest #######################
add_executable (MetricProjectionFunctionalUnitTest 
	${PROJECT_SOURCE_DIR}/src/unittests/TestRunner.cpp 
	MetricProjectionFunctionalUnitTest.cpp) 
target_link_libraries (MetricProjectionFunctionalUnitTest 
	${EXTRA_LIBS}
	${Boost_LIBRARIES}
	"${CPPUNIT_LIBRARY}"
	${LIBRT_LIBRARIES})

add_test (MetricProjectionFunctionalUnitTest 
	${PROJECT_BINARY_DIR}/src/Minimizations/Functions/unittests/MetricProjectionFunctionalUnitTest)

########################## SmoothnessModulusUnitTest #######################
add_executable (SmoothnessModulusUnitTest 
	${PROJECT_SOURCE_DIR}/src/unittests/TestRunner.cpp 
	SmoothnessModulusUnitTest.cpp) 
target_link_libraries (SmoothnessModulusUnitTest 
	${EXTRA_LIBS}
	"${CPPUNIT_LIBRARY}")

add_test (SmoothnessModulusUnitTest 
	${PROJECT_BINARY_DIR}/src/Minimizations/Functions/unittests/SmoothnessModulusUnitTest)

########################## VectorProjectionUnitTest #######################
add_executable (VectorProjectionUnitTest 
	${PROJECT_SOURCE_DIR}/src/unittests/TestRunner.cpp 
	VectorProjectionUnitTest.cpp) 
target_link_libraries (VectorProjectionUnitTest 
	${EXTRA_LIBS}
	${GSL_LIBRARY} 
	${Boost_LIBRARIES}
	${CMAKE_THREAD_LIBS_INIT}
	"${CPPUNIT_LIBRARY}"
	${LIBRT_LIBRARIES}) 

add_test (VectorProjectionUnitTest 
	${PROJECT_BINARY_DIR}/src/Minimizations/Functions/unittests/VectorProjectionUnitTest)

# mark as to fail if required
#set_tests_properties (BregmanFunctionalUnitTest
#	PROPERTIES WILL_FAIL true)
