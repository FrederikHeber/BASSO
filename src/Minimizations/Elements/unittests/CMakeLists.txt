### TESTS ###

########################## LinearDependencyCheckerUnitTest #######################
add_executable (LinearDependencyCheckerUnitTest 
	${PROJECT_SOURCE_DIR}/src/unittests/TestRunner.cpp 
	LinearDependencyCheckerUnitTest.cpp) 
target_link_libraries (LinearDependencyCheckerUnitTest 
	${EXTRA_LIBS}
	${BOOST_CORE_LIBRARIES}
	${CMAKE_THREAD_LIBS_INIT}
	"${CPPUNIT_LIBRARY}"
	${LIBRT_LIBRARIES}) 

add_test (LinearDependencyCheckerUnitTest 
	${PROJECT_BINARY_DIR}/src/Minimizations/Elements/unittests/LinearDependencyCheckerUnitTest)

