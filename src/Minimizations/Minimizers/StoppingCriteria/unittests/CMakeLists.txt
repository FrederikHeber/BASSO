### TESTS ###

########################## StoppingCriteriaFactoryUnitTest #######################
add_executable (StoppingCriteriaFactoryUnitTest 
	${PROJECT_SOURCE_DIR}/src/unittests/TestRunner.cpp 
	StoppingCriteriaFactoryUnitTest.cpp) 
target_link_libraries (StoppingCriteriaFactoryUnitTest 
	${EXTRA_LIBS}
	${BOOST_CORE_LIBRARIES}
	${CMAKE_THREAD_LIBS_INIT}
	"${CPPUNIT_LIBRARY}"
	${LIBRT_LIBRARIES}) 

add_test (StoppingCriteriaFactoryUnitTest 
	${PROJECT_BINARY_DIR}/src/Minimizations/Minimizers/StoppingCriteria/unittests/StoppingCriteriaFactoryUnitTest)

