### TESTS ###

########################## AuxiliaryConstraintsFactoryUnitTest #######################
add_executable (AuxiliaryConstraintsFactoryUnitTest 
	${PROJECT_SOURCE_DIR}/src/unittests/TestRunner.cpp 
	AuxiliaryConstraintsFactoryUnitTest.cpp) 
target_link_libraries (AuxiliaryConstraintsFactoryUnitTest 
	MatrixFactorizerBase
	${EXTRA_LIBS}
	${Boost_LIBRARIES}
	${CMAKE_THREAD_LIBS_INIT}
	"${CPPUNIT_LIBRARY}"
	${LIBRT_LIBRARIES}) 

add_test (AuxiliaryConstraintsFactoryUnitTest 
	${PROJECT_BINARY_DIR}/src/Solvers/AuxiliaryConstraints/unittests/AuxiliaryConstraintsFactoryUnitTest)

