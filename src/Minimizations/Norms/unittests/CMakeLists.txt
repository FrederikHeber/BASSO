### TESTS ###

########################## LpNormUnitTest #######################
add_executable (LpNormUnitTest 
	${PROJECT_SOURCE_DIR}/src/unittests/TestRunner.cpp 
	LpNormUnitTest.cpp) 
target_link_libraries (LpNormUnitTest 
	${EXTRA_LIBS}
	${Boost_LIBRARIES}
	"${CPPUNIT_LIBRARY}") 

add_test (LpNormUnitTest 
	${PROJECT_BINARY_DIR}/src/Minimizations/Norms/unittests/LpNormUnitTest)

########################## L1NormUnitTest #######################
add_executable (L1NormUnitTest
	${PROJECT_SOURCE_DIR}/src/unittests/TestRunner.cpp 
	L1NormUnitTest.cpp) 
target_link_libraries (L1NormUnitTest 
	${EXTRA_LIBS}
	${Boost_LIBRARIES}
	"${CPPUNIT_LIBRARY}") 

add_test (L1NormUnitTest 
	${PROJECT_BINARY_DIR}/src/Minimizations/Norms/unittests/L1NormUnitTest)

########################## LInfinityNormUnitTest #######################
add_executable (LInfinityNormUnitTest 
	${PROJECT_SOURCE_DIR}/src/unittests/TestRunner.cpp 
	LInfinityNormUnitTest.cpp) 
target_link_libraries (LInfinityNormUnitTest 
	${EXTRA_LIBS}
	${Boost_LIBRARIES}
	"${CPPUNIT_LIBRARY}") 

add_test (LInfinityNormUnitTest 
	${PROJECT_BINARY_DIR}/src/Minimizations/Norms/unittests/LInfinityNormUnitTest)

########################## NormFactoryUnitTest #######################
add_executable (NormFactoryUnitTest 
	${PROJECT_SOURCE_DIR}/src/unittests/TestRunner.cpp 
	LInfinityNormUnitTest.cpp) 
target_link_libraries (NormFactoryUnitTest 
	${EXTRA_LIBS}
	${Boost_LIBRARIES}
	"${CPPUNIT_LIBRARY}") 

add_test (NormFactoryUnitTest 
	${PROJECT_BINARY_DIR}/src/Minimizations/Norms/unittests/NormFactoryUnitTest)
