### TESTS ###

########################## LpNormUnitTest #######################
add_executable (LpNormUnitTest 
	${PROJECT_SOURCE_DIR}/src/unittests/TestRunner.cpp 
	LpNormUnitTest.cpp) 
target_link_libraries (LpNormUnitTest 
	${EXTRA_LIBS}
	"${CPPUNIT_LIBRARY}") 

add_test (LpNormUnitTest 
	${PROJECT_BINARY_DIR}/src/Minimizations/Norms/unittests/LpNormUnitTest)
