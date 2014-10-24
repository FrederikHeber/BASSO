### TESTS ###

########################## LpDualityMappingUnitTest #######################
add_executable (LpDualityMappingUnitTest 
	${PROJECT_SOURCE_DIR}/src/unittests/TestRunner.cpp 
	LpDualityMappingUnitTest.cpp) 
target_link_libraries (LpDualityMappingUnitTest 
	${EXTRA_LIBS}
	"${CPPUNIT_LIBRARY}")

add_test (LpDualityMappingUnitTest 
	${PROJECT_BINARY_DIR}/src/Minimizations/Mappings/unittests/LpDualityMappingUnitTest)
