### TESTS ###
add_executable (DualityMappingUnitTest 
	${PROJECT_SOURCE_DIR}/src/unittests/TestRunner.cpp 
	DualityMappingUnitTest.cpp) 
target_link_libraries (DualityMappingUnitTest 
	${EXTRA_LIBS}
	"${CPPUNIT_LIBRARY}")

# does the application run
add_test (DualityMappingUnitTest 
	${PROJECT_BINARY_DIR}/src/Minimizations/unittests/DualityMappingUnitTest)
set_tests_properties (DualityMappingUnitTest
	PROPERTIES WILL_FAIL true)
	