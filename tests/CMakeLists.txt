### TESTS ###
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND}
                  DEPENDS Basso)

# does the application run
add_test (BassoRuns ${PROJECT_BINARY_DIR}/src/Basso)
set_tests_properties (BassoRuns
	PROPERTIES WILL_FAIL true)

# does the usage message work?
add_test (BassoUsage ${PROJECT_BINARY_DIR}/src/Basso --help)
set_tests_properties (BassoUsage
	PROPERTIES WILL_FAIL true)
set_tests_properties (BassoUsage
	PROPERTIES PASS_REGULAR_EXPRESSION "Usage:.*norm")

# does simple test call work?
add_test (BassoSimple 
	${PROJECT_BINARY_DIR}/src/Basso 
		--normx 2 
		--normy 2 
		--powery 2 
		--delta 0.01 
		--matrix ${PROJECT_SOURCE_DIR}/tests/matrix.m 
		--rhs ${PROJECT_SOURCE_DIR}/tests/vector.m)
set_tests_properties (BassoSimple
	PROPERTIES PASS_REGULAR_EXPRESSION "Norm of X.*2")
set_tests_properties (BassoSimple
	PROPERTIES PASS_REGULAR_EXPRESSION "Norm of Y.*2")
set_tests_properties (BassoSimple
	PROPERTIES PASS_REGULAR_EXPRESSION ".*noise.*0.01")
set_tests_properties (BassoSimple
	PROPERTIES PASS_REGULAR_EXPRESSION ".*Ax = y.*1 0")
	
