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

add_subdirectory (simple)

