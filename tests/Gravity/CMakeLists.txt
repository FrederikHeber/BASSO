# does the application run
add_test (GravityRunFails ${PROJECT_BINARY_DIR}/Examples/Gravity/Gravity)
set_tests_properties (GravityRunFails
	PROPERTIES WILL_FAIL true)	# return code 255

add_test (GravityRuns ${PROJECT_BINARY_DIR}/Examples/Gravity/Gravity)
set_tests_properties (GravityRuns
	PROPERTIES PASS_REGULAR_EXPRESSION "There was an error parsing options.*")

# does the usage message work?
add_test (GravityUsage ${PROJECT_BINARY_DIR}/Examples/Gravity/Gravity --help)
set_tests_properties (GravityUsage
	PROPERTIES WILL_FAIL true)	# return code 1
	
add_test (GravityVersionMessage ${PROJECT_BINARY_DIR}/Examples/Gravity/Gravity --version)
set_tests_properties (GravityVersionMessage
	PROPERTIES PASS_REGULAR_EXPRESSION "version")

add_test (GravityUsageMessage ${PROJECT_BINARY_DIR}/Examples/Gravity/Gravity --help)
set_tests_properties (GravityUsageMessage
	PROPERTIES PASS_REGULAR_EXPRESSION "General options:")

add_subdirectory (simple)
