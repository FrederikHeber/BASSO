# does the application run
add_test (NonnegativeMatrixFactorsRunFails ${PROJECT_BINARY_DIR}/Examples/NonnegativeMatrixFactors/NonnegativeMatrixFactors)
set_tests_properties (NonnegativeMatrixFactorsRunFails
	PROPERTIES WILL_FAIL true)	# return code 255

add_test (NonnegativeMatrixFactorsRuns ${PROJECT_BINARY_DIR}/Examples/NonnegativeMatrixFactors/NonnegativeMatrixFactors)
set_tests_properties (NonnegativeMatrixFactorsRuns
	PROPERTIES PASS_REGULAR_EXPRESSION "There was an error parsing options.*")

# does the usage message work?
add_test (NonnegativeMatrixFactorsUsage ${PROJECT_BINARY_DIR}/Examples/NonnegativeMatrixFactors/NonnegativeMatrixFactors --help)
set_tests_properties (NonnegativeMatrixFactorsUsage
	PROPERTIES WILL_FAIL true)	# return code 1
	
add_test (NonnegativeMatrixFactorsUsageMessage ${PROJECT_BINARY_DIR}/Examples/NonnegativeMatrixFactors/NonnegativeMatrixFactors --help)
set_tests_properties (NonnegativeMatrixFactorsUsageMessage
	PROPERTIES PASS_REGULAR_EXPRESSION "NonnegativeMatrixFactors options:")

add_subdirectory (config)

add_subdirectory (simple)
