# does the application run
add_test (BassoRunFails ${PROJECT_BINARY_DIR}/src/Basso)
set_tests_properties (BassoRunFails
	PROPERTIES WILL_FAIL true)	# return code 255

add_test (BassoRuns ${PROJECT_BINARY_DIR}/src/Basso)
set_tests_properties (BassoRuns
	PROPERTIES PASS_REGULAR_EXPRESSION "There was an error parsing options.*")

# does the usage message work?
add_test (BassoUsage ${PROJECT_BINARY_DIR}/src/Basso --help)
set_tests_properties (BassoUsage
	PROPERTIES WILL_FAIL true)	# return code 1
	
add_test (BassoUsageMessage ${PROJECT_BINARY_DIR}/src/Basso --help)
set_tests_properties (BassoUsageMessage
	PROPERTIES PASS_REGULAR_EXPRESSION "General options:")

add_subdirectory (simple)

# tests for Landweber algorithm
add_subdirectory (Landweber)

# tests for SESOP algorithm
add_subdirectory (SESOP)

# tests for RESESOP algorithm
add_subdirectory (RESESOP)

