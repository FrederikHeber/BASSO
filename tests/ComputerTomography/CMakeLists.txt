# does the application run
add_test (ComputerTomographyRunFails ${PROJECT_BINARY_DIR}/src/ComputerTomographyBase/ComputerTomography)
set_tests_properties (ComputerTomographyRunFails
	PROPERTIES WILL_FAIL true)	# return code 255

add_test (ComputerTomographyRuns ${PROJECT_BINARY_DIR}/src/ComputerTomographyBase/ComputerTomography)
set_tests_properties (ComputerTomographyRuns
	PROPERTIES PASS_REGULAR_EXPRESSION "There was an error parsing options.*")

# does the usage message work?
add_test (ComputerTomographyUsage ${PROJECT_BINARY_DIR}/src/ComputerTomographyBase/ComputerTomography --help)
set_tests_properties (ComputerTomographyUsage
	PROPERTIES WILL_FAIL true)	# return code 1
	
add_test (ComputerTomographyUsageMessage ${PROJECT_BINARY_DIR}/src/ComputerTomographyBase/ComputerTomography --help)
set_tests_properties (ComputerTomographyUsageMessage
	PROPERTIES PASS_REGULAR_EXPRESSION "General options:")

add_subdirectory (simple)

# configuration file tests
add_subdirectory (config)

# SheppLogan phantom tests
add_subdirectory (SheppLogan)
