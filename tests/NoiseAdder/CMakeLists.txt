# does the application run
add_test (NoiseAdderRunFails ${PROJECT_BINARY_DIR}/Examples/NoiseAdder/NoiseAdder)
set_tests_properties (NoiseAdderRunFails
	PROPERTIES WILL_FAIL true)	# return code 255

add_test (NoiseAdderRuns ${PROJECT_BINARY_DIR}/Examples/NoiseAdder/NoiseAdder)
set_tests_properties (NoiseAdderRuns
	PROPERTIES PASS_REGULAR_EXPRESSION "There was an error parsing options.*")

# does the usage message work?
add_test (NoiseAdderUsage ${PROJECT_BINARY_DIR}/Examples/NoiseAdder/NoiseAdder --help)
set_tests_properties (NoiseAdderUsage
	PROPERTIES WILL_FAIL true)	# return code 1
	
add_test (NoiseAdderVersionMessage ${PROJECT_BINARY_DIR}/Examples/NoiseAdder/NoiseAdder --version)
set_tests_properties (NoiseAdderVersionMessage
	PROPERTIES PASS_REGULAR_EXPRESSION "version")

add_test (NoiseAdderUsageMessage ${PROJECT_BINARY_DIR}/Examples/NoiseAdder/NoiseAdder --help)
set_tests_properties (NoiseAdderUsageMessage
	PROPERTIES PASS_REGULAR_EXPRESSION "NoiseAdder options:")

add_subdirectory (config)

add_subdirectory (simple)

