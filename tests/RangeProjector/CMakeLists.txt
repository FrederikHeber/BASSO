# does the application run
add_test (RangeProjectorRunFails ${PROJECT_BINARY_DIR}/Examples/RangeProjector/RangeProjector)
set_tests_properties (RangeProjectorRunFails
	PROPERTIES WILL_FAIL true)	# return code 255

add_test (RangeProjectorRuns ${PROJECT_BINARY_DIR}/Examples/RangeProjector/RangeProjector)
set_tests_properties (RangeProjectorRuns
	PROPERTIES PASS_REGULAR_EXPRESSION "There was an error parsing options.*")

# does the usage message work?
add_test (RangeProjectorUsage ${PROJECT_BINARY_DIR}/Examples/RangeProjector/RangeProjector --help)
set_tests_properties (RangeProjectorUsage
	PROPERTIES WILL_FAIL true)	# return code 1
	
add_test (RangeProjectorVersionMessage ${PROJECT_BINARY_DIR}/Examples/RangeProjector/RangeProjector --version)
set_tests_properties (RangeProjectorVersionMessage
	PROPERTIES PASS_REGULAR_EXPRESSION "version")

add_test (RangeProjectorUsageMessage ${PROJECT_BINARY_DIR}/Examples/RangeProjector/RangeProjector --help)
set_tests_properties (RangeProjectorUsageMessage
	PROPERTIES PASS_REGULAR_EXPRESSION "General options:")

add_subdirectory (config)

add_subdirectory (simple)
