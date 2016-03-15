# does the application run
add_test (RadonMatrixWriterRunFails ${PROJECT_BINARY_DIR}/Examples/RadonMatrixWriter/RadonMatrixWriter)
set_tests_properties (RadonMatrixWriterRunFails
	PROPERTIES WILL_FAIL true)	# return code 255

add_test (RadonMatrixWriterRuns ${PROJECT_BINARY_DIR}/Examples/RadonMatrixWriter/RadonMatrixWriter)
set_tests_properties (RadonMatrixWriterRuns
	PROPERTIES PASS_REGULAR_EXPRESSION "There was an error parsing options.*")

# does the usage message work?
add_test (RadonMatrixWriterUsage ${PROJECT_BINARY_DIR}/Examples/RadonMatrixWriter/RadonMatrixWriter --help)
set_tests_properties (RadonMatrixWriterUsage
	PROPERTIES WILL_FAIL true)	# return code 1
	
add_test (RadonMatrixWriterVersionMessage ${PROJECT_BINARY_DIR}/Examples/RadonMatrixWriter/RadonMatrixWriter --version)
set_tests_properties (RadonMatrixWriterVersionMessage
	PROPERTIES PASS_REGULAR_EXPRESSION "version")

add_test (RadonMatrixWriterUsageMessage ${PROJECT_BINARY_DIR}/Examples/RadonMatrixWriter/RadonMatrixWriter --help)
set_tests_properties (RadonMatrixWriterUsageMessage
	PROPERTIES PASS_REGULAR_EXPRESSION "RadonMatrixWriter options:")

add_subdirectory (simple)

# configuration file tests
add_subdirectory (config)
