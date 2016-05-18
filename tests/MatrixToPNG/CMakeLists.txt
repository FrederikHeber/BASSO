# does the application run
add_test (MatrixToPNGRunFails ${PROJECT_BINARY_DIR}/Examples/MatrixToPNG/MatrixToPNG)
set_tests_properties (MatrixToPNGRunFails
	PROPERTIES WILL_FAIL true)	# return code 255

add_test (MatrixToPNGRuns ${PROJECT_BINARY_DIR}/Examples/MatrixToPNG/MatrixToPNG)
set_tests_properties (MatrixToPNGRuns
	PROPERTIES PASS_REGULAR_EXPRESSION "There was an error parsing options.*")

# does the usage message work?
add_test (MatrixToPNGUsage ${PROJECT_BINARY_DIR}/Examples/MatrixToPNG/MatrixToPNG --help)
set_tests_properties (MatrixToPNGUsage
	PROPERTIES WILL_FAIL true)	# return code 1
	
add_test (MatrixToPNGVersionMessage ${PROJECT_BINARY_DIR}/Examples/MatrixToPNG/MatrixToPNG --version)
set_tests_properties (MatrixToPNGVersionMessage
	PROPERTIES PASS_REGULAR_EXPRESSION "version")

add_test (MatrixToPNGUsageMessage ${PROJECT_BINARY_DIR}/Examples/MatrixToPNG/MatrixToPNG --help)
set_tests_properties (MatrixToPNGUsageMessage
	PROPERTIES PASS_REGULAR_EXPRESSION "MatrixToPNG options:")

add_subdirectory (colortables)
