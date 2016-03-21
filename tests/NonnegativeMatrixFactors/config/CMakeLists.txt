# does the application run
add_test (NonnegativeMatrixFactorsConfiguratorFails ${PROJECT_BINARY_DIR}/Examples/NonnegativeMatrixFactors/NonnegativeMatrixFactorsConfigurator)
set_tests_properties (NonnegativeMatrixFactorsConfiguratorFails
	PROPERTIES WILL_FAIL true)	# return code 255

# does the application write config files
add_test (NonnegativeMatrixFactorsConfiguratorWritesConfig 
	${PROJECT_BINARY_DIR}/Examples/NonnegativeMatrixFactors/NonnegativeMatrixFactorsConfigurator
	--config test.conf
	--verbose 2
	--matrix ${CMAKE_CURRENT_SOURCE_DIR}/pre/matrix.mat
	--truncation-dimension 3
	--destination-first-factor spectralmatrix.mat
	--destination-second-factor pixelmatrix.mat
	)
set_tests_properties (NonnegativeMatrixFactorsConfiguratorWritesConfig
	PROPERTIES PASS_REGULAR_EXPRESSION "Writing configuration file test.conf")
