# does the application run
add_test (NonnegativeMatrixFactorsConfiguratorFails ${PROJECT_BINARY_DIR}/Examples/NonnegativeMatrixFactors/NonnegativeMatrixFactorsConfigurator)
set_tests_properties (NonnegativeMatrixFactorsConfiguratorFails
	PROPERTIES WILL_FAIL true)	# return code 255

# does the application write config files
add_test (NonnegativeMatrixFactorsConfiguratorWritesConfig 
	${PROJECT_BINARY_DIR}/Examples/NonnegativeMatrixFactors/NonnegativeMatrixFactorsConfigurator
	--config test.conf
	--verbose 2
	--source-first-factor ${CMAKE_CURRENT_SOURCE_DIR}/pre/spectralmatrix.m
	--source-second-factor ${CMAKE_CURRENT_SOURCE_DIR}/pre/pixelmatrix.m
	--destination-first-factor spectralmatrix.m
	--destination-second-factor pixelmatrix.m
	)
set_tests_properties (NonnegativeMatrixFactorsConfiguratorWritesConfig
	PROPERTIES PASS_REGULAR_EXPRESSION "Writing configuration file test.conf")
