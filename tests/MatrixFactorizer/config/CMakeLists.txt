# does the application run
add_test (MatrixFactorizerConfiguratorFails ${PROJECT_BINARY_DIR}/src/MatrixFactorizerConfigurator)
set_tests_properties (MatrixFactorizerConfiguratorFails
	PROPERTIES WILL_FAIL true)	# return code 255

# does the application run
add_test (MatrixFactorizerConfiguratorFails2 
	${PROJECT_BINARY_DIR}/src/MatrixFactorizerConfigurator
		--normx 2)
set_tests_properties (MatrixFactorizerConfiguratorFails2
	PROPERTIES WILL_FAIL true)	# return code 255

# does the application write config files
add_test (MatrixFactorizerConfiguratorWritesConfig 
	${PROJECT_BINARY_DIR}/src/MatrixFactorizerConfigurator
	--config test.conf
	--normx 10
	--verbose 2)
set_tests_properties (MatrixFactorizerConfiguratorWritesConfig
	PROPERTIES PASS_REGULAR_EXPRESSION "Writing configuration file test.conf")
