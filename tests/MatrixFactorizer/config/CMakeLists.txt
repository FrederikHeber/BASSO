# does the application run
add_test (MatrixFactorizerConfiguratorFails ${PROJECT_BINARY_DIR}/Examples/MatrixFactorizer/MatrixFactorizerConfigurator)
set_tests_properties (MatrixFactorizerConfiguratorFails
	PROPERTIES WILL_FAIL true)	# return code 255

# does the version show
add_test (MatrixFactorizerConfiguratorVersion
	${PROJECT_BINARY_DIR}/Examples/MatrixFactorizer/MatrixFactorizerConfigurator
		--version)
set_tests_properties (MatrixFactorizerConfiguratorVersion
	PROPERTIES WILL_FAIL false)
set_tests_properties (MatrixFactorizerConfiguratorVersion
	PROPERTIES PASS_REGULAR_EXPRESSION "version")

# does the application run
add_test (MatrixFactorizerConfiguratorFails2 
	${PROJECT_BINARY_DIR}/Examples/MatrixFactorizer/MatrixFactorizerConfigurator
		--type-space-x "lp" 
		--px 2)
set_tests_properties (MatrixFactorizerConfiguratorFails2
	PROPERTIES WILL_FAIL true)	# return code 255

# does the application write config files
add_test (MatrixFactorizerConfiguratorWritesConfig 
	${PROJECT_BINARY_DIR}/Examples/MatrixFactorizer/MatrixFactorizerConfigurator
	--config test.conf
	--type-space-x "lp" 
	--px 10
	--verbose 2)
set_tests_properties (MatrixFactorizerConfiguratorWritesConfig
	PROPERTIES PASS_REGULAR_EXPRESSION "Writing configuration file test.conf")
