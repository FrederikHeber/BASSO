# does the application run
add_test (ComputerTomographyConfiguratorFails ${PROJECT_BINARY_DIR}/src/ComputerTomographyBase/ComputerTomographyConfigurator)
set_tests_properties (ComputerTomographyConfiguratorFails
	PROPERTIES WILL_FAIL true)	# return code 255

# does the application run
add_test (ComputerTomographyConfiguratorFails2 
	${PROJECT_BINARY_DIR}/src/ComputerTomographyBase/ComputerTomographyConfigurator
		--type-space-x "lp" 
		--px 2)
set_tests_properties (ComputerTomographyConfiguratorFails2
	PROPERTIES WILL_FAIL true)	# return code 255

# does the application write config files
add_test (ComputerTomographyConfiguratorWritesConfig 
	${PROJECT_BINARY_DIR}/src/ComputerTomographyBase/ComputerTomographyConfigurator
	--config test.conf
	--type-space-x "lp" 
	--px 10
	--verbose 2)
set_tests_properties (ComputerTomographyConfiguratorWritesConfig
	PROPERTIES PASS_REGULAR_EXPRESSION "Writing configuration file test.conf")
