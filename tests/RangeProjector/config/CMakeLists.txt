# does the application run
add_test (RangeProjectorConfiguratorFails ${PROJECT_BINARY_DIR}/Examples/RangeProjector/RangeProjectorConfigurator)
set_tests_properties (RangeProjectorConfiguratorFails
	PROPERTIES WILL_FAIL true)	# return code 255

# does the version show
add_test (RangeProjectorConfiguratorVersion
	${PROJECT_BINARY_DIR}/Examples/RangeProjector/RangeProjectorConfigurator
		--version)
set_tests_properties (RangeProjectorConfiguratorVersion
	PROPERTIES WILL_FAIL false)
set_tests_properties (RangeProjectorConfiguratorVersion
	PROPERTIES PASS_REGULAR_EXPRESSION "version")

# does the application run
add_test (RangeProjectorConfiguratorFails2 
	${PROJECT_BINARY_DIR}/Examples/RangeProjector/RangeProjectorConfigurator
		--type-space-x "lp" 
		--px 2)
set_tests_properties (RangeProjectorConfiguratorFails2
	PROPERTIES WILL_FAIL true)	# return code 255

# does the application write config files
add_test (RangeProjectorConfiguratorWritesConfig 
	${PROJECT_BINARY_DIR}/Examples/RangeProjector/RangeProjectorConfigurator
	--config test.conf
	--type-space-x "lp" 
	--px 10
	--verbose 2)
set_tests_properties (RangeProjectorConfiguratorWritesConfig
	PROPERTIES PASS_REGULAR_EXPRESSION "Writing configuration file test.conf")
