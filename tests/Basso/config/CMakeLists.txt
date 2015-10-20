# does the application run
add_test (BassoConfiguratorFails ${PROJECT_BINARY_DIR}/Examples/Basso/BassoConfigurator)
set_tests_properties (BassoConfiguratorFails
	PROPERTIES WILL_FAIL true)	# return code 255

# does the application run
add_test (BassoConfiguratorFails2 
	${PROJECT_BINARY_DIR}/Examples/Basso/BassoConfigurator
		--type-space-x "lp" 
		--px 2)
set_tests_properties (BassoConfiguratorFails2
	PROPERTIES WILL_FAIL true)	# return code 255

# does the application write config files
add_test (BassoConfiguratorWritesConfig 
	${PROJECT_BINARY_DIR}/Examples/Basso/BassoConfigurator
	--config test.conf
	--type-space-x "lp" 
	--px 10
	--verbose 2)
set_tests_properties (BassoConfiguratorWritesConfig
	PROPERTIES PASS_REGULAR_EXPRESSION "Writing configuration file test.conf")
