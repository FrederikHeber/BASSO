# does the application run
add_test (BassoConfiguratorFails ${PROJECT_BINARY_DIR}/src/BassoBase/BassoConfigurator)
set_tests_properties (BassoConfiguratorFails
	PROPERTIES WILL_FAIL true)	# return code 255

# does the application run
add_test (BassoConfiguratorFails2 
	${PROJECT_BINARY_DIR}/src/BassoBase/BassoConfigurator
		--normx 2)
set_tests_properties (BassoConfiguratorFails2
	PROPERTIES WILL_FAIL true)	# return code 255

# does the application write config files
add_test (BassoConfiguratorWritesConfig 
	${PROJECT_BINARY_DIR}/src/BassoBase/BassoConfigurator
	--config test.conf
	--normx 10
	--verbose 2)
set_tests_properties (BassoConfiguratorWritesConfig
	PROPERTIES PASS_REGULAR_EXPRESSION "Writing configuration file test.conf")
