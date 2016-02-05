# does the application run
add_test (NoiseAdderConfiguratorFails ${PROJECT_BINARY_DIR}/Examples/NoiseAdder/NoiseAdderConfigurator)
set_tests_properties (NoiseAdderConfiguratorFails
	PROPERTIES WILL_FAIL true)	# return code 255

# does the version show
add_test (NoiseAdderConfiguratorVersion
	${PROJECT_BINARY_DIR}/Examples/NoiseAdder/NoiseAdderConfigurator
		--version)
set_tests_properties (NoiseAdderConfiguratorVersion
	PROPERTIES WILL_FAIL false)
set_tests_properties (NoiseAdderConfiguratorVersion
	PROPERTIES PASS_REGULAR_EXPRESSION "version")

# does the application run
add_test (NoiseAdderConfiguratorFails2 
	${PROJECT_BINARY_DIR}/Examples/NoiseAdder/NoiseAdderConfigurator
		--relative-level 2)
set_tests_properties (NoiseAdderConfiguratorFails2
	PROPERTIES WILL_FAIL true)	# return code 255

# does the application write config files
add_test (NoiseAdderConfiguratorWritesConfig 
	${PROJECT_BINARY_DIR}/Examples/NoiseAdder/NoiseAdderConfigurator
		--config test.conf
		--noise-level 0.01 
		--verbose 2)
set_tests_properties (NoiseAdderConfiguratorWritesConfig
	PROPERTIES PASS_REGULAR_EXPRESSION "Writing configuration file test.conf")
