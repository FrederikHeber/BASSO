# does the application run
add_test (GravityConfiguratorFails ${PROJECT_BINARY_DIR}/Examples/Gravity/GravityConfigurator)
set_tests_properties (GravityConfiguratorFails
	PROPERTIES WILL_FAIL true)	# return code 255

# does the version show
add_test (GravityConfiguratorVersion
	${PROJECT_BINARY_DIR}/Examples/Gravity/GravityConfigurator
		--version)
set_tests_properties (GravityConfiguratorVersion
	PROPERTIES WILL_FAIL false)
set_tests_properties (GravityConfiguratorVersion
	PROPERTIES PASS_REGULAR_EXPRESSION "version")

# does the application write config files
add_test (GravityConfiguratorWritesConfig 
	${PROJECT_BINARY_DIR}/Examples/Gravity/GravityConfigurator
	--config test.conf
	--verbose 2
	--type-space-x "lp" 
	--px 2 
	--type-space-y "lp" 
	--py 2 
	--powery 2 
	--maxiter 10
	--delta 0.01 
	--depth 0.5
	--discretization 10
	--density-file density.m)
set_tests_properties (GravityConfiguratorWritesConfig
	PROPERTIES PASS_REGULAR_EXPRESSION "Writing configuration file test.conf")
