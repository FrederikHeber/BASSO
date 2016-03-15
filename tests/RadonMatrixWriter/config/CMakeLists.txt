# does the application run
add_test (RadonMatrixWriterConfiguratorFails ${PROJECT_BINARY_DIR}/Examples/RadonMatrixWriter/RadonMatrixWriterConfigurator)
set_tests_properties (RadonMatrixWriterConfiguratorFails
	PROPERTIES WILL_FAIL true)	# return code 255

# does the version show
add_test (RadonMatrixWriterConfiguratorVersion
	${PROJECT_BINARY_DIR}/Examples/RadonMatrixWriter/RadonMatrixWriterConfigurator
		--version)
set_tests_properties (RadonMatrixWriterConfiguratorVersion
	PROPERTIES WILL_FAIL false)
set_tests_properties (RadonMatrixWriterConfiguratorVersion
	PROPERTIES PASS_REGULAR_EXPRESSION "version")

# does the application run
add_test (RadonMatrixWriterConfiguratorFails2 
	${PROJECT_BINARY_DIR}/Examples/RadonMatrixWriter/RadonMatrixWriterConfigurator
		--num-pixels-x 16 
		--num-pixels-y 16
		--num-angles 10 
		--num-offsets 11)
set_tests_properties (RadonMatrixWriterConfiguratorFails2
	PROPERTIES WILL_FAIL true)	# return code 255

# does the application write config files
add_test (RadonMatrixWriterConfiguratorWritesConfig 
	${PROJECT_BINARY_DIR}/Examples/RadonMatrixWriter/RadonMatrixWriterConfigurator
	--config test.conf
	--num-pixels-x 16 
	--num-pixels-y 16
	--num-angles 10 
	--num-offsets 11
	--verbose 2)
set_tests_properties (RadonMatrixWriterConfiguratorWritesConfig
	PROPERTIES PASS_REGULAR_EXPRESSION "Writing configuration file test.conf")
