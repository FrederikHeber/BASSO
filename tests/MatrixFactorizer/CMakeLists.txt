# does the application run
add_test (MatrixFactorizerRunFails
	${PROJECT_BINARY_DIR}/Examples/MatrixFactorizer/MatrixFactorizer)
set_tests_properties (MatrixFactorizerRunFails
	PROPERTIES WILL_FAIL true)	# return code 255

add_test (MatrixFactorizerRuns 
	${PROJECT_BINARY_DIR}/Examples/MatrixFactorizer/MatrixFactorizer)
set_tests_properties (MatrixFactorizerRuns
	PROPERTIES PASS_REGULAR_EXPRESSION "There was an error parsing options.*")

# does the usage message work?
add_test (MatrixFactorizerUsage 
	${PROJECT_BINARY_DIR}/Examples/MatrixFactorizer/MatrixFactorizer --help)
set_tests_properties (MatrixFactorizerUsage
	PROPERTIES WILL_FAIL true)	# return code 1
	
add_test (MatrixFactorizerUsageMessage 
	${PROJECT_BINARY_DIR}/Examples/MatrixFactorizer/MatrixFactorizer --help)
set_tests_properties (MatrixFactorizerUsageMessage
	PROPERTIES PASS_REGULAR_EXPRESSION "General options:")

# configuration file tests
add_subdirectory (config)

# various symmetric rank decompositions
add_subdirectory (SymmetricRank1)
add_subdirectory (SymmetricRank2)
add_subdirectory (SymmetricRank4)

# test MatrixFactorizer with mpirun
if (MPI_FOUND)
	add_subdirectory (mpirun)
endif (MPI_FOUND)