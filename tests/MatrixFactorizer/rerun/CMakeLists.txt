configure_file(${CMAKE_CURRENT_SOURCE_DIR}/pre/first-factor.m ${CMAKE_CURRENT_BINARY_DIR}/first-factor.m COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/pre/second-factor.m ${CMAKE_CURRENT_BINARY_DIR}/second-factor.m COPYONLY)

# test 1: Matrix factorization with rank 1 but parse in present factors
add_test (
	NAME MatrixFactorizer_Rerunning
	COMMAND ${PROJECT_BINARY_DIR}/Examples/MatrixFactorizer/MatrixFactorizer
	ARGS 
		--type-space-x "lp" 
		--px 2 
		--type-space-y "lp" 
		--py 2 
		--powery 2 
		--delta 0.01 
		--maxiter 50
		--number-directions 1
		--sparse-dim 1
		--data ${CMAKE_CURRENT_SOURCE_DIR}/pre/data.m
		--parse-factors 1 
		--solution-first-factor first-factor.m
		--solution-second-factor second-factor.m
		)
set_tests_properties (MatrixFactorizer_Rerunning
	PROPERTIES WILL_FAIL false)
# we only get 1e-6 due to limited output precision in written files
set_tests_properties (MatrixFactorizer_Rerunning
	PROPERTIES PASS_REGULAR_EXPRESSION "Norm of .*e-06")

	# test 1b: checking results
add_test (MatrixFactorizer_Rerunning_check_first
	${DIFFCMD} ${NDIFF_OPTIONS} ${CMAKE_CURRENT_SOURCE_DIR}/post/first-factor.m first-factor.m)
set_tests_properties (MatrixFactorizer_Rerunning_check_first
	PROPERTIES WILL_FAIL true)
set_tests_properties (MatrixFactorizer_Rerunning_check_first
	PROPERTIES DEPENDS MatrixFactorizer_Rerunning)

# test 1c: checking results
add_test (MatrixFactorizer_Rerunning_check_second
	${DIFFCMD} ${NDIFF_OPTIONS} ${CMAKE_CURRENT_SOURCE_DIR}/post/second-factor.m second-factor.m)
set_tests_properties (MatrixFactorizer_Rerunning_check_second
	PROPERTIES WILL_FAIL true)
set_tests_properties (MatrixFactorizer_Rerunning_check_second
	PROPERTIES DEPENDS MatrixFactorizer_Rerunning)

# test 2: Matrix factorization with rank 1 but do not parse in present factors
add_test (
	NAME MatrixFactorizer_Rerunning_noparse
	COMMAND ${PROJECT_BINARY_DIR}/Examples/MatrixFactorizer/MatrixFactorizer
	ARGS 
		--type-space-x "lp" 
		--px 2 
		--type-space-y "lp" 
		--py 2 
		--powery 2 
		--delta 0.01 
		--maxiter 50
		--number-directions 1
		--sparse-dim 1
		--data ${CMAKE_CURRENT_SOURCE_DIR}/pre/data.m
		--parse-factors 0
		--solution-first-factor first-factor.m
		--solution-second-factor second-factor.m
		)
set_tests_properties (MatrixFactorizer_Rerunning_noparse
	PROPERTIES WILL_FAIL false)
# we only get 1e-6 due to limited output precision in written files
set_tests_properties (MatrixFactorizer_Rerunning_noparse
	PROPERTIES PASS_REGULAR_EXPRESSION "Norm of .*e-08")
	