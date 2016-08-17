configure_file(${CMAKE_CURRENT_SOURCE_DIR}/pre/spectralmatrix.m ${CMAKE_CURRENT_BINARY_DIR}/spectralmatrix_rank1_nonnegative.m COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/pre/pixelmatrix.m ${CMAKE_CURRENT_BINARY_DIR}/pixelmatrix_rank1_nonnegative.m COPYONLY)

# test 1: Matrix factorization with rank 1 but parse in present factors
add_test (
	NAME MatrixFactorizer_Nonnegative_Rerunning
	COMMAND ${PROJECT_BINARY_DIR}/Examples/MatrixFactorizer/MatrixFactorizer
	ARGS 
		--type-space-x "lp" 
		--px 2 
		--type-space-y "lp" 
		--py 2 
		--powery 2 
		--delta 1e-6
		--projection-delta 1e-6
		--maxiter 50
		--orthogonal-directions 1
		--number-directions 1
		--sparse-dim 1
		--residual-threshold 1e-4
		--auxiliary-constraints "Nonnegative"
		--max-sfp-loops 2
		--data ${CMAKE_CURRENT_SOURCE_DIR}/pre/data.m
		--parse-factors 1 
		--solution-first-factor spectralmatrix_rank1_nonnegative.m
		--solution-second-factor pixelmatrix_rank1_nonnegative.m
		)
set_tests_properties (MatrixFactorizer_Nonnegative_Rerunning
	PROPERTIES WILL_FAIL false)
# we only get 1e-6 due to limited output precision in written files
set_tests_properties (MatrixFactorizer_Nonnegative_Rerunning
	PROPERTIES PASS_REGULAR_EXPRESSION "Norm of .*e-06")

	# test 1b: checking results
add_test (MatrixFactorizer_Nonnegative_Rerunning_check_first
	${DIFFCMD} ${NDIFF_OPTIONS} ${CMAKE_CURRENT_SOURCE_DIR}/post/spectralmatrix.m spectralmatrix_rank1_nonnegative.m)
set_tests_properties (MatrixFactorizer_Nonnegative_Rerunning_check_first
	PROPERTIES WILL_FAIL true)
set_tests_properties (MatrixFactorizer_Nonnegative_Rerunning_check_first
	PROPERTIES DEPENDS MatrixFactorizer_Nonnegative_Rerunning)

# test 1c: checking results
add_test (MatrixFactorizer_Nonnegative_Rerunning_check_second
	${DIFFCMD} ${NDIFF_OPTIONS} ${CMAKE_CURRENT_SOURCE_DIR}/post/pixelmatrix.m pixelmatrix_rank1_nonnegative.m)
set_tests_properties (MatrixFactorizer_Nonnegative_Rerunning_check_second
	PROPERTIES WILL_FAIL true)
set_tests_properties (MatrixFactorizer_Nonnegative_Rerunning_check_second
	PROPERTIES DEPENDS MatrixFactorizer_Nonnegative_Rerunning)

# test 2: Matrix factorization with rank 1 but do not parse in present factors
add_test (
	NAME MatrixFactorizer_Nonnegative_Rerunning_noparse
	COMMAND ${PROJECT_BINARY_DIR}/Examples/MatrixFactorizer/MatrixFactorizer
	ARGS 
		--type-space-x "lp" 
		--px 2 
		--type-space-y "lp" 
		--py 2 
		--powery 2 
		--delta 1e-6
		--projection-delta 1e-6
		--maxiter 50
		--orthogonal-directions 1
		--max-loops 10
		--number-directions 1
		--sparse-dim 1
		--residual-threshold 1e-4
		--data ${CMAKE_CURRENT_SOURCE_DIR}/pre/data.m
		--parse-factors 0
		--solution-first-factor spectralmatrix_rank1_nonnegative.m
		--solution-second-factor pixelmatrix_rank1_nonnegative.m
		)
set_tests_properties (MatrixFactorizer_Nonnegative_Rerunning_noparse
	PROPERTIES WILL_FAIL false)
set_tests_properties (MatrixFactorizer_Nonnegative_Rerunning_noparse
	PROPERTIES PASS_REGULAR_EXPRESSION "Norm of .*e-06")
