# test 1: Matrix factorization with rank 1 
add_test (MatrixFactorizer_Nonnegative_SymmetricRank1 
	${PROJECT_BINARY_DIR}/Examples/MatrixFactorizer/MatrixFactorizer 
		--type-space-x "lp" 
		--px 2 
		--type-space-y "lp" 
		--py 2 
		--powery 2 
		--delta 1e-6 
		--maxiter 50
		--max-loops 11
		--orthogonal-directions 1
		--number-directions 1
		--sparse-dim 1
		--residual-threshold 1e-4
		--auxiliary-constraints "Nonnegative"
		--max-sfp-loops 2
		--data ${CMAKE_CURRENT_SOURCE_DIR}/pre/data.m 
		--solution-first-factor spectralmatrix.m
		--solution-second-factor pixelmatrix.m)
set_tests_properties (MatrixFactorizer_Nonnegative_SymmetricRank1
	PROPERTIES WILL_FAIL false)
set_tests_properties (MatrixFactorizer_Nonnegative_SymmetricRank1
	PROPERTIES PASS_REGULAR_EXPRESSION "Norm of .*e-09")

# test 1b: checking results
add_test (MatrixFactorizer_Nonnegative_SymmetricRank1_check_first
	${DIFFCMD} ${NDIFF_OPTIONS} ${CMAKE_CURRENT_SOURCE_DIR}/post/spectralmatrix.m spectralmatrix.m)
set_tests_properties (MatrixFactorizer_Nonnegative_SymmetricRank1_check_first
	PROPERTIES WILL_FAIL false)
set_tests_properties (MatrixFactorizer_Nonnegative_SymmetricRank1_check_first
	PROPERTIES DEPENDS MatrixFactorizer_Nonnegative_SymmetricRank1)

# test 1c: checking results
add_test (MatrixFactorizer_Nonnegative_SymmetricRank1_check_second
	${DIFFCMD} ${NDIFF_OPTIONS} ${CMAKE_CURRENT_SOURCE_DIR}/post/pixelmatrix.m pixelmatrix.m)
set_tests_properties (MatrixFactorizer_Nonnegative_SymmetricRank1_check_second
	PROPERTIES WILL_FAIL false)
set_tests_properties (MatrixFactorizer_Nonnegative_SymmetricRank1_check_second
	PROPERTIES DEPENDS MatrixFactorizer_Nonnegative_SymmetricRank1)

# test 2: Matrix factorization with rank 1 but more than one search direction
add_test (MatrixFactorizer_Nonnegative_SymmetricRank1_fail 
	${PROJECT_BINARY_DIR}/Examples/MatrixFactorizer/MatrixFactorizer 
		--type-space-x "lp" 
		--px 2 
		--type-space-y "lp" 
		--py 2 
		--powery 2 
		--delta 1e-6
		--maxiter 50
		--max-loops 11
		--orthogonal-directions 1
		--number-directions 2
		--sparse-dim 1
		--residual-threshold 1e-4
		--auxiliary-constraints "Nonnegative"
		--max-sfp-loops 2
		--data ${CMAKE_CURRENT_SOURCE_DIR}/pre/data.m)
set_tests_properties (MatrixFactorizer_Nonnegative_SymmetricRank1_fail
	PROPERTIES WILL_FAIL true)
	
# test 3: same as test 1 but using parallel code
if (MPI_FOUND)
	add_test (MatrixFactorizer_Nonnegative_SymmetricRank1_mpi
		${MPIEXEC} -${MPIEXEC_NUMPROC_FLAG} 2
		${PROJECT_BINARY_DIR}/Examples/MatrixFactorizer/MatrixFactorizer
			--verbose 0
			--type-space-x "lp" 
			--px 2 
			--type-space-y "lp" 
			--py 2 
			--powery 2 
			--delta 1e-6
			--maxiter 50
			--max-loops 11
			--orthogonal-directions 1
			--number-directions 1
			--sparse-dim 1
			--residual-threshold 1e-4
			--auxiliary-constraints "Nonnegative"
			--max-sfp-loops 2
			--data ${CMAKE_CURRENT_SOURCE_DIR}/pre/data.m 
			--solution-first-factor spectralmatrix_mpi.m
			--solution-second-factor pixelmatrix_mpi.m)
	set_tests_properties (MatrixFactorizer_Nonnegative_SymmetricRank1_mpi
		PROPERTIES WILL_FAIL false)
	set_tests_properties (MatrixFactorizer_Nonnegative_SymmetricRank1_mpi
		PROPERTIES PASS_REGULAR_EXPRESSION "Norm of .*e-09")

	# test 3b: checking results
	add_test (MatrixFactorizer_Nonnegative_SymmetricRank1_mpi_check_first
		${DIFFCMD} ${NDIFF_OPTIONS} ${CMAKE_CURRENT_SOURCE_DIR}/post/spectralmatrix.m spectralmatrix_mpi.m)
	set_tests_properties (MatrixFactorizer_Nonnegative_SymmetricRank1_mpi_check_first
		PROPERTIES WILL_FAIL false)
	set_tests_properties (MatrixFactorizer_Nonnegative_SymmetricRank1_mpi_check_first
		PROPERTIES DEPENDS MatrixFactorizer_Nonnegative_SymmetricRank1_mpi)
	
	# test 3c: checking results
	add_test (MatrixFactorizer_Nonnegative_SymmetricRank1_mpi_check_second
		${DIFFCMD} ${NDIFF_OPTIONS} ${CMAKE_CURRENT_SOURCE_DIR}/post/pixelmatrix.m pixelmatrix_mpi.m)
	set_tests_properties (MatrixFactorizer_Nonnegative_SymmetricRank1_mpi_check_second
		PROPERTIES WILL_FAIL false)
	set_tests_properties (MatrixFactorizer_Nonnegative_SymmetricRank1_mpi_check_second
		PROPERTIES DEPENDS MatrixFactorizer_Nonnegative_SymmetricRank1_mpi)
endif (MPI_FOUND)
