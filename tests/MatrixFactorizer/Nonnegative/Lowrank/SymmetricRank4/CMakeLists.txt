# does simple test call work?
add_test (MatrixFactorizer_Nonnegative_SymmetricRank4
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
set_tests_properties (MatrixFactorizer_Nonnegative_SymmetricRank4
	PROPERTIES WILL_FAIL false)
set_tests_properties (MatrixFactorizer_Nonnegative_SymmetricRank4
	PROPERTIES PASS_REGULAR_EXPRESSION "Norm of .*1.4")

# test 1b: checking results
add_test (MatrixFactorizer_Nonnegative_SymmetricRank4_check_first
	${DIFFCMD} ${NDIFF_OPTIONS} ${CMAKE_CURRENT_SOURCE_DIR}/post/spectralmatrix.m spectralmatrix.m)
set_tests_properties (MatrixFactorizer_Nonnegative_SymmetricRank4_check_first
	PROPERTIES WILL_FAIL false)
set_tests_properties (MatrixFactorizer_Nonnegative_SymmetricRank4_check_first
	PROPERTIES DEPENDS MatrixFactorizer_Nonnegative_SymmetricRank4)

# test 1c: checking results
add_test (MatrixFactorizer_Nonnegative_SymmetricRank4_check_second
	${DIFFCMD} ${NDIFF_OPTIONS} ${CMAKE_CURRENT_SOURCE_DIR}/post/pixelmatrix.m pixelmatrix.m)
set_tests_properties (MatrixFactorizer_Nonnegative_SymmetricRank4_check_second
	PROPERTIES WILL_FAIL false)
set_tests_properties (MatrixFactorizer_Nonnegative_SymmetricRank4_check_second
	PROPERTIES DEPENDS MatrixFactorizer_Nonnegative_SymmetricRank4)

# test 2: same as test 1 but using parallel code
if (MPI_FOUND)
	add_test (MatrixFactorizer_Nonnegative_SymmetricRank4_mpi
		${MPIEXEC} -${MPIEXEC_NUMPROC_FLAG} 2
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
			--solution-first-factor spectralmatrixpi.m
			--solution-second-factor pixelmatrixpi.m)
	set_tests_properties (MatrixFactorizer_Nonnegative_SymmetricRank4_mpi
		PROPERTIES WILL_FAIL false)
	set_tests_properties (MatrixFactorizer_Nonnegative_SymmetricRank4_mpi
		PROPERTIES PASS_REGULAR_EXPRESSION "Norm of .*1.4")
	
	# test 1b: checking results
	add_test (MatrixFactorizer_Nonnegative_SymmetricRank4_mpi_check_first
		${DIFFCMD} ${NDIFF_OPTIONS} ${CMAKE_CURRENT_SOURCE_DIR}/post/spectralmatrix.m spectralmatrixpi.m)
	set_tests_properties (MatrixFactorizer_Nonnegative_SymmetricRank4_mpi_check_first
		PROPERTIES WILL_FAIL false)
	set_tests_properties (MatrixFactorizer_Nonnegative_SymmetricRank4_mpi_check_first
		PROPERTIES DEPENDS MatrixFactorizer_Nonnegative_SymmetricRank4_mpi)
	
	# test 1c: checking results
	add_test (MatrixFactorizer_Nonnegative_SymmetricRank4_mpi_check_second
		${DIFFCMD} ${NDIFF_OPTIONS} ${CMAKE_CURRENT_SOURCE_DIR}/post/pixelmatrix.m pixelmatrixpi.m)
	set_tests_properties (MatrixFactorizer_Nonnegative_SymmetricRank4_mpi_check_second
		PROPERTIES WILL_FAIL false)
	set_tests_properties (MatrixFactorizer_Nonnegative_SymmetricRank4_mpi_check_second
		PROPERTIES DEPENDS MatrixFactorizer_Nonnegative_SymmetricRank4_mpi)
endif (MPI_FOUND)
