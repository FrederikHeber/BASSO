# test 1: Matrix factorization with rank 1 
add_test (MatrixFactorizer_Indefinite_NonsymmetricRank1 
	${PROJECT_BINARY_DIR}/Examples/MatrixFactorizer/MatrixFactorizer 
		--type-space-x "lp" 
		--px 2 
		--type-space-y "lp" 
		--py 2 
		--powery 2 
		--delta 1e-6
		--maxiter 50
		--orthogonal-directions 1
		--number-directions 1
		--sparse-dim 1
		--residual-threshold 1e-4
		--data ${CMAKE_CURRENT_SOURCE_DIR}/pre/data.m 
		--solution-first-factor first-factor.m
		--solution-second-factor second-factor.m)
set_tests_properties (MatrixFactorizer_Indefinite_NonsymmetricRank1
	PROPERTIES WILL_FAIL false)
set_tests_properties (MatrixFactorizer_Indefinite_NonsymmetricRank1
	PROPERTIES PASS_REGULAR_EXPRESSION "Norm of .*e-09")

# test 1b: checking results
add_test (MatrixFactorizer_Indefinite_NonsymmetricRank1_check_first
	${DIFFCMD} ${NDIFF_OPTIONS} ${CMAKE_CURRENT_SOURCE_DIR}/post/first-factor.m first-factor.m)
set_tests_properties (MatrixFactorizer_Indefinite_NonsymmetricRank1_check_first
	PROPERTIES WILL_FAIL false)
set_tests_properties (MatrixFactorizer_Indefinite_NonsymmetricRank1_check_first
	PROPERTIES DEPENDS MatrixFactorizer_Indefinite_NonsymmetricRank1)

# test 1c: checking results
add_test (MatrixFactorizer_Indefinite_NonsymmetricRank1_check_second
	${DIFFCMD} ${NDIFF_OPTIONS} ${CMAKE_CURRENT_SOURCE_DIR}/post/second-factor.m second-factor.m)
set_tests_properties (MatrixFactorizer_Indefinite_NonsymmetricRank1_check_second
	PROPERTIES WILL_FAIL false)
set_tests_properties (MatrixFactorizer_Indefinite_NonsymmetricRank1_check_second
	PROPERTIES DEPENDS MatrixFactorizer_Indefinite_NonsymmetricRank1)

# test 2: Matrix factorization with rank 1 but more than one search direction
add_test (MatrixFactorizer_Indefinite_NonsymmetricRank1_fail 
	${PROJECT_BINARY_DIR}/Examples/MatrixFactorizer/MatrixFactorizer 
		--type-space-x "lp" 
		--px 2 
		--type-space-y "lp" 
		--py 2 
		--powery 2 
		--delta 1e-6 
		--maxiter 50
		--orthogonal-directions 1
		--number-directions 2
		--sparse-dim 1
		--residual-threshold 1e-4
		--data ${CMAKE_CURRENT_SOURCE_DIR}/pre/data.m)
set_tests_properties (MatrixFactorizer_Indefinite_NonsymmetricRank1_fail
	PROPERTIES WILL_FAIL true)
	
# test 3: same as test 1 but using parallel code
if (MPI_FOUND)
	add_test (MatrixFactorizer_Indefinite_NonsymmetricRank1_mpi
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
			--orthogonal-directions 1
			--max-loops 11
			--number-directions 1
			--sparse-dim 1
			--residual-threshold 1e-4
			--data ${CMAKE_CURRENT_SOURCE_DIR}/pre/data.m 
			--solution-first-factor first-factor_mpi.m
			--solution-second-factor second-factor_mpi.m)
	set_tests_properties (MatrixFactorizer_Indefinite_NonsymmetricRank1
		PROPERTIES WILL_FAIL false)
	set_tests_properties (MatrixFactorizer_Indefinite_NonsymmetricRank1
		PROPERTIES PASS_REGULAR_EXPRESSION "Norm of .*e-09")

	# test 3b: checking results
	add_test (MatrixFactorizer_Indefinite_NonsymmetricRank1_mpi_check_first
		${DIFFCMD} ${NDIFF_OPTIONS} ${CMAKE_CURRENT_SOURCE_DIR}/post/first-factor.m first-factor_mpi.m)
	set_tests_properties (MatrixFactorizer_Indefinite_NonsymmetricRank1_mpi_check_first
		PROPERTIES WILL_FAIL false)
	set_tests_properties (MatrixFactorizer_Indefinite_NonsymmetricRank1_mpi_check_first
		PROPERTIES DEPENDS MatrixFactorizer_Indefinite_NonsymmetricRank1_mpi)
	
	# test 3c: checking results
	add_test (MatrixFactorizer_Indefinite_NonsymmetricRank1_mpi_check_second
		${DIFFCMD} ${NDIFF_OPTIONS} ${CMAKE_CURRENT_SOURCE_DIR}/post/second-factor.m second-factor_mpi.m)
	set_tests_properties (MatrixFactorizer_Indefinite_NonsymmetricRank1_mpi_check_second
		PROPERTIES WILL_FAIL false)
	set_tests_properties (MatrixFactorizer_Indefinite_NonsymmetricRank1_mpi_check_second
		PROPERTIES DEPENDS MatrixFactorizer_Indefinite_NonsymmetricRank1_mpi)
endif (MPI_FOUND)