# does simple test call work?
add_test (MatrixFactorizer_Indefinite_SymmetricRank3
	${PROJECT_BINARY_DIR}/Examples/MatrixFactorizer/MatrixFactorizer 
		--type-space-x "lp" 
		--px 2 
		--type-space-y "lp" 
		--py 2 
		--powery 2 
		--delta 1e-6 
		--projection-delta 1e-6
		--maxiter 50
		--max-loops 11
		--orthogonal-directions 1
		--number-directions 1
		--sparse-dim 3
		--residual-threshold 1e-4
		--data ${CMAKE_CURRENT_SOURCE_DIR}/pre/data.m 
		--solution-first-factor first-factor_rank3_symmetric.m
		--solution-second-factor second-factor_rank3_symmetric.m)
set_tests_properties (MatrixFactorizer_Indefinite_SymmetricRank3
	PROPERTIES WILL_FAIL false)
set_tests_properties (MatrixFactorizer_Indefinite_SymmetricRank3
	PROPERTIES PASS_REGULAR_EXPRESSION "Norm of .*e-05")

# test 1b: checking results
add_test (MatrixFactorizer_Indefinite_SymmetricRank3_check_first
	${DIFFCMD} ${NDIFF_OPTIONS} ${CMAKE_CURRENT_SOURCE_DIR}/post/first-factor.m first-factor_rank3_symmetric.m)
set_tests_properties (MatrixFactorizer_Indefinite_SymmetricRank3_check_first
	PROPERTIES WILL_FAIL false)
set_tests_properties (MatrixFactorizer_Indefinite_SymmetricRank3_check_first
	PROPERTIES DEPENDS MatrixFactorizer_Indefinite_SymmetricRank3)

# test 1c: checking results
add_test (MatrixFactorizer_Indefinite_SymmetricRank3_check_second
	${DIFFCMD} ${NDIFF_OPTIONS} ${CMAKE_CURRENT_SOURCE_DIR}/post/second-factor.m second-factor_rank3_symmetric.m)
set_tests_properties (MatrixFactorizer_Indefinite_SymmetricRank3_check_second
	PROPERTIES WILL_FAIL false)
set_tests_properties (MatrixFactorizer_Indefinite_SymmetricRank3_check_second
	PROPERTIES DEPENDS MatrixFactorizer_Indefinite_SymmetricRank3)

# test 2: same as test 1 but using parallel code
if (MPI_FOUND)
	add_test (MatrixFactorizer_Indefinite_SymmetricRank3_mpi
		${MPIEXEC} -${MPIEXEC_NUMPROC_FLAG} 2
		${PROJECT_BINARY_DIR}/Examples/MatrixFactorizer/MatrixFactorizer 
			--type-space-x "lp" 
			--px 2 
			--type-space-y "lp" 
			--py 2 
			--powery 2 
			--delta 1e-6 
			--projection-delta 1e-6
			--maxiter 50
			--max-loops 11
			--orthogonal-directions 1
			--number-directions 1
			--sparse-dim 3
			--residual-threshold 1e-4
			--data ${CMAKE_CURRENT_SOURCE_DIR}/pre/data.m 
			--solution-first-factor first-factor_mpi_rank3_symmetric.m
			--solution-second-factor second-factor_mpi_rank3_symmetric.m)
	set_tests_properties (MatrixFactorizer_Indefinite_SymmetricRank3_mpi
		PROPERTIES WILL_FAIL false)
	set_tests_properties (MatrixFactorizer_Indefinite_SymmetricRank3_mpi
		PROPERTIES PASS_REGULAR_EXPRESSION "Norm of .*e-05")
	
	# test 1b: checking results
	add_test (MatrixFactorizer_Indefinite_SymmetricRank3_mpi_check_first
		${DIFFCMD} ${NDIFF_OPTIONS} ${CMAKE_CURRENT_SOURCE_DIR}/post/first-factor.m first-factor_mpi_rank3_symmetric.m)
	set_tests_properties (MatrixFactorizer_Indefinite_SymmetricRank3_mpi_check_first
		PROPERTIES WILL_FAIL false)
	set_tests_properties (MatrixFactorizer_Indefinite_SymmetricRank3_mpi_check_first
		PROPERTIES DEPENDS MatrixFactorizer_Indefinite_SymmetricRank3_mpi)
	
	# test 1c: checking results
	add_test (MatrixFactorizer_Indefinite_SymmetricRank3_mpi_check_second
		${DIFFCMD} ${NDIFF_OPTIONS} ${CMAKE_CURRENT_SOURCE_DIR}/post/second-factor.m second-factor_mpi_rank3_symmetric.m)
	set_tests_properties (MatrixFactorizer_Indefinite_SymmetricRank3_mpi_check_second
		PROPERTIES WILL_FAIL false)
	set_tests_properties (MatrixFactorizer_Indefinite_SymmetricRank3_mpi_check_second
		PROPERTIES DEPENDS MatrixFactorizer_Indefinite_SymmetricRank3_mpi)
endif (MPI_FOUND)
