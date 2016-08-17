# This is based on MatrixFactorizer_Indefinite_SymmetricRank3
add_test (MatrixFactorizer_ProjectionDelta_HighLow
	${PROJECT_BINARY_DIR}/Examples/MatrixFactorizer/MatrixFactorizer 
		--type-space-x "lp" 
		--px 2 
		--type-space-y "lp" 
		--py 2 
		--powery 2 
		--delta 1e-2 
		--projection-delta 1e-6
		--maxiter 50
		--max-loops 11
		--orthogonal-directions 1
		--number-directions 1
		--sparse-dim 3
		--residual-threshold 1e-4
		--data ${CMAKE_CURRENT_SOURCE_DIR}/pre/data.m 
		--solution-first-factor first-factor_rank3_symmetric_HighLow.m
		--solution-second-factor second-factor_rank3_symmetric_HighLow.m)
set_tests_properties (MatrixFactorizer_ProjectionDelta_HighLow
	PROPERTIES WILL_FAIL false)
set_tests_properties (MatrixFactorizer_ProjectionDelta_HighLow
	PROPERTIES 
		PASS_REGULAR_EXPRESSION "Iteration was stopped by RelativeChangeResiduum"
		PASS_REGULAR_EXPRESSION "Loop iteration was performed 2 times"
		PASS_REGULAR_EXPRESSION "Norm of .*0.01"
		)

# This is based on MatrixFactorizer_Indefinite_SymmetricRank3
add_test (MatrixFactorizer_ProjectionDelta_LowHigh
	${PROJECT_BINARY_DIR}/Examples/MatrixFactorizer/MatrixFactorizer 
		--type-space-x "lp" 
		--px 2 
		--type-space-y "lp" 
		--py 2 
		--powery 2 
		--delta 1e-6
		--projection-delta 1e-2
		--maxiter 50
		--max-loops 11
		--orthogonal-directions 1
		--number-directions 1
		--sparse-dim 3
		--residual-threshold 1e-4
		--data ${CMAKE_CURRENT_SOURCE_DIR}/pre/data.m 
		--solution-first-factor first-factor_rank3_symmetric_LowHigh.m
		--solution-second-factor second-factor_rank3_symmetric_LowHigh.m)
set_tests_properties (MatrixFactorizer_ProjectionDelta_LowHigh
	PROPERTIES WILL_FAIL false)
set_tests_properties (MatrixFactorizer_ProjectionDelta_LowHigh
	PROPERTIES 
		PASS_REGULAR_EXPRESSION "Iteration was stopped by RelativeChangeResiduum"
		PASS_REGULAR_EXPRESSION "Loop iteration was performed 9 times"
		PASS_REGULAR_EXPRESSION "Norm of .*0.01"
		)

# test 2: same as test 1 but using parallel code
if (MPI_FOUND)
	add_test (MatrixFactorizer_ProjectionDelta_HighLow_mpi
		${MPIEXEC} -${MPIEXEC_NUMPROC_FLAG} 2
		${PROJECT_BINARY_DIR}/Examples/MatrixFactorizer/MatrixFactorizer 
			--type-space-x "lp" 
			--px 2 
			--type-space-y "lp" 
			--py 2 
			--powery 2 
			--delta 1e-2 
			--projection-delta 1e-6
			--maxiter 50
			--max-loops 11
			--orthogonal-directions 1
			--number-directions 1
			--sparse-dim 3
			--residual-threshold 1e-4
			--data ${CMAKE_CURRENT_SOURCE_DIR}/pre/data.m 
			--solution-first-factor first-factor_mpi_rank3_symmetric_HighLow_mpi.m
			--solution-second-factor second-factor_mpi_rank3_symmetric_HighLow_mpi.m)
	set_tests_properties (MatrixFactorizer_ProjectionDelta_HighLow_mpi
		PROPERTIES WILL_FAIL false)
	set_tests_properties (MatrixFactorizer_ProjectionDelta_HighLow_mpi
		PROPERTIES 
			PASS_REGULAR_EXPRESSION "Iteration was stopped by RelativeChangeResiduum"
			PASS_REGULAR_EXPRESSION "Loop iteration was performed 2 times"
		PASS_REGULAR_EXPRESSION "Norm of .*0.01"
			)

	add_test (MatrixFactorizer_ProjectionDelta_LowHigh_mpi
		${MPIEXEC} -${MPIEXEC_NUMPROC_FLAG} 2
		${PROJECT_BINARY_DIR}/Examples/MatrixFactorizer/MatrixFactorizer 
			--type-space-x "lp" 
			--px 2 
			--type-space-y "lp" 
			--py 2 
			--powery 2 
			--delta 1e-6
			--projection-delta 1e-2
			--maxiter 50
			--max-loops 11
			--orthogonal-directions 1
			--number-directions 1
			--sparse-dim 3
			--residual-threshold 1e-4
			--data ${CMAKE_CURRENT_SOURCE_DIR}/pre/data.m 
			--solution-first-factor first-factor_mpi_rank3_symmetric_LowHigh_mpi.m
			--solution-second-factor second-factor_mpi_rank3_symmetric_LowHigh_mpi.m)
	set_tests_properties (MatrixFactorizer_ProjectionDelta_LowHigh_mpi
		PROPERTIES WILL_FAIL false)
	set_tests_properties (MatrixFactorizer_ProjectionDelta_LowHigh_mpi
	PROPERTIES 
			PASS_REGULAR_EXPRESSION "Iteration was stopped by RelativeChangeResiduum"
			PASS_REGULAR_EXPRESSION "Loop iteration was performed 9 times"
		PASS_REGULAR_EXPRESSION "Norm of .*0.01"
			)
endif (MPI_FOUND)
