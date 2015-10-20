# test 1: updateToMostParallel
add_test (Basso_SESOP_RandomMatrix_UpdateToMostParallel
	${PROJECT_BINARY_DIR}/Examples/Basso/Basso 
		--algorithm SESOP
		--number-directions 2
		--update-algorithm 1
		--enforceRandomMapping 0
		--type-space-x "lp" 
		--px 2 
		--type-space-y "lp" 
		--py 0 
		--powerx 2
		--powery 2
		--delta 0.01 
		--maxiter 100 
		--iteration-file SESOP_UpdateToMostParallel.db
		--matrix ${PROJECT_SOURCE_DIR}/tests/Basso/SESOP/VariousSearchspaceUpdates/pre/randommatrix.m 
		--rhs ${PROJECT_SOURCE_DIR}/tests/Basso/SESOP/VariousSearchspaceUpdates/pre/randomvector.m
		--solution solution_mostparallel.m)
set_tests_properties (Basso_SESOP_RandomMatrix_UpdateToMostParallel
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 19 with relative residual of 0.0089.*")

# test 1b: checking results
add_test (Basso_SESOP_RandomMatrix_UpdateToMostParallel_check
	diff ${PROJECT_SOURCE_DIR}/tests/Basso/SESOP/VariousSearchspaceUpdates/post/solution_mostparallel.m  solution_mostparallel.m)
set_tests_properties (Basso_SESOP_RandomMatrix_UpdateToMostParallel_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (Basso_SESOP_RandomMatrix_UpdateToMostParallel_check
	PROPERTIES DEPENDS Basso_SESOP_RandomMatrix_UpdateToMostParallel)

# test 2: updateToMostOrthogonal
add_test (Basso_SESOP_RandomMatrix_UpdateToMostOrthogonal
	${PROJECT_BINARY_DIR}/Examples/Basso/Basso 
		--algorithm SESOP
		--number-directions 2
		--update-algorithm 2
		--enforceRandomMapping 0
		--type-space-x "lp" 
		--px 2 
		--type-space-y "lp" 
		--py 0 
		--powerx 2
		--powery 2
		--delta 0.01 
		--maxiter 100 
		--iteration-file SESOP_UpdateToMostOrthogonal.db
		--matrix ${PROJECT_SOURCE_DIR}/tests/Basso/SESOP/VariousSearchspaceUpdates/pre/randommatrix.m 
		--rhs ${PROJECT_SOURCE_DIR}/tests/Basso/SESOP/VariousSearchspaceUpdates/pre/randomvector.m
		--solution solution_mostorthogonal.m)
set_tests_properties (Basso_SESOP_RandomMatrix_UpdateToMostOrthogonal
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 25 with relative residual of 0.0047.*")

# test 2b: checking results
add_test (Basso_SESOP_RandomMatrix_UpdateToMostOrthogonal_check
	diff ${PROJECT_SOURCE_DIR}/tests/Basso/SESOP/VariousSearchspaceUpdates/post/solution_mostorthogonal.m  solution_mostorthogonal.m)
set_tests_properties (Basso_SESOP_RandomMatrix_UpdateToMostOrthogonal_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (Basso_SESOP_RandomMatrix_UpdateToMostOrthogonal_check
	PROPERTIES DEPENDS Basso_SESOP_RandomMatrix_UpdateToMostOrthogonal)

# test 3: updateToMostParallel_randomMapping
add_test (Basso_SESOP_RandomMatrix_UpdateToMostParallel_randomMapping
	${PROJECT_BINARY_DIR}/Examples/Basso/Basso 
		--algorithm SESOP
		--number-directions 2
		--update-algorithm 1
		--enforceRandomMapping 1
		--type-space-x "lp" 
		--px 2 
		--type-space-y "lp" 
		--py 0 
		--powerx 2
		--powery 2
		--delta 0.01 
		--maxiter 100 
		--iteration-file SESOP_UpdateToMostParallel_randomMapping.db
		--matrix ${PROJECT_SOURCE_DIR}/tests/Basso/SESOP/VariousSearchspaceUpdates/pre/randommatrix.m 
		--rhs ${PROJECT_SOURCE_DIR}/tests/Basso/SESOP/VariousSearchspaceUpdates/pre/randomvector.m
		--solution solution_mostparallel_randommapping.m)
set_tests_properties (Basso_SESOP_RandomMatrix_UpdateToMostParallel_randomMapping
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 21 with relative residual of 0.0083.*")

# test 3b: checking results
add_test (Basso_SESOP_RandomMatrix_UpdateToMostParallel_randomMapping_check
	diff ${PROJECT_SOURCE_DIR}/tests/Basso/SESOP/VariousSearchspaceUpdates/post/solution_mostparallel_randommapping.m  solution_mostparallel_randommapping.m)
set_tests_properties (Basso_SESOP_RandomMatrix_UpdateToMostParallel_randomMapping_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (Basso_SESOP_RandomMatrix_UpdateToMostParallel_randomMapping_check
	PROPERTIES DEPENDS Basso_SESOP_RandomMatrix_UpdateToMostParallel_randomMapping)

# test 4: updateToMostOrthogonal_randomMapping
add_test (Basso_SESOP_RandomMatrix_UpdateToMostOrthogonal_randomMapping
	${PROJECT_BINARY_DIR}/Examples/Basso/Basso 
		--algorithm SESOP
		--number-directions 2
		--update-algorithm 2
		--enforceRandomMapping 1
		--type-space-x "lp" 
		--px 2 
		--type-space-y "lp" 
		--py 0 
		--powerx 2
		--powery 2
		--delta 0.01 
		--maxiter 100 
		--iteration-file SESOP_UpdateToMostOrthogonal_randomMapping.db
		--matrix ${PROJECT_SOURCE_DIR}/tests/Basso/SESOP/VariousSearchspaceUpdates/pre/randommatrix.m 
		--rhs ${PROJECT_SOURCE_DIR}/tests/Basso/SESOP/VariousSearchspaceUpdates/pre/randomvector.m
		--solution solution_mostorthogonal_randommapping.m)
set_tests_properties (Basso_SESOP_RandomMatrix_UpdateToMostOrthogonal_randomMapping
	PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 23 with relative residual of 0.0094.*")

# test 4b: checking results
add_test (Basso_SESOP_RandomMatrix_UpdateToMostOrthogonal_randomMapping_check
	diff ${PROJECT_SOURCE_DIR}/tests/Basso/SESOP/VariousSearchspaceUpdates/post/solution_mostorthogonal_randommapping.m  solution_mostorthogonal_randommapping.m)
set_tests_properties (Basso_SESOP_RandomMatrix_UpdateToMostOrthogonal_randomMapping_check
	PROPERTIES WILL_FAIL false)
set_tests_properties (Basso_SESOP_RandomMatrix_UpdateToMostOrthogonal_randomMapping_check
	PROPERTIES DEPENDS Basso_SESOP_RandomMatrix_UpdateToMostOrthogonal_randomMapping)
