# test 1: Matrix factorization with rank 1 
add_test (MatrixFactorizer_SQLDatabase 
	${PROJECT_BINARY_DIR}/Examples/MatrixFactorizer/MatrixFactorizer 
		--iteration-file sqltest.db
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
		--data ${CMAKE_CURRENT_SOURCE_DIR}/pre/data.m 
		--solution-first-factor first-factor.m
		--solution-second-factor second-factor.m)
set_tests_properties (MatrixFactorizer_SQLDatabase
	PROPERTIES WILL_FAIL false)

if (SQLITE_FOUND)
	add_test (MatrixFactorizer_SQLDatabase_sqliteextract_1
		${SQLITECMD} 
			-csv sqltest.db 
			"SELECT loops, residual FROM loop_overall;")
	set_tests_properties (MatrixFactorizer_SQLDatabase_sqliteextract_1
		PROPERTIES DEPENDS MatrixFactorizer_SQLDatabase)
	set_tests_properties (MatrixFactorizer_SQLDatabase_sqliteextract_1
		PROPERTIES PASS_REGULAR_EXPRESSION "1,1.69.*e-08")

	add_test (MatrixFactorizer_SQLDatabase_sqliteextract_2
		${SQLITECMD} 
			-csv sqltest.db 
			"SELECT p,r,k,n residual FROM parameters;")
	set_tests_properties (MatrixFactorizer_SQLDatabase_sqliteextract_2
		PROPERTIES DEPENDS MatrixFactorizer_SQLDatabase)
	set_tests_properties (MatrixFactorizer_SQLDatabase_sqliteextract_2
		PROPERTIES PASS_REGULAR_EXPRESSION "2.0,2.0,4,4")
endif (SQLITE_FOUND)
