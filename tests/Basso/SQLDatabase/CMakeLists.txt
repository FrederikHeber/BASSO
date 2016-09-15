# does simple test call work with writing sqlite file?
add_test (Basso_SQLDatabase
	${PROJECT_BINARY_DIR}/Examples/Basso/Basso 
		--iteration-file sqltest.db
		--type-space-x "lp" 
		--px 2 
		--type-space-y "lp" 
		--py 2 
		--powery 2 
		--delta 0.01 
		--matrix ${CMAKE_CURRENT_SOURCE_DIR}/pre/matrix.m 
		--rhs ${CMAKE_CURRENT_SOURCE_DIR}/pre/vector.m)
set_tests_properties (Basso_SQLDatabase
	PROPERTIES WILL_FAIL false)

if (SQLITE_FOUND)
	add_test (Basso_SQLDatabase_sqliteextract_1
		${SQLITECMD} 
			-csv sqltest.db 
			"SELECT iterations, residual FROM overall;")
	set_tests_properties (Basso_SQLDatabase_sqliteextract_1
		PROPERTIES DEPENDS Basso_SQLDatabase)
	set_tests_properties (Basso_SQLDatabase_sqliteextract_1
		PROPERTIES PASS_REGULAR_EXPRESSION "1,1.02.*e-08")

	add_test (Basso_SQLDatabase_sqliteextract_2
		${SQLITECMD} 
			-csv sqltest.db 
			"SELECT p,r residual FROM parameters;")
	set_tests_properties (Basso_SQLDatabase_sqliteextract_2
		PROPERTIES DEPENDS Basso_SQLDatabase)
	set_tests_properties (Basso_SQLDatabase_sqliteextract_2
		PROPERTIES PASS_REGULAR_EXPRESSION "2.0,2.0")
endif (SQLITE_FOUND)

# do we deactivate per_iteration with every-nth-tuple equal 0
add_test (Basso_EveryNthTuple_deactivated
	${PROJECT_BINARY_DIR}/Examples/Basso/Basso
		--verbose 1
		--iteration-file sqltest_deactivated.db
		--algorithm Landweber
		--every-nth-tuple 0
		--type-space-x "lp" 
		--px 2 
		--type-space-y "lp" 
		--py 2 
		--powery 2 
		--delta 0.01 
		--matrix ${CMAKE_CURRENT_SOURCE_DIR}/pre/matrix.m 
		--rhs ${CMAKE_CURRENT_SOURCE_DIR}/pre/vector.m)
set_tests_properties (Basso_EveryNthTuple_deactivated
		PROPERTIES PASS_REGULAR_EXPRESSION "Per iteration tuples are not added")

if (SQLITE_FOUND)
	add_test (Basso_EveryNthTuple_deactivated_sqliteextract_1
		${SQLITECMD} 
			-csv sqltest_deactivated.db 
			"SELECT iteration, residual FROM per_iteration;")
	set_tests_properties (Basso_EveryNthTuple_deactivated_sqliteextract_1
		PROPERTIES WILL_FAIL true)
	set_tests_properties (Basso_EveryNthTuple_deactivated_sqliteextract_1
		PROPERTIES DEPENDS Basso_EveryNthTuple_deactivated)
endif (SQLITE_FOUND)

# do we activate per_iteration with every-nth-tuple unequal 0
add_test (Basso_EveryNthTuple_activated
	${PROJECT_BINARY_DIR}/Examples/Basso/Basso
		--iteration-file sqltest_activated.db
		--algorithm Landweber
		--every-nth-tuple 1
		--type-space-x "lp" 
		--px 2 
		--type-space-y "lp" 
		--py 2 
		--powery 2 
		--delta 0.01 
		--matrix ${CMAKE_CURRENT_SOURCE_DIR}/pre/matrix.m 
		--rhs ${CMAKE_CURRENT_SOURCE_DIR}/pre/vector.m)
set_tests_properties (Basso_EveryNthTuple_activated
	PROPERTIES WILL_FAIL false)

if (SQLITE_FOUND)
	add_test (Basso_EveryNthTuple_activated_sqliteextract_1
		${SQLITECMD} 
			-csv sqltest_activated.db 
			"SELECT iteration, residual FROM per_iteration;")
	set_tests_properties (Basso_EveryNthTuple_activated_sqliteextract_1
		PROPERTIES WILL_FAIL false)
	set_tests_properties (Basso_EveryNthTuple_activated_sqliteextract_1
		PROPERTIES DEPENDS Basso_EveryNthTuple_activated)
endif (SQLITE_FOUND)

# do we always add last tuple?
add_test (Basso_EveryNthTuple_lasttuple
	${PROJECT_BINARY_DIR}/Examples/Basso/Basso
		--iteration-file sqltest_lasttuple.db
		--algorithm Landweber
		--every-nth-tuple 10
		--type-space-x "lp" 
		--px 2 
		--type-space-y "lp" 
		--py 2 
		--powery 2 
		--delta 0.01 
		--matrix ${CMAKE_CURRENT_SOURCE_DIR}/pre/matrix.m 
		--rhs ${CMAKE_CURRENT_SOURCE_DIR}/pre/vector.m)
set_tests_properties (Basso_EveryNthTuple_lasttuple
	PROPERTIES WILL_FAIL false)

if (SQLITE_FOUND)
	add_test (Basso_EveryNthTuple_lasttuple_sqliteextract_1
		${SQLITECMD} 
			-csv sqltest_lasttuple.db 
			"SELECT iteration, residual FROM per_iteration;")
	set_tests_properties (Basso_EveryNthTuple_lasttuple_sqliteextract_1
		PROPERTIES PASS_REGULAR_EXPRESSION "1,1.02.*e-08")
	set_tests_properties (Basso_EveryNthTuple_lasttuple_sqliteextract_1
		PROPERTIES DEPENDS Basso_EveryNthTuple_lasttuple)
endif (SQLITE_FOUND)
