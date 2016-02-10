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