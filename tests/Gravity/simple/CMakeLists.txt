# does simple test call work?
add_test (GravitySimple 
	${PROJECT_BINARY_DIR}/Examples/Gravity/Gravity 
		--type-space-x "lp" 
		--px 2 
		--type-space-y "lp" 
		--py 2 
		--powery 2 
		--maxiter 10
		--delta 0.01 
		--depth 0.5
		--discretization 10
		--density-file density.m)
set_tests_properties (GravitySimple
	PROPERTIES PASS_REGULAR_EXPRESSION "olution after 1 with relative residual of")
