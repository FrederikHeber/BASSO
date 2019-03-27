# add the binary tree to the search path for
# include files
include_directories (${CMAKE_CURRENT_BINARY_DIR})

# add the executables
add_executable (random_matrix_inverseproblem_1
	../../doc/listings/random_matrix_inverseproblem_1.cpp)
target_link_libraries ( random_matrix_inverseproblem_1
	${EXTRA_LIBS} 
	${NLOPT_LIBRARIES}
	${GSL_LIBRARY} 
	${CMAKE_THREAD_LIBS_INIT} 
	${BOOST_CORE_LIBRARIES}
	${POCO_SQL_LIBRARIES}
	${LIBRT_LIBRARIES}) 

add_test (random_matrix_inverseproblem_1_test
        random_matrix_inverseproblem_1)
set_tests_properties (random_matrix_inverseproblem_1_test
        PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 50 with relative residual of 0.002.*")

add_executable (random_matrix_inverseproblem_2
	../../doc/listings/random_matrix_inverseproblem_2.cpp)
target_link_libraries ( random_matrix_inverseproblem_2
	${EXTRA_LIBS} 
	${NLOPT_LIBRARIES}
	${GSL_LIBRARY} 
	${CMAKE_THREAD_LIBS_INIT} 
	${BOOST_CORE_LIBRARIES}
	${POCO_SQL_LIBRARIES}
	${LIBRT_LIBRARIES}) 

add_test (random_matrix_inverseproblem_2_test
        random_matrix_inverseproblem_2)

set_tests_properties (random_matrix_inverseproblem_2_test
        PROPERTIES PASS_REGULAR_EXPRESSION "Solution after 50 with relative residual of 0.002.*")


