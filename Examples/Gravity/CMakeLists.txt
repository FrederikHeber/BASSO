# add the binary tree to the search path for
# include files
include_directories (${CMAKE_CURRENT_BINARY_DIR})

# add the executable for matrix factorization
add_executable (Gravity 
	Gravity.cpp
	Options/GravityOptions.cpp)

install (FILES Options/GravityOptions.hpp DESTINATION include)

target_link_libraries ( Gravity
	${EXTRA_LIBS}
	${NLOPT_LIBRARIES}
	${GSL_LIBRARY} 
	${CMAKE_THREAD_LIBS_INIT} 
	${Boost_LIBRARIES}
	${POCO_SQL_LIBRARIES}
	${LIBRT_LIBRARIES}
	${MPI_LIBRARIES})
install (TARGETS Gravity DESTINATION bin)
