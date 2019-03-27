# add the binary tree to the search path for
# include files
include_directories (${CMAKE_CURRENT_BINARY_DIR})

# add the executable for matrix factorization
add_library (GravityBase 
	Options/GravityOptions.cpp)

install (FILES Options/GravityOptions.hpp DESTINATION include)

# add the executable
add_executable (Gravity
	Gravity.cpp)
target_link_libraries ( Gravity
	GravityBase
	${EXTRA_LIBS}
	${NLOPT_LIBRARIES}
	${GSL_LIBRARY} 
	${CMAKE_THREAD_LIBS_INIT} 
	${BOOST_CORE_LIBRARIES}
	${POCO_SQL_LIBRARIES}
	${LIBRT_LIBRARIES}
	${MPI_LIBRARIES})
install (TARGETS Gravity DESTINATION bin)

add_executable (GravityConfigurator 
	GravityConfigurator.cpp)
target_link_libraries ( GravityConfigurator 
	GravityBase
	${EXTRA_LIBS} 
	${NLOPT_LIBRARIES}
	${GSL_LIBRARY} 
	${CMAKE_THREAD_LIBS_INIT} 
	${BOOST_CORE_LIBRARIES}
	${POCO_SQL_LIBRARIES}
	${LIBRT_LIBRARIES}) 
install (TARGETS GravityConfigurator DESTINATION bin)
