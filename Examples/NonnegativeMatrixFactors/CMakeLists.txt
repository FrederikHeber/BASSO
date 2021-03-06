# add the binary tree to the search path for
# include files
include_directories (${CMAKE_CURRENT_BINARY_DIR})

# add the executable for matrix factorization
add_library (NonnegativeMatrixFactorsBase 
	Options/NonnegativeMatrixFactorsOptions.cpp
	TruncatedSVD.cpp)

install (FILES Options/NonnegativeMatrixFactorsOptions.hpp DESTINATION include)
install (FILES TruncatedSVD.hpp DESTINATION include)

# add the executable
add_executable (NonnegativeMatrixFactors
	NonnegativeMatrixFactors.cpp)
target_link_libraries ( NonnegativeMatrixFactors
	NonnegativeMatrixFactorsBase
	${EXTRA_LIBS}
	${ARPACK_LIBRARIES}
	${NLOPT_LIBRARIES}
	${GSL_LIBRARY} 
	${CMAKE_THREAD_LIBS_INIT} 
	${BOOST_CORE_LIBRARIES}
	${POCO_SQL_LIBRARIES}
	${LIBRT_LIBRARIES}
	${MPI_LIBRARIES})
install (TARGETS NonnegativeMatrixFactors DESTINATION bin)

add_executable (NonnegativeMatrixFactorsConfigurator 
	NonnegativeMatrixFactorsConfigurator.cpp)
target_link_libraries ( NonnegativeMatrixFactorsConfigurator 
	NonnegativeMatrixFactorsBase
	${EXTRA_LIBS} 
	${NLOPT_LIBRARIES}
	${GSL_LIBRARY} 
	${CMAKE_THREAD_LIBS_INIT} 
	${BOOST_CORE_LIBRARIES}
	${POCO_SQL_LIBRARIES}
	${LIBRT_LIBRARIES}) 
install (TARGETS NonnegativeMatrixFactorsConfigurator DESTINATION bin)
