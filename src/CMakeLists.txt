# add the executable for matrix factorization
add_executable (MatrixFactorizer MatrixFactorizer.cpp)
target_link_libraries ( MatrixFactorizer 
	${EXTRA_LIBS}
	${NLOPT_LIBRARIES}
	${GSL_LIBRARY} 
	${CMAKE_THREAD_LIBS_INIT} 
	${Boost_LIBRARIES}
	${POCO_SQL_LIBRARIES}
	${LIBRT_LIBRARIES}
	${MPI_LIBRARIES})
# MatrixFactorizer is parallelized with MPI
install (TARGETS MatrixFactorizer DESTINATION bin)

add_executable (MatrixFactorizerConfigurator MatrixFactorizerConfigurator.cpp)
target_link_libraries ( MatrixFactorizerConfigurator 
	${EXTRA_LIBS} 
	${NLOPT_LIBRARIES}
	${GSL_LIBRARY} 
	${CMAKE_THREAD_LIBS_INIT} 
	${Boost_LIBRARIES}
	${POCO_SQL_LIBRARIES}
	${LIBRT_LIBRARIES}) 
install (TARGETS MatrixFactorizerConfigurator DESTINATION bin)
