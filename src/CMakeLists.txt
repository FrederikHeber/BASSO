# add the executable
add_executable (Basso Basso.cpp)
target_link_libraries ( Basso 
	${EXTRA_LIBS} 
	${NLOPT_LIBRARIES}
	${GSL_LIBRARY} 
	${CMAKE_THREAD_LIBS_INIT} 
	${Boost_LIBRARIES}
	${POCO_SQL_LIBRARIES}
	${LIBRT_LIBRARIES}) 

add_executable (BassoConfigurator BassoConfigurator.cpp)
target_link_libraries ( BassoConfigurator 
	${EXTRA_LIBS} 
	${NLOPT_LIBRARIES}
	${GSL_LIBRARY} 
	${CMAKE_THREAD_LIBS_INIT} 
	${Boost_LIBRARIES}
	${POCO_SQL_LIBRARIES}
	${LIBRT_LIBRARIES}) 

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

add_executable (MatrixFactorizerConfigurator MatrixFactorizerConfigurator.cpp)
target_link_libraries ( MatrixFactorizerConfigurator 
	${EXTRA_LIBS} 
	${NLOPT_LIBRARIES}
	${GSL_LIBRARY} 
	${CMAKE_THREAD_LIBS_INIT} 
	${Boost_LIBRARIES}
	${POCO_SQL_LIBRARIES}
	${LIBRT_LIBRARIES}) 

# add the install targets
install (TARGETS Basso DESTINATION bin)
