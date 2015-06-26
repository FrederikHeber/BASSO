# add the executable for matrix factorization
add_executable (MatrixFactorizer 
	MatrixFactorizerBase/MatrixFactorizer.cpp)
target_link_libraries ( MatrixFactorizer 
	${EXTRA_LIBS}
	${NLOPT_LIBRARIES}
	${GSL_LIBRARY} 
	${CMAKE_THREAD_LIBS_INIT} 
	${Boost_LIBRARIES}
	${POCO_SQL_LIBRARIES}
	${LIBRT_LIBRARIES}
	${MPI_LIBRARIES})
install (TARGETS MatrixFactorizer DESTINATION bin)
install (FILES MatrixFactorizerBase/Work/WorkPackage.hpp DESTINATION include)
install (FILES MatrixFactorizerBase/Work/WorkResult.hpp DESTINATION include)

add_executable (MatrixFactorizerConfigurator 
	MatrixFactorizerBase/MatrixFactorizerConfigurator.cpp)
target_link_libraries ( MatrixFactorizerConfigurator 
	${EXTRA_LIBS} 
	${NLOPT_LIBRARIES}
	${GSL_LIBRARY} 
	${CMAKE_THREAD_LIBS_INIT} 
	${Boost_LIBRARIES}
	${POCO_SQL_LIBRARIES}
	${LIBRT_LIBRARIES}) 
install (TARGETS MatrixFactorizerConfigurator DESTINATION bin)
