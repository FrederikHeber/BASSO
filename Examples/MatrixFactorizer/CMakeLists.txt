# add the binary tree to the search path for
# include files
include_directories (${CMAKE_CURRENT_BINARY_DIR})

add_library (MatrixFactorizerBase
	Database/IterationInformation.cpp
	Helpers/detail.cpp
	Options/MatrixFactorizerOptions.cpp
	Solvers/InRangeSolver.cpp
	Solvers/MatrixFactorization.cpp
	Work/Master.cpp
	Work/Slave.cpp)

install (FILES Database/IterationInformation.hpp DESTINATION include)
install (FILES Helpers/detail.hpp DESTINATION include)
install (FILES IterationChecks/MaxIterationsChecker.hpp DESTINATION include)
install (FILES IterationChecks/RelativeResidualChecker.hpp DESTINATION include)
install (FILES IterationChecks/ResidualChecker.hpp DESTINATION include)
install (FILES Options/filesystem_path_serialization.hpp DESTINATION include)
install (FILES Options/MatrixFactorizerOptions.hpp DESTINATION include)
install (FILES ScalingAmbigutiyRemover/MatrixProductEqualizer.hpp DESTINATION include)
install (FILES ScalingAmbigutiyRemover/MatrixProductRenormalizer.hpp DESTINATION include)
install (FILES ScalingAmbigutiyRemover/TraceRenormalizer.hpp DESTINATION include)
install (FILES Solvers/InRangeSolver.hpp DESTINATION include)
install (FILES Work/Master.hpp DESTINATION include)
install (FILES Work/Slave.hpp DESTINATION include)
install (FILES Work/WorkPackage.hpp DESTINATION include)
install (FILES Work/WorkResult.hpp DESTINATION include)

# add the executable for matrix factorization
add_executable (MatrixFactorizer 
	MatrixFactorizer.cpp)
target_link_libraries ( MatrixFactorizer
	MatrixFactorizerBase
	RangeProjectorBase
	${EXTRA_LIBS}
	${NLOPT_LIBRARIES}
	${GSL_LIBRARY} 
	${CMAKE_THREAD_LIBS_INIT} 
	${Boost_LIBRARIES}
	${POCO_SQL_LIBRARIES}
	${LIBRT_LIBRARIES}
	${MPI_LIBRARIES})
install (TARGETS MatrixFactorizer DESTINATION bin)

add_executable (MatrixFactorizerConfigurator 
	MatrixFactorizerConfigurator.cpp)
target_link_libraries ( MatrixFactorizerConfigurator
	MatrixFactorizerBase
	RangeProjectorBase
	${EXTRA_LIBS} 
	${NLOPT_LIBRARIES}
	${GSL_LIBRARY} 
	${CMAKE_THREAD_LIBS_INIT} 
	${Boost_LIBRARIES}
	${POCO_SQL_LIBRARIES}
	${LIBRT_LIBRARIES}) 
install (TARGETS MatrixFactorizerConfigurator DESTINATION bin)