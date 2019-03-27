# add the binary tree to the search path for
# include files
include_directories (${CMAKE_CURRENT_BINARY_DIR})

add_library (MatrixFactorizerBase
	Database/AccumulatedValues.cpp
	Database/InnerProblemDatabase.cpp
	Database/IterationInformation.cpp
	Database/TableDataAccumulator.cpp
	Database/TableDataOperator.cpp
	Helpers/detail.cpp
	Options/MatrixFactorizerOptions.cpp
	Solvers/InRangeSolver.cpp
	Solvers/MatrixFactorization.cpp
	Work/Master.cpp
	Work/Slave.cpp)

install (FILES Database/AccumulatedValues.hpp DESTINATION include)
install (FILES Database/InnerProblemDatabase.hpp DESTINATION include)
install (FILES Database/IterationInformation.hpp DESTINATION include)
install (FILES Database/MaxMinAverage.hpp DESTINATION include)
install (FILES Database/TableDataAccumulator.hpp DESTINATION include)
install (FILES Database/TableDataOperator.hpp DESTINATION include)
install (FILES Helpers/detail.hpp DESTINATION include)
install (FILES Options/MatrixFactorizerOptions.hpp DESTINATION include)
install (FILES ScalingAmbiguityRemover/MatrixProductEqualizer.hpp DESTINATION include)
install (FILES ScalingAmbiguityRemover/MatrixProductRenormalizer.hpp DESTINATION include)
install (FILES ScalingAmbiguityRemover/ScalingAmbiguityMaintainer.hpp DESTINATION include)
install (FILES ScalingAmbiguityRemover/ScalingAmbiguityRemover.hpp DESTINATION include)
install (FILES ScalingAmbiguityRemover/TraceRenormalizer.hpp DESTINATION include)
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
	${BOOST_CORE_LIBRARIES}
	${POCO_SQL_LIBRARIES}
	${LIBRT_LIBRARIES}
	${MPI_LIBRARIES})
# MatrixFactorizer is parallelized with OpenMP
if (USE_OPENMP)
	set_target_properties(MatrixFactorizerBase PROPERTIES
		COMPILE_FLAGS "${OpenMP_CXX_FLAGS}")
	set_target_properties(MatrixFactorizer PROPERTIES
		LINK_FLAGS "${OpenMP_CXX_FLAGS}")
endif  (USE_OPENMP)
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
	${BOOST_CORE_LIBRARIES}
	${POCO_SQL_LIBRARIES}
	${LIBRT_LIBRARIES}) 
install (TARGETS MatrixFactorizerConfigurator DESTINATION bin)
