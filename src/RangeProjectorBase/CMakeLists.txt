# add the binary tree to the search path for
# include files
include_directories (${CMAKE_CURRENT_BINARY_DIR})

add_library (RangeProjectorBase
	RangeProjector/RangeProjector.cpp
	Options/RangeProjectorOptions.cpp)

install (FILES RangeProjector/RangeProjector.hpp DESTINATION include)
install (FILES Options/RangeProjectorOptions.hpp DESTINATION include)

# add the executable
add_executable (RangeProjector
	RangeProjector.cpp)
target_link_libraries ( RangeProjector
	RangeProjectorBase
	${EXTRA_LIBS} 
	${NLOPT_LIBRARIES}
	${GSL_LIBRARY} 
	${CMAKE_THREAD_LIBS_INIT} 
	${Boost_LIBRARIES}
	${POCO_SQL_LIBRARIES}
	${LIBRT_LIBRARIES}) 
install (TARGETS RangeProjector DESTINATION bin)

add_executable (RangeProjectorConfigurator 
	RangeProjectorConfigurator.cpp)
target_link_libraries ( RangeProjectorConfigurator 
	RangeProjectorBase
	${EXTRA_LIBS} 
	${NLOPT_LIBRARIES}
	${GSL_LIBRARY} 
	${CMAKE_THREAD_LIBS_INIT} 
	${Boost_LIBRARIES}
	${POCO_SQL_LIBRARIES}
	${LIBRT_LIBRARIES}) 
install (TARGETS RangeProjectorConfigurator DESTINATION bin)
