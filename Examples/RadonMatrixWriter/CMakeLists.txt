# add the binary tree to the search path for
# include files
include_directories (${CMAKE_CURRENT_BINARY_DIR})

add_library (RadonMatrixWriterBase
	Options/RadonMatrixWriterOptions.cpp)

install (FILES Options/RadonMatrixWriterOptions.hpp DESTINATION include)

# add the executable
add_executable (RadonMatrixWriter
	RadonMatrixWriter.cpp)
target_link_libraries ( RadonMatrixWriter
	RadonMatrixWriterBase
	DiscretizedRadonBase
	${EXTRA_LIBS} 
	${NLOPT_LIBRARIES}
	${GSL_LIBRARY} 
	${CMAKE_THREAD_LIBS_INIT} 
	${Boost_LIBRARIES}
	${POCO_SQL_LIBRARIES}
	${LIBRT_LIBRARIES}) 
install (TARGETS RadonMatrixWriter DESTINATION bin)

add_executable (RadonMatrixWriterConfigurator 
	RadonMatrixWriterConfigurator.cpp)
target_link_libraries ( RadonMatrixWriterConfigurator 
	RadonMatrixWriterBase
	DiscretizedRadonBase
	${EXTRA_LIBS} 
	${NLOPT_LIBRARIES}
	${GSL_LIBRARY} 
	${CMAKE_THREAD_LIBS_INIT} 
	${Boost_LIBRARIES}
	${POCO_SQL_LIBRARIES}
	${LIBRT_LIBRARIES}) 
install (TARGETS RadonMatrixWriterConfigurator DESTINATION bin)
