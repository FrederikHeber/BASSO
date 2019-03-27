# add the binary tree to the search path for
# include files
include_directories (${CMAKE_CURRENT_BINARY_DIR})

add_library (NoiseAdderBase
	Options/NoiseAdderOptions.cpp)

install (FILES Options/NoiseAdderOptions.hpp DESTINATION include)

# add the executable
add_executable (NoiseAdder
	NoiseAdder.cpp)
target_link_libraries ( NoiseAdder
	NoiseAdderBase
	${EXTRA_LIBS} 
	${NLOPT_LIBRARIES}
	${GSL_LIBRARY} 
	${CMAKE_THREAD_LIBS_INIT} 
	${BOOST_CORE_LIBRARIES}
	${POCO_SQL_LIBRARIES}
	${LIBRT_LIBRARIES}) 
install (TARGETS NoiseAdder DESTINATION bin)

add_executable (NoiseAdderConfigurator 
	NoiseAdderConfigurator.cpp)
target_link_libraries ( NoiseAdderConfigurator 
	NoiseAdderBase
	${EXTRA_LIBS} 
	${NLOPT_LIBRARIES}
	${GSL_LIBRARY} 
	${CMAKE_THREAD_LIBS_INIT} 
	${BOOST_CORE_LIBRARIES}
	${POCO_SQL_LIBRARIES}
	${LIBRT_LIBRARIES}) 
install (TARGETS NoiseAdderConfigurator DESTINATION bin)
