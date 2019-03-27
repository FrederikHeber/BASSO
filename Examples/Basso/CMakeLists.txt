# add the binary tree to the search path for
# include files
include_directories (${CMAKE_CURRENT_BINARY_DIR})

add_library (BassoBase
	Options/BassoOptions.cpp
	printCounts.cpp)

install (FILES Options/BassoOptions.hpp DESTINATION include)
install (FILES printCounts.hpp DESTINATION include)

# add the executable
add_executable (Basso
	Basso.cpp)
target_link_libraries ( Basso
	BassoBase
	${EXTRA_LIBS} 
	${NLOPT_LIBRARIES}
	${GSL_LIBRARY} 
	${CMAKE_THREAD_LIBS_INIT} 
	${BOOST_CORE_LIBRARIES}
	${POCO_SQL_LIBRARIES}
	${LIBRT_LIBRARIES}) 
install (TARGETS Basso DESTINATION bin)

add_executable (BassoConfigurator 
	BassoConfigurator.cpp)
target_link_libraries ( BassoConfigurator 
	BassoBase
	${EXTRA_LIBS} 
	${NLOPT_LIBRARIES}
	${GSL_LIBRARY} 
	${CMAKE_THREAD_LIBS_INIT} 
	${BOOST_CORE_LIBRARIES}
	${POCO_SQL_LIBRARIES}
	${LIBRT_LIBRARIES}) 
install (TARGETS BassoConfigurator DESTINATION bin)
