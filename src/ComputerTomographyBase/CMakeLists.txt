# add the binary tree to the search path for
# include files
include_directories (${CMAKE_CURRENT_BINARY_DIR})

add_library (ComputerTomographyBase
	DiscretizedRadon/BackprojectionMatrix.cpp
	DiscretizedRadon/DiscretizedRadonMatrix.cpp
	DiscretizedRadon/point_t.cpp
	Options/ComputerTomographyOptions.cpp)

install (FILES DiscretizedRadon/BackprojectionMatrix.hpp DESTINATION include)
install (FILES DiscretizedRadon/DiscretizedRadonMatrix.hpp DESTINATION include)
install (FILES DiscretizedRadon/point_t.hpp DESTINATION include)
install (FILES Options/ComputerTomographyOptions.hpp DESTINATION include)

# add the executable
add_executable (ComputerTomography
	ComputerTomography.cpp)
target_link_libraries ( ComputerTomography
	ComputerTomographyBase
	${EXTRA_LIBS} 
	${NLOPT_LIBRARIES}
	${GSL_LIBRARY} 
	${CMAKE_THREAD_LIBS_INIT} 
	${Boost_LIBRARIES}
	${POCO_SQL_LIBRARIES}
	${LIBRT_LIBRARIES}) 
install (TARGETS ComputerTomography DESTINATION bin)

add_executable (ComputerTomographyConfigurator 
	ComputerTomographyConfigurator.cpp)
target_link_libraries ( ComputerTomographyConfigurator 
	ComputerTomographyBase
	${EXTRA_LIBS} 
	${NLOPT_LIBRARIES}
	${GSL_LIBRARY} 
	${CMAKE_THREAD_LIBS_INIT} 
	${Boost_LIBRARIES}
	${POCO_SQL_LIBRARIES}
	${LIBRT_LIBRARIES}) 
install (TARGETS ComputerTomographyConfigurator DESTINATION bin)
