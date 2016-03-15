# add the binary tree to the search path for
# include files
include_directories (${CMAKE_CURRENT_BINARY_DIR})

add_library (DiscretizedRadonBase
	DiscretizedRadon/Backprojection.cpp
	DiscretizedRadon/BackprojectionMatrix.cpp
	DiscretizedRadon/DiscretizedRadon.cpp
	DiscretizedRadon/DiscretizedRadonMatrix.cpp
	DiscretizedRadon/Helpers.cpp
	DiscretizedRadon/point_t.cpp)

install (FILES DiscretizedRadon/Backprojection.hpp DESTINATION include)
install (FILES DiscretizedRadon/BackprojectionMatrix.hpp DESTINATION include)
install (FILES DiscretizedRadon/DiscretizedRadon.hpp DESTINATION include)
install (FILES DiscretizedRadon/DiscretizedRadonMatrix.hpp DESTINATION include)
install (FILES DiscretizedRadon/Helpers.hpp DESTINATION include)
install (FILES DiscretizedRadon/point_t.hpp DESTINATION include)

add_library (ComputerTomographyBase
	Options/ComputerTomographyOptions.cpp)

install (FILES Options/ComputerTomographyOptions.hpp DESTINATION include)

# add the executable
add_executable (ComputerTomography
	ComputerTomography.cpp)
target_link_libraries ( ComputerTomography
	ComputerTomographyBase
	DiscretizedRadonBase
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
	DiscretizedRadonBase
	${EXTRA_LIBS} 
	${NLOPT_LIBRARIES}
	${GSL_LIBRARY} 
	${CMAKE_THREAD_LIBS_INIT} 
	${Boost_LIBRARIES}
	${POCO_SQL_LIBRARIES}
	${LIBRT_LIBRARIES}) 
install (TARGETS ComputerTomographyConfigurator DESTINATION bin)
