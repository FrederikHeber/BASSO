# add the binary tree to the search path for
# include files
include_directories (${CMAKE_CURRENT_BINARY_DIR})

add_library (MatrixToPNGBase
	Colortables/ColorTable.cpp
	Options/MatrixToPNGOptions.cpp)

install (FILES Colortables/ColorTable.hpp DESTINATION include)
install (FILES Options/MatrixToPNGOptions.hpp DESTINATION include)

# add the executable
add_executable (MatrixToPNG
	MatrixToPNG.cpp)
target_link_libraries ( MatrixToPNG
	MatrixToPNGBase
	Options
	Log
	${PNG_LIBRARIES}
	${BOOST_CORE_LIBRARIES}
	${BOOST_THREAD_LIBRARY}
	${POCO_SQL_LIBRARIES}
	${LIBRT_LIBRARIES}) 
install (TARGETS MatrixToPNG DESTINATION bin)
