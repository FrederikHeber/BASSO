# add the binary tree to the search path for
# include files
include_directories (${CMAKE_CURRENT_BINARY_DIR})

add_library (Log SHARED
	Logging.cpp 
	Verbose.cpp)

target_link_libraries(Log
	${BOOST_CORE_LIBRARIES}
)

ADD_DEPENDENCIES(Log version)

install (TARGETS Log DESTINATION lib)

install (FILES Logging.hpp DESTINATION include)
install (FILES Verbose.hpp DESTINATION include)

