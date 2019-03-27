# add the binary tree to the search path for
# include files
include_directories (${CMAKE_CURRENT_BINARY_DIR})

add_library (Math SHARED 
	Helpers.cpp)

install (TARGETS Math DESTINATION lib)

install (FILES Helpers.hpp DESTINATION include)
install (FILES MathExceptions.hpp DESTINATION include)
