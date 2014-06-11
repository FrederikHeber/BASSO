# add the binary tree to the search path for
# include files
include_directories (${CMAKE_CURRENT_BINARY_DIR})

add_library (Log Logging.cpp Logging.hpp)

install (FILES Logging.hpp DESTINATION include)
