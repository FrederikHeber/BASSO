# add the binary tree to the search path for
# include files
include_directories (${CMAKE_CURRENT_BINARY_DIR})

add_library (Database Database.cpp Database.hpp)

install (FILES Database.hpp DESTINATION include)
