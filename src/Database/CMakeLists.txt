# add the binary tree to the search path for
# include files
include_directories (${CMAKE_CURRENT_BINARY_DIR})

add_library (Database Database.cpp Table.cpp Table_Tuple.cpp)

install (FILES Database.hpp DESTINATION include)
install (FILES DefaultValue.hpp DESTINATION include)
install (FILES Table.hpp DESTINATION include)
