# add the executable
add_executable (Basso Basso.cpp)
target_link_libraries (Basso ${EXTRA_LIBS})

# include source folder
include_directories ("${PROJECT_SOURCE_DIR}/src")

# add the install targets
install (TARGETS Basso DESTINATION bin)

