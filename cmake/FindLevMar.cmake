#
# Find the LevMar includes and library
#
# This module defines
# LEVMAR_INCLUDE_DIR, where to find tiff.h, etc.
# LEVMAR_LIBRARIES, the libraries to link against to use LevMar.
# LEVMAR_FOUND, If false, do not try to use LevMar.

# also defined, but not for general use are
# LEVMAR_LIBRARY, where to find the LevMar library.
# LEVMAR_DEBUG_LIBRARY, where to find the LevMar library in debug mode.

# taken from cmake mailing list by Anton Deguet anton.deguet at jhu.edu 

FIND_PATH(LEVMAR_INCLUDE_DIR levmar.h
  /usr/local/include
  /usr/include
)

# With Win32, important to have both
IF(WIN32)
  FIND_LIBRARY(LEVMAR_LIBRARY levmar
               ${LEVMAR_INCLUDE_DIR}/../lib
               /usr/local/lib
               /usr/lib)
  FIND_LIBRARY(LEVMAR_DEBUG_LIBRARY levmard
               ${LEVMAR_INCLUDE_DIR}/../lib
               /usr/local/lib
               /usr/lib)
ELSE(WIN32)
  # On unix system, debug and release have the same name
  FIND_LIBRARY(LEVMAR_LIBRARY levmar
               ${LEVMAR_INCLUDE_DIR}/../lib
               /usr/local/lib
               /usr/lib)
  FIND_LIBRARY(LEVMAR_DEBUG_LIBRARY levmar
               ${LEVMAR_INCLUDE_DIR}/../lib
               /usr/local/lib
               /usr/lib)
ENDIF(WIN32)

IF(LEVMAR_INCLUDE_DIR)
  IF(LEVMAR_LIBRARY)
    SET(LEVMAR_FOUND "YES")
    SET(LEVMAR_LIBRARIES ${LEVMAR_LIBRARY} ${CMAKE_DL_LIBS})
    SET(LEVMAR_DEBUG_LIBRARIES ${LEVMAR_DEBUG_LIBRARY}
${CMAKE_DL_LIBS})
  ENDIF(LEVMAR_LIBRARY)
ENDIF(LEVMAR_INCLUDE_DIR)

