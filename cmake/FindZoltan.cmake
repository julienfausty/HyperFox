# - Try to find Zoltan
# Once done this will define
#
#  Zoltan_FOUND        - system has Zoltan
#  Zoltan_INCLUDES     - the Zoltan include directories
#  Zoltan_LIBRARIES    - Link these to use Zoltan
#  Zoltan_DEFINITIONS  - Compiler switches for using Zoltan
#
#  Usage:
#  find_package(Zoltan)
#
# Setting these changes the behavior of the search
#  Zoltan_PREFIX - directory in which Zoltan resides (mandatory)
#
# This find is largely inspired from the SCOREC/core version of FindZoltan.cmake
#

set(ZOLTAN_PREFIX "$ENV{ZOLTAN_PREFIX}" CACHE STRING "Zoltan install directory")
if(ZOLTAN_PREFIX)
  message(STATUS "ZOLTAN_PREFIX ${ZOLTAN_PREFIX}")
else()
  message(FATAL_ERROR "ZOLTAN_PREFIX is empty, it must be defined as an environment variable (vanilla Zoltan installation is /usr/local)")
endif()

find_path(ZOLTAN_INCLUDE_DIR zoltan.h PATHS "${ZOLTAN_PREFIX}/include")

find_library(ZOLTAN_LIBRARY zoltan PATHS "${ZOLTAN_PREFIX}/lib")

set(ZOLTAN_LIBRARIES ${ZOLTAN_LIBRARY} )
set(ZOLTAN_INCLUDES ${ZOLTAN_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set ZOLTAN_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(
    Zoltan
    DEFAULT_MSG
    ZOLTAN_LIBRARY ZOLTAN_INCLUDE_DIR
)

mark_as_advanced(ZOLTAN_INCLUDES ZOLTAN_LIBRARIES)
