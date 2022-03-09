# print header
# --------------------------------------------------------------------
message(STATUS "**************************************************************************")
message(STATUS "Starting configuration of ${PROJECT_NAME} version ${${PROJECT_NAME}_VERSION}...")

# print CMake version at the start so we can use it to diagnose issues even if the configure fails
message(STATUS "CMake Version: ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}")

# fix search path so that libraries from the install tree are not used
# --------------------------------------------------------------------
list(REMOVE_ITEM CMAKE_SYSTEM_PREFIX_PATH "${CMAKE_INSTALL_PREFIX}")

# eliminate extraneous install messages
# --------------------------------------------------------------------
set(CMAKE_INSTALL_MESSAGE LAZY)

# configure module path
# --------------------------------------------------------------------
list(
  APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}
                          "${CMAKE_CURRENT_LIST_DIR}/ThirdPartyTools"
    )

# prevent obliteration of the old build system's makefiles
# --------------------------------------------------------------------

if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
	message(FATAL_ERROR "You are building in the source directory.  ${PROJECT_NAME} does not support this, since it would obliterate the Makefile build system.")
endif()
# Basic Utilities
# These files CANNOT use any sort of compile checks or system introspection
#  because no languages are enabled yet.
include(Policies NO_POLICY_SCOPE)
include(Utils)
include(Shorthand)
include(ColorMessage)
include(DebugCpptrajCmake)

# get install directories
include(InstallDirs)
