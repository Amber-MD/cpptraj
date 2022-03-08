# print header
# --------------------------------------------------------------------
message(STATUS "**************************************************************************")
message(STATUS "Starting configuration of ${PROJECT_NAME} version ${${PROJECT_NAME}_VERSION}...")

# configure cmake module path
# --------------------------------------------------------------------
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})

# Basic Utilities
include(Utils)
include(Shorthand)
