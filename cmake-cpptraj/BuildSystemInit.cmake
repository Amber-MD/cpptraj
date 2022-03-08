# print header
# --------------------------------------------------------------------
message(STATUS "**************************************************************************")
message(STATUS "Starting configuration of ${PROJECT_NAME} version ${${PROJECT_NAME}_VERSION}...")

# configure cmake module path
# --------------------------------------------------------------------
list(
  APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}
                          "${CMAKE_CURRENT_LIST_DIR}/ThirdPartyTools"
    )

# Basic Utilities
# These files CANNOT use any sort of compile checks or system introspection
#  because no languages are enabled yet.
include(Policies NO_POLICY_SCOPE)
include(Utils)
include(Shorthand)
include(ColorMessage)
include(DebugCpptrajCmake)
