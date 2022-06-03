#CMake buildfile for CPPTRAJ Structure subdirectory.
target_sources(cpptraj_common_obj PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/Sugar.cpp
  ${CMAKE_CURRENT_LIST_DIR}/SugarLinkAtoms.cpp
  ${CMAKE_CURRENT_LIST_DIR}/SugarToken.cpp
)
