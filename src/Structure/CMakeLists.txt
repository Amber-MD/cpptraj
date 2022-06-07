#CMake buildfile for CPPTRAJ Structure subdirectory.
target_sources(cpptraj_common_obj PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/Disulfide.cpp
  ${CMAKE_CURRENT_LIST_DIR}/FunctionalGroup.cpp
  ${CMAKE_CURRENT_LIST_DIR}/FxnGroupBuilder.cpp
  ${CMAKE_CURRENT_LIST_DIR}/StructureRoutines.cpp
  ${CMAKE_CURRENT_LIST_DIR}/SugarBuilder.cpp
  ${CMAKE_CURRENT_LIST_DIR}/Sugar.cpp
  ${CMAKE_CURRENT_LIST_DIR}/SugarLinkAtoms.cpp
  ${CMAKE_CURRENT_LIST_DIR}/SugarToken.cpp
)
