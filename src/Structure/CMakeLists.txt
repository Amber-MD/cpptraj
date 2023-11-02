#CMake buildfile for CPPTRAJ Structure subdirectory.
target_sources(cpptraj_common_obj PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/Builder.cpp
  ${CMAKE_CURRENT_LIST_DIR}/Chirality.cpp
  ${CMAKE_CURRENT_LIST_DIR}/Disulfide.cpp
  ${CMAKE_CURRENT_LIST_DIR}/FunctionalGroup.cpp
  ${CMAKE_CURRENT_LIST_DIR}/FxnGroupBuilder.cpp
  ${CMAKE_CURRENT_LIST_DIR}/HisProt.cpp
  ${CMAKE_CURRENT_LIST_DIR}/InternalCoords.cpp
  ${CMAKE_CURRENT_LIST_DIR}/Model.cpp
  ${CMAKE_CURRENT_LIST_DIR}/StructureEnum.cpp
  ${CMAKE_CURRENT_LIST_DIR}/StructureRoutines.cpp
  ${CMAKE_CURRENT_LIST_DIR}/SugarBuilder.cpp
  ${CMAKE_CURRENT_LIST_DIR}/Sugar.cpp
  ${CMAKE_CURRENT_LIST_DIR}/SugarLinkAtoms.cpp
  ${CMAKE_CURRENT_LIST_DIR}/SugarToken.cpp
  ${CMAKE_CURRENT_LIST_DIR}/Zmatrix.cpp
)
