#CMake buildfile for CPPTRAJ Parm subdirectory.
target_sources(cpptraj_common_obj PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/AssignParams.cpp
  ${CMAKE_CURRENT_LIST_DIR}/Bond_Params.cpp
  ${CMAKE_CURRENT_LIST_DIR}/GB_Params.cpp
  ${CMAKE_CURRENT_LIST_DIR}/GetParams.cpp
  ${CMAKE_CURRENT_LIST_DIR}/LJ1264_Params.cpp
  ${CMAKE_CURRENT_LIST_DIR}/Merge.cpp
  ${CMAKE_CURRENT_LIST_DIR}/ParameterSet.cpp
  ${CMAKE_CURRENT_LIST_DIR}/ParmEnum.cpp
)
