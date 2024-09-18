#CMake buildfile for CPPTRAJ Energy subdirectory.
target_sources(cpptraj_common_obj PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/EnergyDecomp_Ewald.cpp
  ${CMAKE_CURRENT_LIST_DIR}/EnergyDecomposer.cpp
  ${CMAKE_CURRENT_LIST_DIR}/ErfcFxn.cpp
  ${CMAKE_CURRENT_LIST_DIR}/EwaldParams.cpp
)
