#CMake buildfile for CPPTRAJ Energy subdirectory.
target_sources(cpptraj_common_obj PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/EnergyDecomposer.cpp
  ${CMAKE_CURRENT_LIST_DIR}/ErfcFxn.cpp
  ${CMAKE_CURRENT_LIST_DIR}/EwaldCalc_LJPME.cpp
  ${CMAKE_CURRENT_LIST_DIR}/EwaldCalc_PME.cpp
  ${CMAKE_CURRENT_LIST_DIR}/EwaldParams.cpp
  ${CMAKE_CURRENT_LIST_DIR}/EwaldParams_LJPME.cpp
  ${CMAKE_CURRENT_LIST_DIR}/EwaldParams_PME.cpp
  ${CMAKE_CURRENT_LIST_DIR}/PME_Recip.cpp
  ${CMAKE_CURRENT_LIST_DIR}/VDW_LongRange_Correction.cpp
)
