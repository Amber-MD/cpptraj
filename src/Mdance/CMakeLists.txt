#CMake buildfile for CPPTRAJ Mdance subdirectory.
target_sources(cpptraj_common_obj PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/bts.cpp
  ${CMAKE_CURRENT_LIST_DIR}/cluster.cpp
  ${CMAKE_CURRENT_LIST_DIR}/esim.cpp
  ${CMAKE_CURRENT_LIST_DIR}/hc_utils.cpp
  ${CMAKE_CURRENT_LIST_DIR}/helm.cpp
  ${CMAKE_CURRENT_LIST_DIR}/KMeans.cpp
  ${CMAKE_CURRENT_LIST_DIR}/result_utils.cpp
  ${CMAKE_CURRENT_LIST_DIR}/scores.cpp
)
