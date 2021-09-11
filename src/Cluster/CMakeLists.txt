#CMake buildfile for CPPTRAJ Cluster subdirectory.
set(CLUSTER_SOURCES
  Algorithm_DBscan.cpp
  Algorithm_DPeaks.cpp
  Algorithm_HierAgglo.cpp
  Algorithm_Kmeans.cpp
  BestReps.cpp
  Centroid_Coord.cpp
  Cframes.cpp
  Cmatrix_Binary.cpp
  Cmatrix_NC.cpp
  Control.cpp
  DrawGraph.cpp
  DynamicMatrix.cpp
  List.cpp
  Metric_Data.cpp
  Metric_Data_Euclid.cpp
  Metric_Data_Manhattan.cpp
  Metric_DME.cpp
  Metric_RMS.cpp
  Metric_SRMSD.cpp
  Node.cpp
  Output.cpp
  PairwiseMatrix.cpp
  Results_Coords.cpp
  Sieve.cpp
)

#------------------------------------------------------------------------------------------
#I assume that this library should be optimized
add_compile_options(${OPT_CFLAGS})

#if("${CMAKE_C_COMPILER_ID}" STREQUAL "MSVC")
#	# xdrfile is not set up to build as a Windows shared library.  Someone would need to add dll exports to all of the symbols (bluh)
#	add_library(xdrfile STATIC ${CLUSTER_SOURCES})
#else()
	add_library(Cluster ${CLUSTER_SOURCES})
#endif()

target_include_directories(Cluster PUBLIC .)
install_libraries(Cluster)
