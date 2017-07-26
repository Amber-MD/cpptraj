add_subdirectory(xdrfile)

#Source files which use the USE_SANDERLIB or LIBCPPTRAJ definitions.
#These are compiled seperately for libcpptraj and cpptraj.
set(CPPTRAJ_LIBDIFFERENT_CXX_SOURCES
	Action_Esander.cpp
	Cpptraj.cpp
	Command.cpp
	Energy_Sander.cpp
	ReadLine.cpp)

#Source files which are only added to the cpptraj executable, not the library
set(CPPTRAJ_EXEONLY_CXX_SOURCES main.cpp)

#regular cpptraj C++ sources, added to cpptraj and libcpptraj 
set(CPPTRAJ_CXX_SOURCES 
	ActionFrameCounter.cpp 
	ActionList.cpp 
	Action_Align.cpp 
	Action_Angle.cpp 
	Action_AreaPerMol.cpp 
	Action_AtomMap.cpp 
	Action_AtomicCorr.cpp 
	Action_AtomicFluct.cpp 
	Action_AutoImage.cpp 
	Action_Average.cpp 
	Action_Bounds.cpp 
	Action_Box.cpp 
	Action_Center.cpp 
	Action_Channel.cpp 
	Action_CheckChirality.cpp 
	Action_CheckStructure.cpp 
	Action_Closest.cpp 
	Action_ClusterDihedral.cpp 
	Action_Contacts.cpp 
	Action_CreateCrd.cpp 
	Action_CreateReservoir.cpp 
	Action_DNAionTracker.cpp 
	Action_DSSP.cpp 
	Action_Density.cpp 
	Action_Diffusion.cpp 
	Action_Dihedral.cpp 
	Action_Dipole.cpp 
	Action_Distance.cpp 
	Action_DistRmsd.cpp 
	Action_Energy.cpp 
	Action_FilterByData.cpp 
	Action_FixAtomOrder.cpp 
	Action_FixImagedBonds.cpp 
	Action_GIST.cpp 
	Action_Grid.cpp 
	Action_GridFreeEnergy.cpp 
	Action_HydrogenBond.cpp 
	Action_Image.cpp 
	Action_Jcoupling.cpp 
	Action_LESsplit.cpp 
	Action_LIE.cpp 
	Action_MakeStructure.cpp 
	Action_Mask.cpp 
	Action_Matrix.cpp 
	Action_MinImage.cpp 
	Action_Molsurf.cpp 
	Action_MultiDihedral.cpp 
	Action_MultiVector.cpp 
	Action_NAstruct.cpp 
	Action_NativeContacts.cpp 
	Action_NMRrst.cpp 
	Action_OrderParameter.cpp 
	Action_Outtraj.cpp 
	Action_PairDist.cpp 
	Action_Pairwise.cpp 
	Action_Principal.cpp 
	Action_Projection.cpp 
	Action_Pucker.cpp 
	Action_Radgyr.cpp 
	Action_Radial.cpp 
	Action_RandomizeIons.cpp 
	Action_Remap.cpp 
	Action_ReplicateCell.cpp 
	Action_Rmsd.cpp 
	Action_Rotate.cpp 
	Action_RunningAvg.cpp 
	Action_STFC_Diffusion.cpp 
	Action_Scale.cpp 
	Action_SetVelocity.cpp 
	Action_Spam.cpp 
	Action_Strip.cpp 
	Action_Surf.cpp 
	Action_SymmetricRmsd.cpp 
	Action_Temperature.cpp 
	Action_Translate.cpp 
	Action_Unstrip.cpp 
	Action_Unwrap.cpp 
	Action_Vector.cpp 
	Action_VelocityAutoCorr.cpp 
	Action_Volmap.cpp 
	Action_Volume.cpp 
	Action_Watershell.cpp 
	AnalysisList.cpp 
	Analysis_AmdBias.cpp 
	Analysis_AutoCorr.cpp 
	Analysis_Average.cpp 
	Analysis_Clustering.cpp 
	Analysis_Corr.cpp 
	Analysis_CrankShaft.cpp 
	Analysis_CrdFluct.cpp 
	Analysis_CrossCorr.cpp 
	Analysis_CurveFit.cpp 
	Analysis_Divergence.cpp 
	Analysis_FFT.cpp 
	Analysis_Hist.cpp 
	Analysis_Integrate.cpp 
	Analysis_IRED.cpp 
	Analysis_KDE.cpp 
	Analysis_Lifetime.cpp 
	Analysis_LowestCurve.cpp 
	Analysis_Matrix.cpp 
	Analysis_MeltCurve.cpp 
	Analysis_Modes.cpp 
	Analysis_Multicurve.cpp 
	Analysis_MultiHist.cpp 
	Analysis_Overlap.cpp 
	Analysis_PhiPsi.cpp 
	Analysis_Regression.cpp 
	Analysis_RemLog.cpp 
	Analysis_Rms2d.cpp 
	Analysis_RmsAvgCorr.cpp 
	Analysis_Rotdif.cpp 
	Analysis_RunningAvg.cpp 
	Analysis_Spline.cpp 
	Analysis_State.cpp 
	Analysis_Statistics.cpp 
	Analysis_TI.cpp 
	Analysis_Timecorr.cpp 
	Analysis_VectorMath.cpp 
	Analysis_Wavelet.cpp 
	ArgList.cpp 
	Array1D.cpp 
	AssociatedData.cpp 
	Atom.cpp 
	AtomMap.cpp 
	AtomMask.cpp 
	AxisType.cpp 
	BondSearch.cpp 
	Box.cpp 
	BufferedFrame.cpp 
	BufferedLine.cpp 
	ByteRoutines.cpp 
	CharMask.cpp 
	CIFfile.cpp 
	ClusterDist.cpp 
	ClusterList.cpp 
	ClusterMap.cpp 
	ClusterMatrix.cpp 
	ClusterNode.cpp 
	ClusterSieve.cpp 
	Cluster_DBSCAN.cpp 
	Cluster_DPeaks.cpp 
	Cluster_HierAgglo.cpp 
	Cluster_Kmeans.cpp 
	Cluster_ReadInfo.cpp 
	Cmd.cpp 
	CmdInput.cpp 
	CmdList.cpp 
	ComplexArray.cpp 
	Constraints.cpp 
	CoordinateInfo.cpp 
	Corr.cpp 
	CpptrajFile.cpp 
	CpptrajState.cpp 
	CpptrajStdio.cpp 
	CurveFit.cpp 
	DataFile.cpp 
	DataFileList.cpp 
	DataIO.cpp 
	DataIO_CCP4.cpp 
	DataIO_CharmmFastRep.cpp 
	DataIO_CharmmRepLog.cpp 
	DataIO_Cmatrix.cpp 
	DataIO_Evecs.cpp 
	DataIO_Gnuplot.cpp 
	DataIO_Grace.cpp 
	DataIO_OpenDx.cpp 
	DataIO_Mdout.cpp 
	DataIO_NC_Cmatrix.cpp 
	DataIO_RemLog.cpp 
	DataIO_Std.cpp 
	DataIO_VecTraj.cpp 
	DataIO_Xplor.cpp 
	DataIO_XVG.cpp 
	DataSet.cpp 
	DataSetList.cpp 
	DataSet_1D.cpp 
	DataSet_3D.cpp 
	DataSet_Cmatrix.cpp 
	DataSet_Cmatrix_DISK.cpp 
	DataSet_Cmatrix_MEM.cpp 
	DataSet_Cmatrix_NOMEM.cpp 
	DataSet_Coords.cpp 
	DataSet_Coords_CRD.cpp 
	DataSet_Coords_REF.cpp 
	DataSet_Coords_TRJ.cpp 
	DataSet_GridDbl.cpp 
	DataSet_GridFlt.cpp 
	DataSet_Mat3x3.cpp 
	DataSet_MatrixDbl.cpp 
	DataSet_MatrixFlt.cpp 
	DataSet_Mesh.cpp 
	DataSet_Modes.cpp 
	DataSet_RemLog.cpp 
	DataSet_Topology.cpp 
	DataSet_Vector.cpp 
	DataSet_double.cpp 
	DataSet_float.cpp 
	DataSet_integer.cpp 
	DataSet_string.cpp 
	Deprecated.cpp 
	DihedralSearch.cpp 
	DispatchObject.cpp 
	DistRoutines.cpp 
	Energy.cpp 
	EnsembleIn.cpp 
	EnsembleIn_Multi.cpp 
	EnsembleIn_Single.cpp 
	EnsembleNavigator.cpp 
	EnsembleOut.cpp 
	EnsembleOut_Multi.cpp 
	EnsembleOut_Single.cpp 
	EnsembleOutList.cpp 
	Ewald.cpp 
	Exec_Analyze.cpp 
	Exec_Calc.cpp 
	Exec_Change.cpp 
	Exec_ClusterMap.cpp 
	Exec_CombineCoords.cpp 
	Exec_Commands.cpp 
	Exec_CompareTop.cpp 
	Exec_CrdAction.cpp 
	Exec_CrdOut.cpp 
	Exec_DataFile.cpp 
	Exec_DataFilter.cpp 
	Exec_DataSetCmd.cpp 
	Exec_GenerateAmberRst.cpp 
	Exec_Help.cpp 
	Exec_LoadCrd.cpp 
	Exec_LoadTraj.cpp 
	Exec_ParmBox.cpp 
	Exec_ParmSolvent.cpp 
	Exec_ParmStrip.cpp 
	Exec_ParmWrite.cpp 
	Exec_PermuteDihedrals.cpp 
	Exec_Precision.cpp 
	Exec_PrintData.cpp 
	Exec_ReadData.cpp 
	Exec_ReadInput.cpp 
	Exec_RotateDihedral.cpp 
	Exec_RunAnalysis.cpp 
	Exec_ScaleDihedralK.cpp 
	Exec_SequenceAlign.cpp 
	Exec_System.cpp 
	Exec_Top.cpp 
	Exec_Traj.cpp 
	Exec_ViewRst.cpp 
	FileIO_Bzip2.cpp 
	FileIO_Gzip.cpp 
	FileIO_Mpi.cpp 
	FileIO_Std.cpp 
	FileName.cpp 
	FileTypes.cpp 
	Frame.cpp 
	GridAction.cpp 
	HistBin.cpp 
	Hungarian.cpp 
	ImageRoutines.cpp 
	InputTrajCommon.cpp 
	KDE.cpp 
	MapAtom.cpp 
	MaskToken.cpp 
	Matrix_3x3.cpp 
	MetaData.cpp 
	Mol2File.cpp 
	NameType.cpp 
	NC_Cmatrix.cpp 
	NC_Routines.cpp 
	NetcdfFile.cpp 
	OutputTrajCommon.cpp 
	PDBfile.cpp 
	PairList.cpp 
	Parallel.cpp 
	ParallelNetcdf.cpp 
	ParmFile.cpp 
	Parm_Amber.cpp 
	Parm_CharmmPsf.cpp 
	Parm_CIF.cpp 
	Parm_Gromacs.cpp 
	Parm_Mol2.cpp 
	Parm_PDB.cpp 
	Parm_SDF.cpp 
	Parm_Tinker.cpp 
	ProgressBar.cpp 
	ProgressTimer.cpp 
	PubFFT.cpp 
	Random.cpp 
	Range.cpp 
	RPNcalc.cpp 
	ReferenceAction.cpp 
	Residue.cpp 
	SDFfile.cpp 
	SimplexMin.cpp 
	Spline.cpp 
	StringRoutines.cpp 
	StructureMapper.cpp 
	SymmetricRmsdCalc.cpp 
	TextFormat.cpp 
	Timer.cpp 
	TinkerFile.cpp 
	TopInfo.cpp 
	Topology.cpp 
	TorsionRoutines.cpp 
	Traj_AmberCoord.cpp 
	Traj_AmberNetcdf.cpp 
	Traj_AmberRestart.cpp 
	Traj_AmberRestartNC.cpp 
	Traj_Binpos.cpp 
	Traj_CharmmCor.cpp 
	Traj_CharmmDcd.cpp 
	Traj_CharmmRestart.cpp 
	Traj_CIF.cpp 
	Traj_Conflib.cpp 
	Traj_GmxTrX.cpp 
	Traj_GmxXtc.cpp 
	Traj_Gro.cpp 
	Traj_Mol2File.cpp 
	Traj_NcEnsemble.cpp 
	Traj_PDBfile.cpp 
	Traj_SDF.cpp 
	Traj_SQM.cpp 
	Traj_Tinker.cpp 
	TrajectoryFile.cpp 
	TrajectoryIO.cpp 
	TrajFrameCounter.cpp 
	TrajinList.cpp 
	Trajin_Multi.cpp 
	Trajin_Single.cpp 
	TrajIOarray.cpp 
	Trajout_Single.cpp 
	TrajoutList.cpp 
	Vec3.cpp 
	ViewRst.cpp)
        
set(CPPTRAJ_C_SOURCES molsurf.c)

set(CPPTRAJ_Fortran_SOURCES pub_fft.F90)

set(AMBPDB_CXX_SOURCES
	AmbPDB.cpp 
	ActionFrameCounter.cpp 
	ArgList.cpp 
	Atom.cpp 
	AtomMask.cpp 
	BondSearch.cpp 
	Box.cpp 
	BufferedFrame.cpp 
	BufferedLine.cpp 
	ByteRoutines.cpp 
	CharMask.cpp 
	CIFfile.cpp 
	CoordinateInfo.cpp 
	CpptrajFile.cpp 
	CpptrajStdio.cpp 
	DistRoutines.cpp 
	FileIO_Bzip2.cpp 
	FileIO_Gzip.cpp 
	FileIO_Std.cpp 
	FileName.cpp 
	FileTypes.cpp 
	Frame.cpp 
	InputTrajCommon.cpp 
	MaskToken.cpp 
	Matrix_3x3.cpp 
	Mol2File.cpp 
	NameType.cpp 
	NC_Routines.cpp 
	NetcdfFile.cpp 
	OutputTrajCommon.cpp 
	ParmFile.cpp 
	Parm_Amber.cpp 
	Parm_CIF.cpp 
	Parm_CharmmPsf.cpp 
	Parm_Gromacs.cpp 
	Parm_Mol2.cpp 
	Parm_PDB.cpp 
	Parm_SDF.cpp 
	Parm_Tinker.cpp 
	PDBfile.cpp 
	ProgressBar.cpp 
	Range.cpp 
	SDFfile.cpp 
	StringRoutines.cpp 
	TextFormat.cpp 
	TinkerFile.cpp 
	Timer.cpp 
	Topology.cpp 
	TrajectoryFile.cpp 
	TrajFrameCounter.cpp 
	Trajin_Single.cpp 
	Trajout_Single.cpp 
	Residue.cpp 
	Traj_AmberCoord.cpp 
	Traj_AmberRestart.cpp 
	Traj_CIF.cpp 
	Traj_Mol2File.cpp 
	Traj_SDF.cpp 
	Traj_AmberNetcdf.cpp 
	Traj_Binpos.cpp 
	Traj_Conflib.cpp 
	Traj_SQM.cpp 
	Traj_AmberRestartNC.cpp 
	Traj_CharmmCor.cpp 
	Traj_CharmmDcd.cpp 
	Traj_CharmmRestart.cpp 
	Traj_GmxTrX.cpp 
	Traj_GmxXtc.cpp 
	Traj_Gro.cpp 
	Traj_PDBfile.cpp 
	Traj_Tinker.cpp 
	Traj_NcEnsemble.cpp 
	Vec3.cpp)
	
set_property(SOURCE ${CPPTRAJ_Fortran_SOURCES} PROPERTY COMPILE_FLAGS "${OPT_FFLAGS_SPC} ${OpenMP_Fortran_FLAGS} ${MPI_Fortran_COMPILE_FLAGS}")
set_property(SOURCE ${CPPTRAJ_CXX_SOURCES} ${AMBPDB_CXX_SOURCES} ${CPPTRAJ_LIBDIFFERENT_CXX_SOURCES} ${CPPTRAJ_EXEONLY_CXX_SOURCES}
	PROPERTY COMPILE_FLAGS "${OPT_CXXFLAGS_SPC} ${OpenMP_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
set_property(SOURCE ${CPPTRAJ_C_SOURCES} PROPERTY COMPILE_FLAGS "${OPT_CFLAGS_SPC} ${OpenMP_C_FLAGS} ${MPI_C_COMPILE_FLAGS}")

include_directories(${AMBERTOOLS_INC_DIR})

if(MPI)
	include_directories(${MPI_CXX_INCLUDE_PATH})
endif()
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
#add the common object library

#concatenate all the source files
set(CPPTRAJ_COMMON_SOURCES ${CPPTRAJ_CXX_SOURCES} ${CPPTRAJ_C_SOURCES} ${CPPTRAJ_Fortran_SOURCES})

add_library(cpptraj_common_obj OBJECT ${CPPTRAJ_COMMON_SOURCES})
set_property(TARGET cpptraj_common_obj PROPERTY POSITION_INDEPENDENT_CODE ${SHARED})

#normally this would be applied by target_link_libraries, but since we use that intermediary object library, we have to apply it manually

# NOTE: there is a CMake bug where if we were to set these as a directory-scope includes, the CUDA build would fail (on some platforms only???)
# it turns out that CMake's cuda library passes the include paths after the first one from each of these generator expressions to nvcc without the -I flag
# This causes the error "A single input file is required for a non-link phase when an outputfile is specified"
target_include_directories(cpptraj_common_obj PRIVATE $<TARGET_PROPERTY:xdrfile,INTERFACE_INCLUDE_DIRECTORIES> $<TARGET_PROPERTY:netcdf,INTERFACE_INCLUDE_DIRECTORIES>)
if(fftw_ENABLED)
	target_include_directories(cpptraj_common_obj PRIVATE $<TARGET_PROPERTY:fftw,INTERFACE_INCLUDE_DIRECTORIES>)
endif()
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------

add_executable(cpptraj $<TARGET_OBJECTS:cpptraj_common_obj> ${CPPTRAJ_LIBDIFFERENT_CXX_SOURCES} ${CPPTRAJ_EXEONLY_CXX_SOURCES})
if(BUILD_SANDER_API)
	#add the sander-specific definitions and libraries
	set_property(TARGET cpptraj PROPERTY COMPILE_DEFINITIONS USE_SANDERLIB)
	target_link_libraries(cpptraj libsander)
endif()

set_property(TARGET cpptraj PROPERTY LINK_FLAGS ${OpenMP_CXX_FLAGS})
target_link_libraries(cpptraj netcdf netlib)

if(MPI)
	set_property(TARGET cpptraj PROPERTY OUTPUT_NAME cpptraj.MPI)
elseif(OPENMP)
	set_property(TARGET cpptraj PROPERTY OUTPUT_NAME cpptraj.OMP)
endif()

#------------------------------------------------------------------------------------------

add_executable(ambpdb ${AMBPDB_CXX_SOURCES})
target_link_libraries(ambpdb netcdf)
install(TARGETS cpptraj ambpdb DESTINATION ${BINDIR})
set_property(TARGET ambpdb PROPERTY LINK_FLAGS ${OpenMP_CXX_FLAGS})

#libcpptraj
#------------------------------------------------------------------------------------------

add_library(libcpptraj $<TARGET_OBJECTS:cpptraj_common_obj> ${CPPTRAJ_LIBDIFFERENT_CXX_SOURCES})
set_property(TARGET libcpptraj PROPERTY COMPILE_DEFINITIONS LIBCPPTRAJ)
target_link_libraries(libcpptraj netcdf netlib)
remove_prefix(libcpptraj)
install_libraries(libcpptraj)
set_property(TARGET libcpptraj PROPERTY LINK_FLAGS ${OpenMP_CXX_FLAGS})

#tell others where to find the cpptraj includes
target_include_directories(libcpptraj INTERFACE .)

#Deal with external libraries
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------

# NOTE: you CANNOT set any directory-scope include directories that use generator expressions here.
# These expressions get propagated down into the cuda_kernels subdir, and trigger an undocumented bug in old versions of CMake
# where they get passed to nvcc without a -I prefix, breaking everything.  So please don't do it.

if(libbz2_ENABLED)
	add_definitions(-DHASBZ2)
	include_directories(${BZIP2_INCLUDE_DIR})
	targets_link_libraries(cpptraj libcpptraj ambpdb LIBRARIES bzip2)
endif()

if(zlib_ENABLED)
	add_definitions(-DHASGZ)
	include_directories(${ZLIB_INCLUDE_DIRS})
	targets_link_libraries(cpptraj libcpptraj ambpdb LIBRARIES ${ZLIB_LIBRARIES})
endif()

#readline
if(readline_ENABLED)
	targets_link_libraries(cpptraj libcpptraj LIBRARIES readline)
else()
	target_compile_definitions(cpptraj PRIVATE NO_READLINE)
	target_compile_definitions(libcpptraj PRIVATE NO_READLINE)
endif()

if(MPI)
	#don't give the MPI definition to ambpdb
	target_compile_definitions(cpptraj PRIVATE MPI)
	target_compile_definitions(libcpptraj PRIVATE MPI)
	target_compile_definitions(cpptraj_common_obj PRIVATE MPI)
		
	if(pnetcdf_ENABLED)
		targets_link_libraries(cpptraj libcpptraj LIBRARIES pnetcdf)
		target_include_directories(cpptraj_common_obj PUBLIC $<TARGET_PROPERTY:pnetcdf,INTERFACE_INCLUDE_DIRECTORIES>)
		add_definitions(-DHAS_PNETCDF)
	endif()
		
	set_property(TARGET cpptraj PROPERTY OUTPUT_FILENAME cpptraj.MPI)	
endif()
link_mpi(cpptraj CXX)
link_mpi(libcpptraj CXX)

enable_openmp(cpptraj CXX)
enable_openmp(libcpptraj CXX)

#xdrfile
targets_link_libraries(cpptraj libcpptraj ambpdb LIBRARIES xdrfile)


# CUDA
if(CUDA)
	add_subdirectory(cuda_kernels)
	
	add_definitions(-DCUDA)
	include_directories(${CUDA_INCLUDE_DIRS})
	targets_link_libraries(cpptraj libcpptraj LIBRARIES cpptraj_cuda)
endif()

# FFTW
if(fftw_ENABLED)	
	set_property(SOURCE PubFFT.cpp PROPERTY COMPILE_DEFINITIONS FFTW_FFT)
	targets_link_libraries(cpptraj libcpptraj LIBRARIES fftw)
endif()

# arpack
if(arpack_DISABLED)
	add_definitions(-DNO_ARPACK)
endif()