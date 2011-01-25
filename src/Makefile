include ../config.h



SOURCE=main.cpp AmberNetcdf.cpp AmberParm.cpp AmberTraj.cpp TrajFile.cpp \
       Frame.cpp PtrajState.cpp ArgList.cpp ptrajmask.c Action_Distance.cpp \
       ParmFileList.cpp TrajinList.cpp TrajoutList.cpp ReferenceList.cpp \
       CoordFileList.cpp PtrajActionList.cpp DataSet.cpp \
       dDataSet.cpp DataSetList.cpp vectormath.cpp FrameList.cpp Action_Rmsd.cpp \
       Action_Dihedral.cpp DataFile.cpp DataFileList.cpp PDBfile.cpp \
       PtrajFile.cpp PDBfileRoutines.cpp AtomMap.cpp BaseFileIO.cpp StdFile.cpp \
       GzipFile.cpp AmberRestart.cpp Action_Angle.cpp PtrajMpi.c MpiFile.cpp \
       Action_Strip.cpp Action.cpp RemdTraj.cpp AtomMask.cpp mapDataSet.cpp \
       Action_DSSP.cpp stringDataSet.cpp Action_Center.cpp Bzip2File.cpp \
       Action_Hbond.cpp intDataSet.cpp Action_Image.cpp Action_Surf.cpp \
       Action_Radgyr.cpp Conflib.cpp Action_Mask.cpp Action_Closest.cpp \
       NetcdfRoutines.cpp AmberRestartNC.cpp CpptrajStdio.cpp Mol2File.cpp \
       Mol2FileRoutines.cpp Action_NAstruct.cpp DistRoutines.cpp AxisType.cpp

OBJECTS=main.o AmberNetcdf.o AmberParm.o AmberTraj.o TrajFile.o \
        Frame.o PtrajState.o ArgList.o ptrajmask.o Action_Distance.o \
        ParmFileList.o TrajinList.o TrajoutList.o ReferenceList.o \
        CoordFileList.o PtrajActionList.o DataSet.o \
        dDataSet.o DataSetList.o vectormath.o FrameList.o Action_Rmsd.o \
        Action_Dihedral.o DataFile.o DataFileList.o PDBfile.o \
        PtrajFile.o PDBfileRoutines.o AtomMap.o BaseFileIO.o StdFile.o \
        GzipFile.o AmberRestart.o Action_Angle.o PtrajMpi.o MpiFile.o \
        Action_Strip.o Action.o RemdTraj.o AtomMask.o mapDataSet.o \
        Action_DSSP.o stringDataSet.o Action_Center.o Bzip2File.o \
        Action_Hbond.o intDataSet.o Action_Image.o Action_Surf.o \
        Action_Radgyr.o Conflib.o Action_Mask.o Action_Closest.o \
        NetcdfRoutines.o AmberRestartNC.o CpptrajStdio.o Mol2File.o \
        Mol2FileRoutines.o Action_NAstruct.o DistRoutines.o AxisType.o

HEADERS=AmberNetcdf.h AmberParm.h AmberTraj.h TrajFile.h \
        Frame.h PtrajState.h ArgList.h ptrajmask.h Action.h Action_Distance.h \
        ParmFileList.h TrajinList.h TrajoutList.h ReferenceList.h \
        CoordFileList.h PtrajActionList.h DataSet.h \
        dDataSet.h DataSetList.h vectormath.h FrameList.h Action_Rmsd.h \
        Action_Dihedral.h DataFile.h DataFileList.h PDBfile.h \
        PtrajFile.h PDBfileRoutines.h AtomMap.h BaseFileIO.h StdFile.h \
        GzipFile.h AmberRestart.h Action_Angle.h PtrajMpi.h MpiFile.h \
        Action_Strip.h RemdTraj.h AtomMask.h mapDataSet.h \
        Action_DSSP.h stringDataSet.h Action_Center.h Bzip2File.h \
        Action_Hbond.h intDataSet.h Action_Image.h Action_Surf.h \
        Action_Radgyr.h Conflib.h Action_Mask.h Action_Closest.h \
        NetcdfRoutines.h AmberRestartNC.h CpptrajStdio.h Mol2File.h \
        Mol2FileRoutines.h Action_NAstruct.h DistRoutines.h AxisType.h

all: cpptraj$(SFX)

install: cpptraj$(SFX) 
	mv cpptraj$(SFX) $(CPPTRAJBIN) 

findDepend: FindDepend.o 
	$(CXX) -Wall -o findDepend FindDepend.o

depend: findDepend
	./findDepend $(SOURCE) > cpptrajdepend

dependclean:
	-/bin/rm FindDepend.o
	-/bin/rm findDepend

cpptraj$(SFX): $(OBJECTS)
	$(CXX) $(LDFLAGS) -o cpptraj$(SFX) $(OBJECTS) $(NETCDFLIB) $(ZLIB) $(BZLIB)

clean:
	/bin/rm -f $(OBJECTS) cpptraj$(SFX)

# Header dependencies
include cpptrajdepend
