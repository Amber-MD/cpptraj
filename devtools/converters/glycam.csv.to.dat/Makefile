include ../../../config.h

all: glycam.csv.to.dat ReadComponents

clean:
	-/bin/rm glycam.csv.to.dat ReadComponents 

# -----------------------------------------------
CPPTRAJ_SOURCE=../../../src/ArgList.cpp ../../../src/StringRoutines.cpp ../../../src/CpptrajStdio.cpp

glycam.csv.to.dat: glycam.csv.to.dat.cpp ../../../src/ArgList.h ../../../src/StringRoutines.h $(CPPTRAJ_SOURCE)
	$(CXX) -o glycam.csv.to.dat glycam.csv.to.dat.cpp $(CPPTRAJ_SOURCE)

# -----------------------------------------------
#SOURCES=../../src/ArgList.cpp CIFfile.cpp \
#        FileIO_Bzip2.cpp FileIO_MpiShared.cpp \
#        StringRoutines.cpp Atom.cpp CpptrajFile.cpp \
#        FileIO_Gzip.cpp FileIO_Std.cpp BufferedLine.cpp \
#        CpptrajStdio.cpp FileIO_Mpi.cpp ReadComponents.cpp \
#        FileName.cpp NameType.cpp
#OBJECTS=$(SOURCES:.cpp=.o)
LIBS=$(CPPTRAJLIB)/libcpptraj_file.a $(CPPTRAJLIB)/libcpptraj_core.a

ReadComponents: ReadComponents.cpp $(OBJECTS) libs
	$(CXX) -I../../../src -o ReadComponents ReadComponents.cpp GlycamPdbResMap.cpp $(LIBS) $(LDFLAGS)

libs:
	cd ../../../src && $(MAKE) libstatic 

