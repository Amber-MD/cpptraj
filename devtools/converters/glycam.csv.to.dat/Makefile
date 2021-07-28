include ../../../config.h

CPPTRAJ_SOURCE=../../../src/ArgList.cpp ../../../src/StringRoutines.cpp ../../../src/CpptrajStdio.cpp

glycam.csv.to.dat: glycam.csv.to.dat.cpp ../../../src/ArgList.h $(CPPTRAJ_SOURCE)
	$(CXX) -o glycam.csv.to.dat glycam.csv.to.dat.cpp $(CPPTRAJ_SOURCE)

clean:
	-/bin/rm glycam.csv.to.dat
