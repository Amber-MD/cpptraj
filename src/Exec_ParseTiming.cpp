#include "Exec_ParseTiming.h"
#include "CpptrajStdio.h"
#include "FileName.h"
#include "BufferedLine.h"
#include <cstring> // strncmp

/** CONSTRUCTOR */
Exec_ParseTiming::Exec_ParseTiming() :
  Exec(GENERAL)
{
  SetHidden( true );
}

/** Hold run timing and details. */
class Exec_ParseTiming::RunTiming {
  public:
    RunTiming() : version_major_(-1), version_minor_(-1), version_revision_(-1),
                  isMPI_(false), isOpenMP_(false), isCUDA_(false),
                  nprocs_(-1), t_total_(-1) {}
    RunTiming(std::string const& vstr, bool isM, bool isO, bool isC) :
      isMPI_(isM), isOpenMP_(isO), isCUDA_(isC), nprocs_(-1), t_total_(-1)
    {
      ArgList varg( vstr, "V." );
      version_major_ = varg.getNextInteger(-1);
      version_minor_ = varg.getNextInteger(-1);
      version_revision_ = varg.getNextInteger(-1);
    }

    void SetNprocs(int n) { nprocs_ = n; }
    void SetTotalTime(double t) { t_total_ = t; }

    void Print() const {
      mprintf("Version %i.%i.%i mpi=%i omp=%i cuda=%i nprocs=%i t_total=%g\n",
              version_major_, version_minor_, version_revision_,
              (int)isMPI_, (int)isOpenMP_, (int)isCUDA_, nprocs_, t_total_);
    }
  private:
    int version_major_;
    int version_minor_;
    int version_revision_;
    bool isMPI_;
    bool isOpenMP_;
    bool isCUDA_;
    int nprocs_;
    double t_total_;
};

int Exec_ParseTiming::read_cpptraj_output(std::string const& fname) {
  BufferedLine infile;

  if (infile.OpenFileRead( fname )) {
    mprinterr("Error: Could not open '%s'\n", fname.c_str());
    return 1;
  }
  RunTiming thisRun;
  const char* ptr = infile.Line();
  while (ptr != 0) {
    if (ptr[0] == 'C') {
      if (strncmp(ptr, "CPPTRAJ: Trajectory Analysis.", 29)==0) {
        mprintf("DEBUG: Title: %s\n", ptr);
        ArgList title( ptr+29, " " );
        mprintf("DEBUG: TitleArg: %s\n", title.ArgLine());
        // Version should be the first arg
        std::string versionStr = title.GetStringNext();
        bool isMPI = title.hasKey("MPI");
        bool isOpenMP = title.hasKey("OpenMP");
        bool isCUDA = title.hasKey("CUDA");
        thisRun = RunTiming(versionStr, isMPI, isOpenMP, isCUDA);
      }
    } else if (ptr[0] == '|') {
      if (strncmp(ptr, "| Running on", 12)==0) {
        ArgList procs( ptr+12, " " );
        // processes should be the next arg
        thisRun.SetNprocs( procs.getNextInteger(-1) );
      }
    } else if (ptr[0] == 'T') {
      if (strncmp(ptr, "TIME:", 5) == 0) {
        ArgList timeArg( ptr+5, " " );
        if (timeArg.Nargs() == 5 && timeArg[0] == "Total" && timeArg[1] == "execution")
          thisRun.SetTotalTime( timeArg.getKeyDouble("time:", -1) );
      } // END TIME:
    }
    ptr = infile.Line();
  }
  thisRun.Print();

  infile.CloseFile();
  return 0;
}
// Exec_ParseTiming::Help()
void Exec_ParseTiming::Help() const
{
  mprintf("\t<filename args> ...\n"
          "  Parse cpptraj output for timing data.\n"
         );
}


// Exec_ParseTiming::Execute()
Exec::RetType Exec_ParseTiming::Execute(CpptrajState& State, ArgList& argIn)
{
  File::NameArray FileNameList;

  std::string filearg = argIn.GetStringNext();
  while (!filearg.empty()) {
    File::NameArray fnames = File::ExpandToFilenames( filearg );
    if (fnames.empty()) {
      mprintf("Warning: No files matching '%s'\n", filearg.c_str());
    } else {
      for (File::NameArray::const_iterator it = fnames.begin(); it != fnames.end(); ++it)
        FileNameList.push_back( *it );
    }
    filearg = argIn.GetStringNext();
  }

  for (File::NameArray::const_iterator it = FileNameList.begin(); it != FileNameList.end(); ++it) {
    mprintf("\t%s\n", it->full());
    if (read_cpptraj_output( it->Full() )) {
      mprinterr("Error: Could not read cpptraj output '%s'\n", it->full());
      return CpptrajState::ERR;
    }
  }

  return CpptrajState::OK;
}
