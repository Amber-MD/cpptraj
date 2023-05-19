#include "Exec_ParseTiming.h"
#include "CpptrajStdio.h"
#include "FileName.h"
#include "BufferedLine.h"
#include <cstring> // strncmp
#include <algorithm> // std::sort

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
                  nprocs_(-1), nthreads_(-1), t_total_(-1) {}
    RunTiming(std::string const& name, std::string const& vstr, bool isM, bool isO, bool isC) :
      name_(name), isMPI_(isM), isOpenMP_(isO), isCUDA_(isC), nprocs_(-1), nthreads_(-1), t_total_(-1)
    {
      ArgList varg( vstr, "V." );
      version_major_ = varg.getNextInteger(-1);
      version_minor_ = varg.getNextInteger(-1);
      version_revision_ = varg.getNextInteger(-1);
    }

    void SetNprocs(int n) { nprocs_ = n; }
    void SetNthreads(int n) { nthreads_ = n; }
    void SetTotalTime(double t) { t_total_ = t; }

    double TotalTime() const { return t_total_; }

    int TotalCores() const {
      if (isMPI_ && isOpenMP_)
        return nprocs_ * nthreads_;
      else if (isMPI_)
        return nprocs_;
      else if (isOpenMP_)
        return nthreads_;
      else
        return 1;
    }

    void Print() const {
      mprintf("%s Version %i.%i.%i mpi=%i omp=%i cuda=%i nprocs=%i nthreads=%i ncores=%i t_total=%g\n",
              name_.c_str(), version_major_, version_minor_, version_revision_,
              (int)isMPI_, (int)isOpenMP_, (int)isCUDA_, nprocs_, nthreads_, TotalCores(), t_total_);
    }

    bool IsBad() const {
      if (name_.empty()) { mprinterr("Error: Empty run name.\n"); return true; }
      if (t_total_ < 0) { mprinterr("Error: Empty run time.\n"); return true; }
      if (isMPI_ && nprocs_ < 1) { mprinterr("Error: MPI && procs < 1.\n"); return true; }
      if (isOpenMP_ && nthreads_ < 1) { mprinterr("Error: OpenMP && threads < 1.\n"); return true; }
      return false;
    }

    /** Sort by total time, longest first. */
    struct sort_by_total_time {
      inline bool operator()(RunTiming const& first, RunTiming const& second) const {
        return (first.TotalTime() > second.TotalTime());
      }
    };

    /** Sort by total cores, then time. */
    struct sort_by_cores {
      inline bool operator()(RunTiming const& first, RunTiming const& second) const {
        int first_cores = first.TotalCores();
        int second_cores = second.TotalCores();
        if (first_cores == second_cores)
          return (first.TotalTime() > second.TotalTime());
        else
          return (first_cores < second_cores);
      }
    };
  private:
    std::string name_;
    int version_major_;
    int version_minor_;
    int version_revision_;
    bool isMPI_;
    bool isOpenMP_;
    bool isCUDA_;
    int nprocs_;
    int nthreads_;
    double t_total_;
};

/** Execute command. */
Exec_ParseTiming::RunTiming Exec_ParseTiming::read_cpptraj_output(std::string const& fname) {
  BufferedLine infile;


  RunTiming thisRun;
  if (infile.OpenFileRead( fname )) {
    mprinterr("Error: Could not open '%s'\n", fname.c_str());
    return thisRun;
  }
  const char* ptr = infile.Line();
  while (ptr != 0) {
    if (ptr[0] == 'C') {
      if (strncmp(ptr, "CPPTRAJ: Trajectory Analysis.", 29)==0) {
        //mprintf("DEBUG: Title: %s\n", ptr);
        ArgList titleArg( ptr+29, " " );
        //mprintf("DEBUG: TitleArg: %s\n", title.ArgLine());
        // Version should be the first arg
        std::string versionStr = titleArg.GetStringNext();
        bool isMPI = titleArg.hasKey("MPI");
        bool isOpenMP = titleArg.hasKey("OpenMP");
        bool isCUDA = titleArg.hasKey("CUDA");
        thisRun = RunTiming(fname, versionStr, isMPI, isOpenMP, isCUDA);
      }
    } else if (ptr[0] == '|') {
      ArgList infoArg(ptr+1, " ");
      if (infoArg.Nargs() == 4 && infoArg[0] == "Running" && infoArg[1] == "on") {
        // processes should be the next integer arg
        thisRun.SetNprocs( infoArg.getNextInteger(-1) );
      } else if (infoArg.Nargs() == 4 && infoArg[1] == "OpenMP" &&  infoArg[2] == "threads") {
        mprintf("DEBUG: OpenMP threads. %s\n", ptr);
        // threads should be the next integer arg
        thisRun.SetNthreads( infoArg.getNextInteger(-1) );
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
  //thisRun.Print();
  infile.CloseFile();
  return thisRun;
}
// Exec_ParseTiming::Help()
void Exec_ParseTiming::Help() const
{
  mprintf("\t<filename args> ... [sortby {time|cores}]\n"
          "  Parse cpptraj output for timing data.\n"
         );
}


// Exec_ParseTiming::Execute()
Exec::RetType Exec_ParseTiming::Execute(CpptrajState& State, ArgList& argIn)
{
  enum SortType { SORT_T_TOTAL = 0, SORT_CORES };
  const char* SortTypeStr[] = {"time", "total # cores"};
  SortType sort;
  std::string sortarg = argIn.GetStringKey("sortby");
  if (sortarg.empty())
    sort = SORT_T_TOTAL;
  else if (sortarg == "time")
    sort = SORT_T_TOTAL;
  else if (sortarg == "cores")
    sort = SORT_CORES;
  mprintf("\tSort by %s\n", SortTypeStr[sort]);

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

  typedef std::vector<RunTiming> RunArray;
  RunArray Runs;
  Runs.reserve( FileNameList.size() );
  for (File::NameArray::const_iterator it = FileNameList.begin(); it != FileNameList.end(); ++it) {
    //mprintf("\t%s\n", it->full());
    Runs.push_back( read_cpptraj_output( it->Full() ) );
    if (Runs.back().IsBad()) {
      mprinterr("Error: Problem reading cpptraj output from '%s'\n", it->full());
      return CpptrajState::ERR;
    }
    //Runs.back().Print();
    //mprintf("\n");
  }

  switch (sort) {
    case SORT_T_TOTAL : std::sort(Runs.begin(), Runs.end(), RunTiming::sort_by_total_time()); break;
    case SORT_CORES   : std::sort(Runs.begin(), Runs.end(), RunTiming::sort_by_cores()); break;
  }

  for (RunArray::const_iterator it = Runs.begin(); it != Runs.end(); ++it) {
    it->Print();
    mprintf("\n");
  }

  return CpptrajState::OK;
}
