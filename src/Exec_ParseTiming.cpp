#include "Exec_ParseTiming.h"
#include "CpptrajStdio.h"
#include "FileName.h"
#include "BufferedLine.h"
#include "DataSet_Mesh.h"
#include "StringRoutines.h"
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
      filename_(name), isMPI_(isM), isOpenMP_(isO), isCUDA_(isC), nprocs_(-1), nthreads_(-1), t_total_(-1)
    {
      ArgList varg( vstr, "V." );
      version_major_ = varg.getNextInteger(-1);
      version_minor_ = varg.getNextInteger(-1);
      version_revision_ = varg.getNextInteger(-1);
    }

    void SetNprocs(int n) { nprocs_ = n; }
    void SetNthreads(int n) { nthreads_ = n; }
    void SetTotalTime(double t) { t_total_ = t; }
    void AddTrajReadTime(double t) { t_traj_read_.push_back( t ); }
    void AddActionFrameTime(double t) { t_action_frame_.push_back( t ); }

    double TotalTime() const { return t_total_; }

    double TrajReadTime() const {
      double t = 0;
      for (Darray::const_iterator it = t_traj_read_.begin(); it != t_traj_read_.end(); ++it)
        if (*it > t) t = *it;
      return t;
    }

    double ActionFrameTime() const {
      double t = 0;
      for (Darray::const_iterator it = t_action_frame_.begin(); it != t_action_frame_.end(); ++it)
        if (*it > t) t = *it;
      return t;
    }

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

    bool HasDetailedTiming() const { return !(t_traj_read_.empty() || t_action_frame_.empty()); }

    FileName const& Filename() const { return filename_; }

    std::string Name() const {
      std::string out;
      if (isMPI_ && isOpenMP_)
        out.assign("H" + integerToString(nprocs_) + "x" + integerToString(nthreads_));
      else if (isMPI_)
        out.assign("M" + integerToString(nprocs_));
      else if (isOpenMP_)
        out.assign("O" + integerToString(nthreads_));
      else
        out.assign("S");
      if (isCUDA_) out.append("(C)");
      return out;
    }

    void Print() const {
      mprintf("%s Version %i.%i.%i mpi=%i omp=%i cuda=%i nprocs=%i nthreads=%i ncores=%i t_total=%g\n",
              filename_.full(), version_major_, version_minor_, version_revision_,
              (int)isMPI_, (int)isOpenMP_, (int)isCUDA_, nprocs_, nthreads_, TotalCores(), t_total_);
      for (Darray::const_iterator it = t_traj_read_.begin(); it != t_traj_read_.end(); ++it)
        mprintf(" %g", *it);
      for (Darray::const_iterator it = t_action_frame_.begin(); it != t_action_frame_.end(); ++it)
        mprintf(" %g", *it);
      mprintf("\n");
    }

    bool IsBad() {
      errMsg_.clear();
      if (filename_.empty()) { errMsg_.assign("Empty run name.\n"); return true; }
      if (t_total_ < 0) { errMsg_.assign("Empty run time.\n"); return true; }
      if (isMPI_ && nprocs_ < 1) { errMsg_.assign("MPI && procs < 1.\n"); return true; }
      if (isOpenMP_ && nthreads_ < 1) { errMsg_.assign("OpenMP && threads < 1.\n"); return true; }
      return false;
    }

    /** \return last error message set by IsBad() */
    std::string const& ErrMsg() const { return errMsg_; }

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
    FileName filename_;
    int version_major_;
    int version_minor_;
    int version_revision_;
    bool isMPI_;
    bool isOpenMP_;
    bool isCUDA_;
    int nprocs_;
    int nthreads_;
    double t_total_;
    Darray t_traj_read_;
    Darray t_action_frame_;
    std::string errMsg_;
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
        //mprintf("DEBUG: OpenMP threads. %s\n", ptr);
        // threads should be the next integer arg
        thisRun.SetNthreads( infoArg.getNextInteger(-1) );
      }
    } else if (ptr[0] == 'T') {
      if (strncmp(ptr, "TIME:", 5) == 0) {
        ArgList timeArg( ptr, " \t" );
        //timeArg.PrintDebug();
        if (timeArg.Nargs() == 6 && timeArg[1] == "Total" && timeArg[2] == "execution")
          thisRun.SetTotalTime( timeArg.getKeyDouble("time:", -1) );
        else if (timeArg.Nargs() > 3 && timeArg[1] == "Trajectory" && timeArg[2] == "read:")
          thisRun.AddTrajReadTime( timeArg.getNextDouble(-1) );
        else if (timeArg.Nargs() > 4 && timeArg[1] == "Action" && timeArg[2] == "frame" && timeArg[3] == "processing:")
          thisRun.AddActionFrameTime( timeArg.getNextDouble(-1) );
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
  mprintf("\t<filename args> ... [sortby {time|cores}] [out <file>] [name <setname>]\n"
          "\t[type {trajread|actframe}]\n"
          "  Parse cpptraj output for timing data.\n"
         );
}


// Exec_ParseTiming::Execute()
Exec::RetType Exec_ParseTiming::Execute(CpptrajState& State, ArgList& argIn)
{
  enum SortType { SORT_T_TOTAL = 0, SORT_CORES };
  const char* SortTypeStr[] = {"time", "total # cores"};
  SortType sort = SORT_T_TOTAL;
  std::string sortarg = argIn.GetStringKey("sortby");
  if (!sortarg.empty()) {
    if (sortarg == "time")
      sort = SORT_T_TOTAL;
    else if (sortarg == "cores")
      sort = SORT_CORES;
    else {
      mprinterr("Error: Unrecognized sort: %s\n", sortarg.c_str());
      return CpptrajState::ERR;
    }
  }
  mprintf("\tSort by %s\n", SortTypeStr[sort]);

  DataFile* outfile = State.DFL().AddDataFile( argIn.GetStringKey("out"), argIn );
  if (outfile != 0)
    mprintf("\tOutput to file '%s'\n", outfile->DataFilename().full());

  std::string dsname = argIn.GetStringKey("name");
  if (dsname.empty())
    dsname = State.DSL().GenerateDefaultName("TIMING");
  mprintf("\tDataSet name: %s\n", dsname.c_str());

  // Which variables to plot
  enum Xtype { X_INDEX=0, X_CORES };
  enum Ytype { Y_T_TOTAL=0, Y_T_TRAJREAD, Y_T_ACTFRAME };

  Xtype xvar;
  Ytype yvar = Y_T_TOTAL;
  Dimension Xdim;
  //Dimension Ydim;

  if (sort == SORT_CORES) {
    xvar = X_CORES;
    Xdim.SetLabel("Cores");
  } else {
    xvar = X_INDEX;
    Xdim.SetLabel("Run");
  }

  bool needsDetailedTiming = false;
  std::string typearg = argIn.GetStringKey("type");
  if (typearg.empty())
    yvar = Y_T_TOTAL;
  else if (typearg == "trajread") {
    yvar = Y_T_TRAJREAD;
    needsDetailedTiming = true;
  } else if (typearg == "actframe") {
    yvar = Y_T_ACTFRAME;
    needsDetailedTiming = true;
  }
  //Ydim.SetLabel("TotalTime");

  // Only file name args below here
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
    RunTiming thisRun = read_cpptraj_output( it->Full() );
    if (thisRun.IsBad()) {
      mprintf("Warning: Problem reading cpptraj output from '%s' : %s\n", it->full(), thisRun.ErrMsg().c_str());
      //return CpptrajState::ERR;
    } else {
      Runs.push_back( thisRun );
    }
    //Runs.back().Print();
    //mprintf("\n");
  }

  switch (sort) {
    case SORT_T_TOTAL : std::sort(Runs.begin(), Runs.end(), RunTiming::sort_by_total_time()); break;
    case SORT_CORES   : std::sort(Runs.begin(), Runs.end(), RunTiming::sort_by_cores()); break;
  }

  DataSet* ds = State.DSL().AddSet( DataSet::XYMESH, MetaData(dsname) );
  if (ds == 0) {
    mprinterr("Error: Could not allocate output set.\n");
    return CpptrajState::ERR;
  }
  if (outfile != 0) outfile->AddDataSet( ds );
  DataSet_Mesh& outset = static_cast<DataSet_Mesh&>( *ds );
  outset.Allocate( DataSet::SizeArray(1, Runs.size()) );
  outset.SetDim(Dimension::X, Xdim);
  //outset.SetDim(Dimension::Y, Ydim);

  DataSet* nameSet = State.DSL().AddSet(DataSet::STRING, MetaData(dsname, "name"));
  if (nameSet == 0) {
    mprinterr("Error: Could not allocate output names set.\n");
    return CpptrajState::ERR;
  }
  if (outfile != 0) outfile->AddDataSet( nameSet );
  nameSet->Allocate( DataSet::SizeArray(1, Runs.size()) );
  nameSet->SetDim(Dimension::X, Xdim);

  DataSet* dirNameSet = State.DSL().AddSet(DataSet::STRING, MetaData(dsname, "dir"));
  if (dirNameSet == 0) {
    mprinterr("Error: Could not allocate directory names set.\n");
    return CpptrajState::ERR;
  }
  if (outfile != 0) outfile->AddDataSet( dirNameSet );
  dirNameSet->Allocate( DataSet::SizeArray(1, Runs.size()) );
  dirNameSet->SetDim(Dimension::X, Xdim);

  for (RunArray::const_iterator it = Runs.begin(); it != Runs.end(); ++it) {
    if (needsDetailedTiming && !it->HasDetailedTiming()) {
      mprinterr("Error: Detailed timing requested but output %s does not contain detailed timing.\n",
                it->Filename().full());
      it->Print();
      return CpptrajState::ERR;
    }
    double X = 0;
    switch (xvar) {
      case X_INDEX : X = (double)(it - Runs.begin()); break;
      case X_CORES : X = (double)it->TotalCores(); break;
    }
    double Y = 0;
    switch (yvar) {
      case Y_T_TOTAL : Y = it->TotalTime(); break;
      case Y_T_TRAJREAD : Y = it->TrajReadTime(); break;
      case Y_T_ACTFRAME : Y = it->ActionFrameTime(); break;
    }
    outset.AddXY(X, Y);
    nameSet->Add(it - Runs.begin(), it->Name().c_str());
    dirNameSet->Add(it - Runs.begin(), it->Filename().DirPrefix().c_str());
    it->Print();
  }

  return CpptrajState::OK;
}
