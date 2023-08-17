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

// -----------------------------------------------------------------------------
/** Hold run timing and details. */
class Exec_ParseTiming::RunTiming {
  public:
    /// CONSTRUCTOR - blank
    RunTiming() : version_major_(-1), version_minor_(-1), version_revision_(-1),
                  isMPI_(false), isOpenMP_(false), isCUDA_(false),
                  nprocs_(-1), nthreads_(-1), t_total_(-1), t_traj_proc_(-1) {}
    /// CONSTRUCTOR - filename
    RunTiming(std::string const& name) : filename_(name),
                  version_major_(-1), version_minor_(-1), version_revision_(-1),
                  isMPI_(false), isOpenMP_(false), isCUDA_(false),
                  nprocs_(-1), nthreads_(-1), t_total_(-1), t_traj_proc_(-1)
    {
      // Get directory prefix.
      std::string dirprefix( filename_.DirPrefix_NoSlash() );
      FileName fname( filename_ );
      while (!dirprefix.empty()) {
        fname.SetFileName_NoExpansion( fname.DirPrefix_NoSlash() );
        dirprefix = fname.DirPrefix_NoSlash();
      }
      prefix_ = fname.Base();
      //FileName dirprefix2( dirprefix1.Full() );
      //mprintf("DEBUG: '%s' Dirprefix= '%s' name= '%s'\n", filename_.full(), fname.full(), fname.base());
    }

    void SetVersionAndType(std::string const& vstr, bool isM, bool isO, bool isC) {
      isMPI_ = isM;
      isOpenMP_ = isO;
      isCUDA_ = isC;
      ArgList varg( vstr, "V." );
      version_major_ = varg.getNextInteger(-1);
      version_minor_ = varg.getNextInteger(-1);
      version_revision_ = varg.getNextInteger(-1);
    }

    void SetNprocs(int n) { nprocs_ = n; }
    void SetNthreads(int n) { nthreads_ = n; }
    void SetTotalTime(double t) { t_total_ = t; }
    void SetTrajProcTime(double t) { t_traj_proc_ = t; }
    void AddTrajReadTime(double t) { t_traj_read_.push_back( t ); }
    void AddActionFrameTime(double t) { t_action_frame_.push_back( t ); }

    int Nthreads() const { return nthreads_; }
    double TotalTime() const { return t_total_; }
    double TrajProcTime() const { return t_traj_proc_; }

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
    std::string const& Prefix() const { return prefix_; }

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
      if (isCUDA_) {
        std::string tmp = out;
        out.assign( "G(" + tmp + ")" );
      }
      return out;
    }

    void Print() const {
      mprintf("%s Version %i.%i.%i mpi=%i omp=%i cuda=%i nprocs=%i nthreads=%i ncores=%i t_total=%g t_traj_proc=%g",
              filename_.full(), version_major_, version_minor_, version_revision_,
              (int)isMPI_, (int)isOpenMP_, (int)isCUDA_, nprocs_, nthreads_, TotalCores(), t_total_, t_traj_proc_);
      if (!t_traj_read_.empty()) {
        mprintf(" t_traj_read={");
        for (Darray::const_iterator it = t_traj_read_.begin(); it != t_traj_read_.end(); ++it)
          mprintf(" %g", *it);
        mprintf("}");
      }
      if (!t_action_frame_.empty()) {
        mprintf(" t_action_frame={");
        for (Darray::const_iterator it = t_action_frame_.begin(); it != t_action_frame_.end(); ++it)
          mprintf(" %g", *it);
        mprintf("}");
      }
      mprintf("\n");
    }

    bool IsBad() {
      errMsg_.clear();
      if (filename_.empty()) { errMsg_.assign("Empty run name."); return true; }
      if (t_total_ < 0) { errMsg_.assign("Empty run time."); return true; }
      if (isMPI_ && nprocs_ < 1) { errMsg_.assign("MPI && procs < 1."); return true; }
      if (isOpenMP_ && nthreads_ < 1) { errMsg_.assign("OpenMP && threads < 1."); return true; }
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

    /** Sort by file name. */
    struct sort_by_filename {
      inline bool operator()(RunTiming const& first, RunTiming const& second) const {
        if (first.filename_.Full() == second.filename_.Full())
          return (first.TotalTime() > second.TotalTime());
        else
          return (first.filename_.Full() < second.filename_.Full());
      }
    };
  private:
    FileName filename_;
    std::string prefix_; ///< Filename prefix
    int version_major_;
    int version_minor_;
    int version_revision_;
    bool isMPI_;
    bool isOpenMP_;
    bool isCUDA_;
    int nprocs_;
    int nthreads_;
    double t_total_;
    double t_traj_proc_; 
    Darray t_traj_read_;
    Darray t_action_frame_;
    std::string errMsg_;
};
// -----------------------------------------------------------------------------

/** Execute command. */
Exec_ParseTiming::RunTiming Exec_ParseTiming::read_cpptraj_output(std::string const& fname) {
  BufferedLine infile;


  RunTiming thisRun( fname );
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
        thisRun.SetVersionAndType(versionStr, isMPI, isOpenMP, isCUDA);
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
        else if (timeArg.Nargs() >= 6 && timeArg[1] == "Trajectory" && timeArg[2] == "Process")
          thisRun.SetTrajProcTime( timeArg.getKeyDouble(":", -1) );
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

/** \return Y value according to yvar */
double Exec_ParseTiming::YVAL(Ytype yvar, RunTiming const& run)
{
  double Y = 0;
  switch (yvar) {
    case Y_T_TOTAL    : Y = run.TotalTime(); break;
    case Y_T_TRAJPROC : Y = run.TrajProcTime(); break;
    case Y_T_TRAJREAD : Y = run.TrajReadTime(); break;
    case Y_T_ACTFRAME : Y = run.ActionFrameTime(); break;
  }
  return Y;
}

/** Make an output data set for RunArray */
int Exec_ParseTiming::create_output_set(RunArray const& Runs, DataSetList& DSL,
                                        DataFile* outfile, std::string const& dsname,
                                        Dimension const& Xdim, Xtype xvar, Ytype yvar)
const
{
  // Create DataSet

  DataSet* ds = DSL.AddSet( DataSet::XYMESH, MetaData(dsname) );
  if (ds == 0) {
    mprinterr("Error: Could not allocate output set.\n");
    return 1;
  }
  if (outfile != 0) outfile->AddDataSet( ds );
  DataSet_Mesh& outset = static_cast<DataSet_Mesh&>( *ds );
  outset.Allocate( DataSet::SizeArray(1, Runs.size()) );
  outset.SetDim(Dimension::X, Xdim);
  //outset.SetDim(Dimension::Y, Ydim);

  DataSet* nameSet = DSL.AddSet(DataSet::STRING, MetaData(dsname, "name"));
  if (nameSet == 0) {
    mprinterr("Error: Could not allocate output names set.\n");
    return 1;
  }
  if (outfile != 0) outfile->AddDataSet( nameSet );
  nameSet->Allocate( DataSet::SizeArray(1, Runs.size()) );
  nameSet->SetDim(Dimension::X, Xdim);

  DataSet* dirNameSet = DSL.AddSet(DataSet::STRING, MetaData(dsname, "dir"));
  if (dirNameSet == 0) {
    mprinterr("Error: Could not allocate directory names set.\n");
    return 1;
  }
  if (outfile != 0) outfile->AddDataSet( dirNameSet );
  dirNameSet->Allocate( DataSet::SizeArray(1, Runs.size()) );
  dirNameSet->SetDim(Dimension::X, Xdim);

  for (RunArray::const_iterator it = Runs.begin(); it != Runs.end(); ++it) {
    double X = 0;
    switch (xvar) {
      case X_INDEX : X = (double)(it - Runs.begin()); break;
      case X_CORES : X = (double)it->TotalCores(); break;
    }
    double Y = YVAL( yvar, *it );

    outset.AddXY(X, Y);
    nameSet->Add(it - Runs.begin(), it->Name().c_str());
    dirNameSet->Add(it - Runs.begin(), it->Filename().DirPrefix().c_str());
    //it->Print();
  }
  return 0;
}

/** Write a run group to file. */
void Exec_ParseTiming::write_to_file(CpptrajFile& outfile, RunArray const& Runs, Xtype xvar, Ytype yvar,
                                     double refy, double refcores)
const
{
  //double refy = 0;
  //double refcores = 0;
  for (RunArray::const_iterator it = Runs.begin(); it != Runs.end(); ++it) {
    double X = 0;
    switch (xvar) {
      case X_INDEX : X = (double)(it - Runs.begin()); break;
      case X_CORES : X = (double)it->TotalCores(); break;
    }
    double Y = YVAL( yvar, *it );
    
    //if (it == Runs.begin()) {
    //  refy = Y;
    //  refcores = (double)it->TotalCores();
    //}
    double speedup = refy / Y;
    int totalCores = it->TotalCores();
    double ideal_speedup = (double)totalCores / refcores;
    double efficiency = speedup / ideal_speedup;
    outfile.Printf("%12.4f %12.4f %6.2f %6i %6.2f %12s %s\n", X, Y, speedup, totalCores, efficiency, it->Name().c_str(), it->Filename().DirPrefix().c_str());
  }
}

// Exec_ParseTiming::Help()
void Exec_ParseTiming::Help() const
{
  mprintf("\t<filename args> ... [out <file>] [name <setname>]\n"
          "\t[sortby {time|cores|filename}] [includebad] [showdetails]\n"
          "\t[type {trajproc|trajread|actframe}] [reverse]\n"
          "\t[groupout <file> [grouptype {prefix|name|kind}]]\n"
          "  Parse cpptraj output for timing data.\n"
         );
}

// Exec_ParseTiming::Execute()
Exec::RetType Exec_ParseTiming::Execute(CpptrajState& State, ArgList& argIn)
{
  enum SortType { SORT_T_TOTAL = 0, SORT_CORES, SORT_FILENAME };
  const char* SortTypeStr[] = {"time", "total # cores", "file name"};
  SortType sort = SORT_T_TOTAL;
  std::string sortarg = argIn.GetStringKey("sortby");
  if (!sortarg.empty()) {
    if (sortarg == "time")
      sort = SORT_T_TOTAL;
    else if (sortarg == "cores")
      sort = SORT_CORES;
    else if (sortarg == "filename")
      sort = SORT_FILENAME;
    else {
      mprinterr("Error: Unrecognized sort: %s\n", sortarg.c_str());
      return CpptrajState::ERR;
    }
  }
  mprintf("\tSort by %s.\n", SortTypeStr[sort]);

  bool reverse_sort = argIn.hasKey("reverse");
  if (reverse_sort)
    mprintf("\tPerforming reverse sort.\n");

  bool showdetails = argIn.hasKey("showdetails");
  if (showdetails)
    mprintf("\tWill print details for each run to STDOUT.\n");

  bool include_bad = argIn.hasKey("includebad");
  if (include_bad)
    mprintf("\tIncluding bad results.\n");

  DataFile* outfile = State.DFL().AddDataFile( argIn.GetStringKey("out"), argIn );
  if (outfile != 0)
    mprintf("\tOutput to file '%s'\n", outfile->DataFilename().full());

  std::string groupOutArg = argIn.GetStringKey("groupout");
  CpptrajFile* groupout = 0;
  GroupType groupOpt = GROUPBY_PREFIX;
  if (!groupOutArg.empty()) {
    groupout = State.DFL().AddCpptrajFile(groupOutArg, "Timing by group", DataFileList::TEXT);
    if (groupout == 0) {
      mprinterr("Error: Could not add file '%s'\n", groupOutArg.c_str());
      return CpptrajState::ERR;
    }
    mprintf("\tGroup output to '%s'\n", groupout->Filename().full());
    // NOTE: groupby is already a keyword for standard data file output
    std::string grouptypeArg = argIn.GetStringKey("grouptype");
    if (!grouptypeArg.empty()) {
      if (grouptypeArg == "prefix")
        groupOpt = GROUPBY_PREFIX;
      else if (grouptypeArg == "name")
        groupOpt = GROUPBY_NAME;
      else if (grouptypeArg == "kind")
        groupOpt = GROUPBY_KIND;
      else {
        mprinterr("Error: Unrecognized option for 'grouptype': %s\n", grouptypeArg.c_str());
        return CpptrajState::ERR;
      }
    }
    static const char* GroupTypeStr[] = { "directory prefix", "run type name", "run kind" };
    mprintf("\tGroup by %s.\n", GroupTypeStr[groupOpt]);
  }

  std::string dsname = argIn.GetStringKey("name");
  if (dsname.empty())
    dsname = State.DSL().GenerateDefaultName("TIMING");
  mprintf("\tDataSet name: %s\n", dsname.c_str());

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
  const char* YtypeStr[] = { "Total", "Trajectory Process", "Trajectory Read", "Action Frame" };
  std::string typearg = argIn.GetStringKey("type");
  if (typearg.empty())
    yvar = Y_T_TOTAL;
  else if (typearg == "trajproc")
    yvar = Y_T_TRAJPROC;
  else if (typearg == "trajread") {
    yvar = Y_T_TRAJREAD;
    needsDetailedTiming = true;
  } else if (typearg == "actframe") {
    yvar = Y_T_ACTFRAME;
    needsDetailedTiming = true;
  }
  mprintf("\tUsing %s time.\n", YtypeStr[yvar]);
  //Ydim.SetLabel("TotalTime");

  // ----- Only file name args below here --------
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
      if (include_bad)
        Runs.push_back( thisRun );
      //return CpptrajState::ERR;
    } else {
      Runs.push_back( thisRun );
    }
    //Runs.back().Print();
    //mprintf("\n");
  }
  mprintf("\t%zu runs.\n", Runs.size());

  // Check Runs
  for (RunArray::const_iterator it = Runs.begin(); it != Runs.end(); ++it) {
    if (needsDetailedTiming && !it->HasDetailedTiming()) {
      mprinterr("Error: Detailed timing requested but output %s does not contain detailed timing.\n",
                it->Filename().full());
      //it->Print();
      return CpptrajState::ERR;
    }
  }

  // Sort runs
  switch (sort) {
    case SORT_T_TOTAL  : std::sort(Runs.begin(), Runs.end(), RunTiming::sort_by_total_time()); break;
    case SORT_CORES    : std::sort(Runs.begin(), Runs.end(), RunTiming::sort_by_cores()); break;
    case SORT_FILENAME : std::sort(Runs.begin(), Runs.end(), RunTiming::sort_by_filename()); break;
  }

  if (reverse_sort)
    std::reverse( Runs.begin(), Runs.end() );

  // Group runs - TODO use pointers?
  if (groupout != 0) {
    //groupout.OpenWrite("");
    typedef std::map<std::string, RunArray> RunMap;
    RunMap groupByPrefix;
    // Reference Y for speedup will be the largest Y
    double refY = 0;
    double refCores = 0;
    for (RunArray::const_iterator it = Runs.begin(); it != Runs.end(); ++it)
    {
      // Determine reference
      double runY = YVAL(yvar, *it);
      if (runY > refY) {
        refY = runY;
        refCores = (double)it->TotalCores();
      }
      // Insert into map by desired group type
      std::string key;
      switch (groupOpt) {
        case GROUPBY_PREFIX : key = it->Prefix(); break;
        case GROUPBY_NAME   : key = it->Name(); break;
        case GROUPBY_KIND   :
          key = it->Name()[0];
          if (key[0] == 'H')
            key.append(integerToString(it->Nthreads())); 
          break;
      }
      RunMap::iterator jt = groupByPrefix.lower_bound( key );
      if (jt == groupByPrefix.end() || jt->first != key )
      {
        // New group
        jt = groupByPrefix.insert( jt, std::pair<std::string,RunArray>( key, RunArray() ) );
      }
      jt->second.push_back( *it );
    }
    mprintf("\tFound %zu groups.\n", groupByPrefix.size());
    for (RunMap::const_iterator jt = groupByPrefix.begin(); jt != groupByPrefix.end(); ++jt)
    {
      groupout->Printf("#\t\t%s\n", jt->first.c_str());
      write_to_file(*groupout, jt->second, xvar, yvar, refY, refCores);
      //for (RunArray::const_iterator it = jt->second.begin(); it != jt->second.end(); ++it)
      //  mprintf("\t\t\t%s\n", it->Filename().DirPrefix().c_str());
    }
  }

  // Create DataSet
  if (create_output_set(Runs, State.DSL(), outfile, dsname, Xdim, xvar, yvar))
    return CpptrajState::ERR;
  // Print details
  if (showdetails) {
    mprintf("\tRun details:\n");
    for (RunArray::const_iterator it = Runs.begin(); it != Runs.end(); ++it) {
      mprintf("\t\t");
      it->Print();
    }
  }
  return CpptrajState::OK;
}
