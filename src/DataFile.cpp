#include "DataFile.h"
#include "CpptrajStdio.h"
#ifdef TIMER
# include "Timer.h"
#endif
#ifdef MPI
# include "Parallel.h"
#endif
// All DataIO classes go here
#include "DataIO_Std.h"
#include "DataIO_Grace.h"
#include "DataIO_Gnuplot.h"
#include "DataIO_Xplor.h"
#include "DataIO_OpenDx.h"
#include "DataIO_RemLog.h"
#include "DataIO_Mdout.h"
#include "DataIO_Evecs.h"
#include "DataIO_VecTraj.h"
#include "DataIO_XVG.h"

// CONSTRUCTOR
DataFile::DataFile() :
  debug_(0),
  member_(-1), 
  dimension_(-1),
  dfType_(DATAFILE),
  dflWrite_(true),
  setDataSetPrecision_(false), //TODO: Just use default_width_ > -1?
  default_width_(-1),
  default_precision_(-1),
  dataio_(0),
  defaultDim_(3), // default to X/Y/Z dims
  minIsSet_(3, false)
{}

// DESTRUCTOR
DataFile::~DataFile() { if (dataio_ != 0) delete dataio_; }

// ----- STATIC VARS / ROUTINES ------------------------------------------------
// NOTE: Must be in same order as DataFormatType
const FileTypes::AllocToken DataFile::DF_AllocArray[] = {
  { "Standard Data File", DataIO_Std::ReadHelp,    DataIO_Std::WriteHelp,    DataIO_Std::Alloc    },
  { "Grace File",         0,                       DataIO_Grace::WriteHelp,  DataIO_Grace::Alloc  },
  { "Gnuplot File",       0,                       DataIO_Gnuplot::WriteHelp,DataIO_Gnuplot::Alloc},
  { "Xplor File",         0,                       0,                        DataIO_Xplor::Alloc  },
  { "OpenDX File",        0,                       DataIO_OpenDx::WriteHelp, DataIO_OpenDx::Alloc },
  { "Amber REM log",      DataIO_RemLog::ReadHelp, 0,                        DataIO_RemLog::Alloc },
  { "Amber MDOUT file",   0,                       0,                        DataIO_Mdout::Alloc  },
  { "Evecs file",         DataIO_Evecs::ReadHelp,  0,                        DataIO_Evecs::Alloc  },
  { "Vector pseudo-traj", 0,                       DataIO_VecTraj::WriteHelp,DataIO_VecTraj::Alloc},
  { "XVG file",           0,                       0,                        DataIO_XVG::Alloc    },
  { "Unknown Data file",  0,                       0,                        0                    }
};

const FileTypes::KeyToken DataFile::DF_KeyArray[] = {
  { DATAFILE,     "dat",    ".dat"   },
  { XMGRACE,      "grace",  ".agr"   },
  { XMGRACE,      "grace",  ".xmgr"  },
  { GNUPLOT,      "gnu",    ".gnu"   },
  { XPLOR,        "xplor",  ".xplor" },
  { XPLOR,        "xplor",  ".grid"  },
  { OPENDX,       "opendx", ".dx"    },
  { REMLOG,       "remlog", ".log"   },
  { MDOUT,        "mdout",  ".mdout" },
  { EVECS,        "evecs",  ".evecs" },
  { VECTRAJ,      "vectraj",".vectraj" },
  { XVG,          "xvg",    ".xvg"   },
  { UNKNOWN_DATA, 0,        0        }
};

void DataFile::WriteHelp() {
  mprintf("\t[<format keyword>]\n"
          "\t[{xlabel|ylabel|zlabel} <label>] [{xmin|ymin|zmin} <min>]\n"
          "\t[{xstep|ystep|zstep} <step>] [time <dt>] [prec <width>[.<precision>]]\n");
}

// DataFile::DetectFormat()
DataIO* DataFile::DetectFormat(FileName const& fname, DataFormatType& ftype) {
  CpptrajFile file;
  if (file.SetupRead(fname, 0) == 0) {
    for (int i = 0; i < (int)UNKNOWN_DATA; i++) {
      ftype = (DataFormatType)i;
      DataIO* IO = (DataIO*)FileTypes::AllocIO(DF_AllocArray, ftype, true );
      if (IO != 0 && IO->ID_DataFormat( file ))
        return IO;
      delete IO; 
    }
  }
  ftype = UNKNOWN_DATA;
  return 0;
}

// -----------------------------------------------------------------------------
// DataFile::SetDebug()
void DataFile::SetDebug(int debugIn) {
  debug_ = debugIn;
  if (debug_ > 0) mprintf("\tDataFile debug level set to %i\n", debug_);
}

inline int Error(const char* msg) { mprinterr(msg); return 1; }

// DataFile::ReadDataIn()
// TODO: Should this read to internal DataSetList?
int DataFile::ReadDataIn(FileName const& fnameIn, ArgList const& argListIn, 
                         DataSetList& datasetlist)
{
  if (fnameIn.empty()) return Error("Error: No input data file name given.\n"); 
  ArgList argIn = argListIn;
  if (dataio_ != 0) delete dataio_;
  dataio_ = 0;
  if (!File::Exists(fnameIn)) {
    mprinterr("Error: File '%s' does not exist.\n", fnameIn.full());
    return 1;
  }
  filename_ = fnameIn;
  // 'as' keyword specifies a format
  std::string as_arg = argIn.GetStringKey("as");
  if (!as_arg.empty()) {
    dfType_ = (DataFormatType)FileTypes::GetFormatFromString( DF_KeyArray, as_arg, UNKNOWN_DATA );
    if (dfType_ == UNKNOWN_DATA) {
      mprinterr("Error: DataFile format '%s' not recognized.\n", as_arg.c_str());
      return 1;
    }
    dataio_ = (DataIO*)FileTypes::AllocIO( DF_AllocArray, dfType_, false );
  } else
    dataio_ = DetectFormat( filename_, dfType_ );
  // Default to detection by extension.
  if (dataio_ == 0) {
    dfType_ = (DataFormatType)FileTypes::GetTypeFromExtension(DF_KeyArray, filename_.Ext(), 
                                                              DATAFILE);
    dataio_ = (DataIO*)FileTypes::AllocIO( DF_AllocArray, dfType_, false );
    if (dataio_ == 0) return Error("Error: DataIO allocation failed.\n"); 
  }
  dataio_->SetDebug( debug_ );
  // Check if user specifed DataSet name; otherwise use filename base.
  std::string dsname = argIn.GetStringKey("name");
  if (dsname.empty()) dsname = filename_.Base();
  mprintf("\tReading '%s' as %s with name '%s'\n", filename_.full(), 
          FileTypes::FormatDescription(DF_AllocArray,dfType_), dsname.c_str());
  // Read data
# ifdef TIMER
  Timer dftimer;
  dftimer.Start();
# endif
  int err = dataio_->processReadArgs(argIn);
  if (err == 0) {
    err += dataio_->ReadData( filename_, datasetlist, dsname );
    // Treat any remaining arguments as file names.
    std::string nextFile = argIn.GetStringNext();
    while (!nextFile.empty()) {
      if (filename_.SetFileName( nextFile )) return 1;
      err += dataio_->ReadData( filename_, datasetlist, dsname );
      nextFile = argIn.GetStringNext();
    }
  }
  if (err)
    mprinterr("Error: reading datafile %s\n", filename_.Full().c_str());
# ifdef TIMER
  dftimer.Stop();
  mprintf("TIME: DataFile read took %.4f seconds.\n", dftimer.Total());
# endif
  return err;
}

// DataFile::ReadDataOfType()
int DataFile::ReadDataOfType(FileName const& fnameIn, DataFormatType typeIn,
                             DataSetList& datasetlist)
{
  if (fnameIn.empty()) return Error("Error: No input data file name given.\n");
  if (dataio_ != 0) delete dataio_;
  dataio_ = 0;
  if (!File::Exists( fnameIn )) {
    mprinterr("Error: File '%s' does not exist.\n", fnameIn.full());
    return 1;
  }
  filename_ = fnameIn;
  dataio_ = (DataIO*)FileTypes::AllocIO( DF_AllocArray, typeIn, false );
  if (dataio_ == 0) return 1;
  dataio_->SetDebug( debug_ );
  return dataio_->ReadData( filename_, datasetlist, filename_.Full() );
}

// -----------------------------------------------------------------------------
// DataFile::SetupDatafile()
int DataFile::SetupDatafile(FileName const& fnameIn, ArgList& argIn, int debugIn) {
  return SetupDatafile(fnameIn, argIn, UNKNOWN_DATA, debugIn);
}

int DataFile::SetupDatafile(FileName const& fnameIn, ArgList& argIn,
                            DataFormatType typeIn, int debugIn)
{
  SetDebug( debugIn );
  if (fnameIn.empty()) return Error("Error: No data file name specified.\n");
  filename_ = fnameIn;
  dfType_ = typeIn;
  // If unknown, first look for keyword, then guess from extension.
  if (dfType_ == UNKNOWN_DATA)
    dfType_ = (DataFormatType)FileTypes::GetFormatFromArg(DF_KeyArray, argIn, UNKNOWN_DATA );
  if (dfType_ == UNKNOWN_DATA)
    dfType_ = (DataFormatType)FileTypes::GetTypeFromExtension(DF_KeyArray, filename_.Ext(),
                                                              DATAFILE);
  // Set up DataIO based on format.
  dataio_ = (DataIO*)FileTypes::AllocIO( DF_AllocArray, dfType_, false );
  if (dataio_ == 0) return Error("Error: Data file allocation failed.\n");
  if (!argIn.empty())
    ProcessArgs( argIn );
  return 0;
}

int DataFile::SetupStdout(ArgList& argIn, int debugIn) {
  SetDebug( debugIn );
  filename_.clear();
  dataio_ = (DataIO*)FileTypes::AllocIO( DF_AllocArray, DATAFILE, false );
  if (dataio_ == 0) return Error("Error: Data file allocation failed.\n");
  if (!argIn.empty())
    ProcessArgs( argIn );
  return 0;
}

// DataFile::AddDataSet()
int DataFile::AddDataSet(DataSet* dataIn) {
  if (dataIn == 0) return 1;
  if (dataio_ == 0) {
    mprinterr("Internal Error: Attempting to add set to DataFile that is not set up.\n");
    return 1;
  }
  if (SetList_.empty()) {
    dimension_ = dataIn->Ndim();
    // If current format not valid for first set, attempt to find valid format
    if (!dataio_->CheckValidFor(*dataIn)) {
      delete dataio_;
      dataio_ = 0;
      for (int dft = 0; dft != (int)UNKNOWN_DATA; dft++) {
        dfType_ = (DataFormatType)dft;
        dataio_ = (DataIO*)FileTypes::AllocIO( DF_AllocArray, dfType_, false );
        if (dataio_ == 0) return Error("Error: Data file allocation failed.\n");
        if (dataio_->CheckValidFor(*dataIn)) break;
        delete dataio_;
        dataio_ = 0;
      }
      if (dataio_ == 0) return Error("Error: Data file allocation failed.\n");
      mprintf("\tChanged DataFile '%s' type to %s for set %s\n", filename_.base(),
              FileTypes::FormatDescription(DF_AllocArray, dfType_),
              dataIn->legend());
    }
  } else {
    if ((int)dataIn->Ndim() != dimension_) {
      mprinterr("Error: DataSets in DataFile %s have dimension %i\n" 
                "Error: Attempting to add set %s of dimension %u\n", 
                filename_.base(), dimension_,
                dataIn->legend(), dataIn->Ndim());
      return Error("Error: Adding DataSets with different dimensions to same file"
                   " is currently unsupported.\n");
    }
    if (!dataio_->CheckValidFor(*dataIn)) {
      mprinterr("Error: DataSet '%s' is not valid for DataFile '%s' format.\n",
                 dataIn->legend(), filename_.base());
      return 1;
    }
  }
  // Set default width.precision
  if (setDataSetPrecision_)
    dataIn->SetupFormat().SetFormatWidthPrecision( default_width_, default_precision_ );
  // Set default label/min/step
  for (unsigned int nd = 0; nd != std::min(defaultDim_.size(), dataIn->Ndim()); nd++) {
    Dimension dim = dataIn->Dim(nd); 
    if (!defaultDim_[nd].label_.empty()) dim.SetLabel( defaultDim_[nd].label_ );
    if (defaultDim_[nd].step_ != 0.0)    dim.ChangeStep( defaultDim_[nd].step_ );
    if (minIsSet_[nd])                   dim.ChangeMin( defaultDim_[nd].min_ );
    dataIn->SetDim(nd, dim);
  }
  // Add copy of set to this DataFile
  SetList_.AddCopyOfSet( dataIn );
  // Reset dflWrite status
  dflWrite_ = true;
  return 0;
}

// DataFile::RemoveDataSet()
int DataFile::RemoveDataSet(DataSet* dataIn) {
  if (dataIn == 0) return 1;
  SetList_.RemoveSet( dataIn );
  return 0;
}

// DataFile::ProcessArgs() // FIXME make WriteArgs
int DataFile::ProcessArgs(ArgList &argIn) {
  if (dataio_==0) return 1;
  // Dimension labels 
  defaultDim_[0].label_ = argIn.GetStringKey("xlabel");
  defaultDim_[1].label_ = argIn.GetStringKey("ylabel");
  defaultDim_[2].label_ = argIn.GetStringKey("zlabel");
  // Dimension mins
  if (argIn.Contains("xmin")) {
    defaultDim_[0].min_ = argIn.getKeyDouble("xmin",1.0);
    minIsSet_[0] = true;
  }
  if (argIn.Contains("ymin")) {
    defaultDim_[1].min_ = argIn.getKeyDouble("ymin",1.0);
    minIsSet_[1] = true;
  }
  if (argIn.Contains("zmin")) {
    defaultDim_[2].min_ = argIn.getKeyDouble("zmin",1.0);
    minIsSet_[2] = true;
  }
  // Dimension steps
  defaultDim_[0].step_ = argIn.getKeyDouble("xstep", 0.0);
  defaultDim_[1].step_ = argIn.getKeyDouble("ystep", 0.0);
  defaultDim_[2].step_ = argIn.getKeyDouble("zstep", 0.0);
  // ptraj 'time' keyword
  if (argIn.Contains("time")) {
    defaultDim_[0].step_ = argIn.getKeyDouble("time", 1.0);
    defaultDim_[0].min_ = defaultDim_[0].step_;;
    minIsSet_[0] = true;
  }
  // Default DataSet width/precision
  std::string prec_str = argIn.GetStringKey("prec");
  if (!prec_str.empty()) {
    ArgList prec_arg(prec_str, ".");
    default_width_ = prec_arg.getNextInteger(-1);
    if (default_width_ < 0) {
      mprinterr("Error: Invalid width in prec arg '%s'\n", prec_str.c_str());
      return 1;
    }
    default_precision_ = prec_arg.getNextInteger(0);
    setDataSetPrecision_ = true;
  } 
  if (dataio_->processWriteArgs(argIn)==1) return 1;
  //if (debug_ > 0) argIn.CheckForMoreArgs();
  return 0;
}

// DataFile::ProcessArgs()
int DataFile::ProcessArgs(std::string const& argsIn) {
  if (argsIn.empty()) return 1;
  ArgList args(argsIn);
  return ProcessArgs(args);
}

// DataFile::WriteDataOut()
void DataFile::WriteDataOut() {
# ifdef MPI
  if (Parallel::Trajin() && !Parallel::World().Master()) {
    rprintf("DEBUG: Parallel 'trajin' mode: skipping data file write on this rank.\n");
    return;
  }
# endif
  rprintf("DEBUG: Writing file '%s'\n", DataFilename().full());
  //mprintf("DEBUG:\tFile %s has %i sets, dimension=%i, maxFrames=%i\n", dataio_->FullFileStr(),
  //        SetList_.size(), dimenison_, maxFrames);
  // Loop over all sets, decide which ones should be written.
  // All DataSets should have same dimension (enforced by AddDataSet()).
  DataSetList setsToWrite;
  for (unsigned int idx = 0; idx < SetList_.size(); ++idx) {
    DataSet& ds = static_cast<DataSet&>( *SetList_[idx] );
    // Check if set has no data.
    if ( ds.Empty() ) {
      mprintf("Warning: Set '%s' contains no data.\n", ds.legend());
      continue;
    }
    // Setup formats with a leading space initially. Maintains backwards compat. 
    ds.SetupFormat().SetFormatAlign(TextFormat::LEADING_SPACE);
    // Ensure current DataIO is valid for this set. May not be needed right now
    // but useful in case file format can be changed later on.
    if (!dataio_->CheckValidFor( ds )) {
      mprinterr("Error: DataSet '%s' is not valid for DataFile '%s' format.\n",
                 ds.legend(), filename_.base());
      continue;
    }
    setsToWrite.AddCopyOfSet( SetList_[idx] );
  }
  if (setsToWrite.empty()) {
    mprintf("Warning: File '%s' has no sets containing data.\n", filename_.base());
    return;
  }
#ifdef TIMER
  Timer dftimer;
  dftimer.Start();
#endif
  int err = dataio_->WriteData(filename_, setsToWrite); 
#ifdef TIMER
  dftimer.Stop();
  mprintf("TIME: DataFile %s Write took %.4f seconds.\n", filename_.base(),
          dftimer.Total());
#endif
  if (err > 0) 
    mprinterr("Error writing %iD Data to %s\n", dimension_, filename_.base());
}

// DataFile::SetDataFilePrecision()
/** Set precision for all DataSets in file to width.precision. */
void DataFile::SetDataFilePrecision(int widthIn, int precisionIn) {
  setDataSetPrecision_ = true;
  default_width_ = widthIn;
  default_precision_ = precisionIn;
  SetList_.SetPrecisionOfDataSets("*", widthIn, precisionIn);
}

// DataFile::DataSetNames()
/** Store DataSet names in one line. If the number of datasets is greater 
  * than 10 just print the first and last 4 data sets.
  */
std::string DataFile::DataSetNames() const {
  std::string setNames;
  DataSetList::const_iterator set = SetList_.begin();
  if (SetList_.size() > 10) {
    int setnum = 0;
    while (setnum < 4) {
      setNames.append(" " + (*set)->Meta().Legend());
      ++setnum;
      ++set;
    }
    setNames.append(" ...");
    set = SetList_.end() - 4;
    setnum = 0;
    while (setnum < 4) {
      setNames.append(" " + (*set)->Meta().Legend());
      ++setnum;
      ++set;
    }
  } else {
    for (; set != SetList_.end(); set++)
      setNames.append(" " + (*set)->Meta().Legend());
  }
  return setNames;
}
