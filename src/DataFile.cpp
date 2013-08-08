#ifdef TIMER
#include "Timer.h" 
#endif
#include "DataFile.h"
#include "CpptrajStdio.h"
// All DataIO classes go here
#include "DataIO_Std.h"
#include "DataIO_Grace.h"
#include "DataIO_Gnuplot.h"
#include "DataIO_Xplor.h"
#include "DataIO_OpenDx.h"
#include "DataIO_RemLog.h"
#include "DataIO_Mdout.h"

// TODO: Support these args:
//       - xlabel, xmin, xstep, time (all dimensions).
// CONSTRUCTOR
DataFile::DataFile() :
  debug_(0),
  dimension_(-1),
  dfType_(DATAFILE),
  dflWrite_(true),
  isInverted_(false),
  setDataSetPrecision_(false), //TODO: Just use default_width_ > -1?
  default_width_(-1),
  default_precision_(-1),
  dataio_(0),
  defaultDim_(3) // default to X/Y/Z dims
{
  for (std::vector<Dimension>::iterator dim = defaultDim_.begin(); 
                                        dim!=defaultDim_.end(); ++dim)
  {
    (*dim).SetMin(1.0);
    (*dim).SetStep(1.0);
  }
}

// DESTRUCTOR
DataFile::~DataFile() {
  if (dataio_ != 0) delete dataio_;
}

// ----- STATIC VARS / ROUTINES ------------------------------------------------
const DataFile::DataFileToken DataFile::DataFileArray[] = {
  { DATAFILE,     "dat",    "Standard Data File", ".dat",   DataIO_Std::Alloc     },
  { XMGRACE,      "grace",  "Grace File",         ".agr",   DataIO_Grace::Alloc   },
  { GNUPLOT,      "gnu",    "Gnuplot File",       ".gnu",   DataIO_Gnuplot::Alloc },
  { XPLOR,        "xplor",  "Xplor File",         ".xplor", DataIO_Xplor::Alloc   },
  { XPLOR,        "xplor",  "Xplor File",         ".grid",  DataIO_Xplor::Alloc   },
  { OPENDX,       "opendx", "OpenDx File",        ".dx",    DataIO_OpenDx::Alloc  },
  { REMLOG,       "remlog", "Amber REM log",      ".log",   DataIO_RemLog::Alloc  },
  { MDOUT,        "mdout",  "Amber MDOUT file",   ".mdout", DataIO_Mdout::Alloc   },
  { UNKNOWN_DATA, 0,        "Unknown",            0,        0                     }
};

// DataFile::GetFormatFromArg()
/** Given an ArgList, search for one of the file format keywords.
  */
DataFile::DataFormatType DataFile::GetFormatFromArg(ArgList& argIn) 
{
  for (TokenPtr token = DataFileArray; token->Type != UNKNOWN_DATA; ++token)
    if (argIn.hasKey( token->Key )) return token->Type;
  return UNKNOWN_DATA;
}

// DataFile::GetFormatFromString()
DataFile::DataFormatType DataFile::GetFormatFromString(std::string const& fmt)
{
  for (TokenPtr token = DataFileArray; token->Type != UNKNOWN_DATA; ++token)
    if ( fmt.compare( token->Key )==0 ) return token->Type;
  return DATAFILE;
}

// DataFile::GetExtensionForType()
std::string DataFile::GetExtensionForType(DataFormatType typeIn) {
  for (TokenPtr token = DataFileArray; token->Type != UNKNOWN_DATA; ++token)
    if ( token->Type == typeIn )
      return std::string( token->Extension );
  return std::string();
}

// DataFile::GetTypeFromExtension()
DataFile::DataFormatType DataFile::GetTypeFromExtension( std::string const& extIn)
{
  for (TokenPtr token = DataFileArray; token->Type != UNKNOWN_DATA; ++token)
    if ( extIn.compare( token->Extension ) == 0 ) return token->Type;
  // Default to DATAFILE
  return DATAFILE;
}

// DataFile::FormatString()
const char* DataFile::FormatString() const {
  TokenPtr token;
  for (token = DataFileArray; token->Type != UNKNOWN_DATA; ++token)
    if ( token->Type == dfType_ ) return token->Description;
  return token->Description; // Should be at UNKNOWN
}
// -----------------------------------------------------------------------------

// DataFile::SetDebug()
void DataFile::SetDebug(int debugIn) {
  debug_ = debugIn;
  if (debug_ > 0) mprintf("\tDataFile debug level set to %i\n", debug_);
}

// ----- DATA FILE ALLOCATION / DETECTION ROUTINES -----------------------------
// DataFile::AllocDataIO()
DataIO* DataFile::AllocDataIO(DataFormatType tformat) {
  for (TokenPtr token = DataFileArray; token->Type != UNKNOWN_DATA; ++token) {
    if (token->Type == tformat) {
      if (token->Alloc == 0) {
        mprinterr("Error: CPPTRAJ was compiled without support for %s files.\n",
                  token->Description);
        return 0;
      } else
        return (DataIO*)token->Alloc();
    }
  }
  return 0;
}

// DataFile::DetectFormat()
DataIO* DataFile::DetectFormat(std::string const& fname, DataFormatType& ftype) {
  CpptrajFile file;
  if (file.SetupRead(fname, 0)) return 0;
  for (TokenPtr token = DataFileArray; token->Type != UNKNOWN_DATA; ++token) {
    if (token->Alloc != 0) {
      DataIO* io = (DataIO*)token->Alloc();
      if ( io->ID_DataFormat( file ) ) {
        ftype = token->Type;
        return io;
      }
      delete io;
    }
  }
  ftype = UNKNOWN_DATA;
  return 0;
}

// DataFile::DataFormat()
/*DataFile::DataFormatType DataFile::DataFormat(std::string const& fname) {
  CpptrajFile file;
  if (file.SetupRead(fname, 0)) return UNKNOWN_DATA;
  for (TokenPtr token = DataFileArray; token->Type != UNKNOWN_DATA; ++token) {
    if (token->Alloc != 0) {
      DataIO* io = (DataIO*)token->Alloc();
      if ( io->ID_DataFormat( file ) ) {
        delete io;
        return token->Type;
      }
      delete io;
    }
  }
  return UNKNOWN_DATA;
}*/
// -----------------------------------------------------------------------------

// DataFile::ReadDataIn()
// TODO: Should this read to internal DataSetList?
int DataFile::ReadDataIn(std::string const& fnameIn, ArgList const& argListIn, 
                         DataSetList& datasetlist)
{
  if (fnameIn.empty()) {
    mprinterr("Error: No input data file name given.\n");
    return 1;
  }
  ArgList argIn = argListIn;
  if (dataio_ != 0) delete dataio_;
  dataio_ = 0;
  filename_.SetFileNameWithExpansion( fnameIn );
  // 'as' keyword specifies a format
  std::string as_arg = argIn.GetStringKey("as");
  if (!as_arg.empty()) {
    dfType_ = GetFormatFromString( as_arg );
    if (dfType_ == UNKNOWN_DATA) {
      mprinterr("Error: DataFile format %s not recognized.\n", as_arg.c_str());
      return 1;
    }
    dataio_ = AllocDataIO( dfType_ );
  } else
    dataio_ = DetectFormat( filename_.Full(), dfType_ );
  // Default to detection by extension.
  if (dataio_ == 0) {
    dfType_ = GetTypeFromExtension(filename_.Ext());
    dataio_ = AllocDataIO( dfType_ );
  }
  mprintf("\tReading %s as %s\n", filename_.full(), FormatString());
  // Check if user specifed DataSet name; otherwise use filename base.
  std::string dsname = argIn.GetStringKey("name");
  if (dsname.empty()) dsname = filename_.Base();
  // Read data
# ifdef TIMER
  Timer dftimer;
  dftimer.Start();
# endif
  int err = dataio_->ReadData( filename_.Full(), argIn, datasetlist, dsname );
  if (err)
    mprinterr("Error reading datafile %s\n", filename_.Full().c_str());
# ifdef TIMER
  dftimer.Stop();
  mprintf("TIME: DataFile read took %.4f seconds.\n", dftimer.Total());
# endif
  return err;
}

inline int Error(const char* msg) { mprinterr(msg); return 1; }

// DataFile::SetupDatafile()
int DataFile::SetupDatafile(std::string const& fnameIn, ArgList& argIn, int debugIn) {
  SetDebug( debugIn );
  if (fnameIn.empty()) return Error("Error: No data file name specified.\n");
  filename_.SetFileName( fnameIn );
  dfType_ = GetFormatFromArg( argIn );
  if (dfType_ == UNKNOWN_DATA)
    dfType_ = GetTypeFromExtension(filename_.Ext());
  // Set up DataIO based on format.
  dataio_ = AllocDataIO( dfType_ );
  if (dataio_ == 0) return Error("Error: Data file allocation failed.\n");
  if (!argIn.empty())
    ProcessArgs( argIn );
  return 0;
}

// DataFile::AddSet()
int DataFile::AddSet(DataSet* dataIn) {
  if (dataIn == 0) return 1;
  if (SetList_.empty())
    dimension_ = dataIn->Ndim();
  else if ((int)dataIn->Ndim() != dimension_) {
    mprinterr("Error: DataSets in DataFile %s have dimension %i\n" 
              "Error: Attempting to add set %s of dimension %u\n", 
              filename_.base(), dimension_,
              dataIn->Legend().c_str(), dataIn->Ndim());
    return Error("Error: Adding DataSets with different dimensions to same file"
                 " is currently unsupported.\n");
  }
  // Set default width.precision
  if (setDataSetPrecision_)
    dataIn->SetPrecision( default_width_, default_precision_ );
  SetList_.AddCopyOfSet( dataIn );
  // Reset dflWrite status
  dflWrite_ = true;
  return 0;
}

// DataFile::ProcessArgs()
int DataFile::ProcessArgs(ArgList &argIn) {
  if (dataio_==0) return 1;
  if (argIn.hasKey("invert")) {
    isInverted_ = true;
    // Currently GNUPLOT files cannot be inverted.
    if (dfType_ == GNUPLOT) {
      mprintf("Warning: (%s) Gnuplot files cannot be inverted.\n",filename_.base());
      isInverted_ = false;;
    }
  }
  // Axis args.
  defaultDim_[0].SetLabel( argIn.GetStringKey("xlabel") );
  defaultDim_[1].SetLabel( argIn.GetStringKey("ylabel") );
  // Axis min/step
  defaultDim_[0].SetMin( argIn.getKeyDouble("xmin",1.0) );
  defaultDim_[1].SetMin( argIn.getKeyDouble("ymin",1.0) );
  defaultDim_[0].SetStep( argIn.getKeyDouble("xstep", 1.0) );
  defaultDim_[1].SetStep( argIn.getKeyDouble("ystep", 1.0) );
  // ptraj 'time' keyword
  if (argIn.Contains("time")) {
    defaultDim_[0].SetStep( argIn.getKeyDouble("time", 1.0) );
    defaultDim_[0].SetMin( defaultDim_[0].Step() );
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
  if (debug_ > 0) argIn.CheckForMoreArgs();
  return 0;
}

// DataFile::ProcessArgs()
int DataFile::ProcessArgs(std::string const& argsIn) {
  if (argsIn.empty()) return 1;
  ArgList args(argsIn);
  return ProcessArgs(args);
}

// DataFile::WriteData()
void DataFile::WriteData() {
  // Remove data sets that do not contain data.
  // All DataSets should have same dimension (enforced by AddSet()).
  DataSetList::const_iterator Dset = SetList_.end();
  while (Dset != SetList_.begin()) {
    --Dset;
    // Check if set has no data.
    if ( (*Dset)->Empty() ) {
      // If set has no data, remove it
      mprintf("Warning: Set %s contains no data. Skipping.\n",(*Dset)->Legend().c_str());
      SetList_.erase( Dset );
      Dset = SetList_.end();
    } else {
      // If set has data, set its format to right-aligned initially.
      if ( (*Dset)->SetDataSetFormat(false) ) {
        mprinterr("Error: could not set format string for set %s. Skipping.\n", 
                  (*Dset)->Legend().c_str());
        SetList_.erase( Dset );
        Dset = SetList_.end();
      } 
    }
  }
  // If all data sets are empty no need to write
  if (SetList_.empty()) {
    mprintf("Warning: file %s has no sets containing data.\n", filename_.base());
    return;
  }
  //mprintf("DEBUG:\tFile %s has %i sets, dimension=%i, maxFrames=%i\n", dataio_->FullFileStr(),
  //        SetList_.size(), dimenison_, maxFrames);
  // Loop over all sets, set default min, step, label if not already set.
  // Set default min and step for all dimensions if not already set.
  for (unsigned int idx = 0; idx < SetList_.size(); ++idx) {
    DataSet& ds = static_cast<DataSet&>( *SetList_[idx] );
    for (unsigned int nd = 0; nd < ds.Ndim(); ++nd) {
      Dimension& dim = ds.Dim(nd);
      double min, step;
      if (nd < 3) {
        min = defaultDim_[nd].Min();
        step = defaultDim_[nd].Step();
      } else {
        min = 1.0;
        step = 1.0;
      }
      if (!dim.MinIsSet()) dim.SetMin( min );
      if (dim.Step() < 0 ) dim.SetStep( step );
    }
    // Set x label if not already set
    if (ds.Ndim()==1 && ds.Dim(0).Label().empty()) ds.Dim(0).SetLabel("Frame");
  }
#ifdef TIMER
  Timer dftimer;
  dftimer.Start();
#endif
  int err = 0;
  if ( dimension_ == 1 ) {       // One-dimensional
    if (!isInverted_)
      err = dataio_->WriteData(filename_.Full(), SetList_);
    else
      err = dataio_->WriteDataInverted(filename_.Full(), SetList_);
  } else if ( dimension_ == 2) { // Two-dimensional
    for ( DataSetList::const_iterator set = SetList_.begin();
                                      set != SetList_.end(); ++set)
      err += dataio_->WriteData2D(filename_.Full(), *(*set) );
  } else if ( dimension_ == 3) { // Three-dimensional
    for ( DataSetList::const_iterator set = SetList_.begin();
                                      set != SetList_.end(); ++set)
      err += dataio_->WriteData3D(filename_.Full(), *(*set) );
  } else {
    mprinterr("Error: %iD writes not yet supported.\n", dimension_);
    err = 1;
  }
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
/** Print Dataset names to one line. If the number of datasets is greater 
  * than 10 just print the first and last 4 data sets.
  */
void DataFile::DataSetNames() const {
  DataSetList::const_iterator set = SetList_.begin();
  if (SetList_.size() > 10) {
    int setnum = 0;
    while (setnum < 4) {
      mprintf(" %s",(*set)->Legend().c_str());
      ++setnum;
      ++set;
    }
    mprintf(" ...");
    set = SetList_.end() - 4;
    setnum = 0;
    while (setnum < 4) {
      mprintf(" %s",(*set)->Legend().c_str());
      ++setnum;
      ++set;
    }
  } else {
    for (; set != SetList_.end(); set++)
      mprintf(" %s",(*set)->Legend().c_str());
  }
}
