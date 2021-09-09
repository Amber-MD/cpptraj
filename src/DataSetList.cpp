#include "DataSetList.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // DigitWidth
#include "Constants.h" // For AddOrAppendSets
#include "ArgList.h"
// Data types go here
#include "DataSet_double.h"
#include "DataSet_float.h"
#include "DataSet_integer_mem.h"
#ifdef BINTRAJ
#include "DataSet_integer_disk.h"
#endif
#include "DataSet_string.h"
#include "DataSet_MatrixDbl.h"
#include "DataSet_MatrixFlt.h"
#include "DataSet_Coords_CRD.h"
#include "DataSet_Vector.h"
#include "DataSet_Modes.h"
#include "DataSet_GridFlt.h"
#include "DataSet_RemLog.h"
#include "DataSet_Mesh.h"
#include "DataSet_Coords_TRJ.h"
#include "DataSet_Coords_REF.h"
#include "DataSet_Mat3x3.h"
#include "DataSet_Topology.h"
#include "DataSet_GridDbl.h"
#include "DataSet_pH.h"
#include "DataSet_PHREMD_Explicit.h"
#include "DataSet_PHREMD_Implicit.h"
#include "DataSet_Parameters.h"
#include "DataSet_PairwiseCache_MEM.h"
#include "DataSet_PairwiseCache_NC.h"
#include "DataSet_Tensor.h"
#include "DataSet_StringVar.h"
#include "DataSet_Vector_Scalar.h"
#include "DataSet_unsignedInt.h"

bool DataSetList::useDiskCache_ = false;

/** Master data set allocation routine. */
DataSet* DataSetList::NewSet(DataSet::DataType typeIn) {
  DataSet* ds = 0;
  bool cannotUseDiskCache = true; // Most cannot.
  switch (typeIn) {
    case DataSet::UNKNOWN_DATA :
      mprinterr("Internal Error: DataSetList::NewSet() called with no type.\n"); break;
    case DataSet::DOUBLE  : ds = DataSet_double::Alloc(); break;
    case DataSet::FLOAT   : ds = DataSet_float::Alloc(); break;
    case DataSet::INTEGER :
#     ifdef BINTRAJ
      if (useDiskCache_) {
        ds = DataSet_integer_disk::Alloc();
        cannotUseDiskCache = false;
      } else
        ds = DataSet_integer_mem::Alloc();
#     else
      if (useDiskCache_)
        mprintf("Warning: Integer data set disk cache requires NetCDF. Using memory.\n");
      ds = DataSet_integer_mem::Alloc();
#     endif
      break;
    case DataSet::UNSIGNED_INTEGER : ds = DataSet_unsignedInt::Alloc(); break;
    case DataSet::STRING  : ds = DataSet_string::Alloc(); break;
    case DataSet::MATRIX_DBL : ds = DataSet_MatrixDbl::Alloc(); break;
    case DataSet::MATRIX_FLT : ds = DataSet_MatrixFlt::Alloc(); break;
    case DataSet::COORDS     : ds = DataSet_Coords_CRD::Alloc(); break;
    case DataSet::VECTOR     : ds = DataSet_Vector::Alloc() ; break;
    case DataSet::MODES      : ds = DataSet_Modes::Alloc(); break;
    case DataSet::GRID_FLT   : ds = DataSet_GridFlt::Alloc(); break;
    case DataSet::GRID_DBL   : ds = DataSet_GridDbl::Alloc(); break;
    case DataSet::REMLOG     : ds = DataSet_RemLog::Alloc(); break;
    case DataSet::XYMESH     : ds = DataSet_Mesh::Alloc(); break;
    case DataSet::TRAJ       : ds = DataSet_Coords_TRJ::Alloc(); break;
    case DataSet::REF_FRAME  : ds = DataSet_Coords_REF::Alloc(); break;
    case DataSet::MAT3X3     : ds = DataSet_Mat3x3::Alloc(); break;
    case DataSet::TOPOLOGY   : ds = DataSet_Topology::Alloc(); break;
    case DataSet::PH            : ds = DataSet_pH::Alloc(); break;
    case DataSet::PH_EXPL       : ds = DataSet_PHREMD_Explicit::Alloc(); break;
    case DataSet::PH_IMPL       : ds = DataSet_PHREMD_Implicit::Alloc(); break;
    case DataSet::PARAMETERS    : ds = DataSet_Parameters::Alloc(); break;
    // TODO useDiskCache
    case DataSet::PMATRIX_MEM   : ds = DataSet_PairwiseCache_MEM::Alloc(); break;
    case DataSet::PMATRIX_NC    :
#     ifdef BINTRAJ
      ds = DataSet_PairwiseCache_NC::Alloc();
#     else
      mprinterr("Error: NetCDF pairwise disk cache requires NetCDF.\n");
      ds = 0;
#     endif
      break;
    case DataSet::TENSOR        : ds = DataSet_Tensor::Alloc(); break;
    case DataSet::STRINGVAR     : ds = DataSet_StringVar::Alloc(); break;
    case DataSet::VECTOR_SCALAR : ds = DataSet_Vector_Scalar::Alloc(); break;
    // Sanity check
    default:
      mprinterr("Internal Error: No allocator for DataSet type '%s'\n",
                DataSet::description(typeIn));
  }
  if (ds == 0)
    mprinterr("Error: Could not allocate DataSet type '%s'\n", DataSet::description(typeIn));
  else if (useDiskCache_ && cannotUseDiskCache)
    mprintf("Warning: Use disk cache specified, but DataSet type '%s' cannot use disk caching.\n",
            DataSet::description(typeIn));
  return ds;
}

// CONSTRUCTOR
DataSetList::DataSetList() :
  activeRef_(0),
  maxFrames_(-1),
  debug_(0),
  ensembleNum_(-1),
  hasCopies_(false),
  dataSetsPending_(false)
# ifdef MPI
  , newSetsNeedSync_(false)
# endif
{}

// DESTRUCTOR
DataSetList::~DataSetList() { ClearAll(); }

// DataSetList::ClearAll()
void DataSetList::ClearAll() {
  if (!hasCopies_)
    for (DataListType::iterator ds = DataList_.begin(); ds != DataList_.end(); ++ds) 
      delete *ds;
  DataList_.clear();
  RefList_.clear();
  TopList_.clear();
  hasCopies_ = false;
  dataSetsPending_ = false;
# ifdef MPI
  newSetsNeedSync_ = false;
# endif
  activeRef_ = 0;
} 

// DataSetList::Clear()
void DataSetList::Clear() {
  DataListType setsToErase;
  DataListType setsToKeep;
  setsToKeep.reserve( RefList_.size() + TopList_.size() );
  for (DataListType::const_iterator ds = DataList_.begin(); ds != DataList_.end(); ++ds)
    if ( (*ds)->Type() == DataSet::REF_FRAME || (*ds)->Type() == DataSet::TOPOLOGY )
      setsToKeep.push_back( *ds );
    else
      setsToErase.push_back( *ds );
  if (!hasCopies_)
    for (DataListType::iterator ds = setsToErase.begin(); ds != setsToErase.end(); ++ds)
      delete *ds;
  DataList_ = setsToKeep;
}

void DataSetList::Sort() {
  std::sort( DataList_.begin(), DataList_.end(), DataSet::DS_PtrCmp() );
}

// DataSetList::Push_Back()
void DataSetList::Push_Back(DataSet* ds) {
  DataList_.push_back( ds );
  if (!hasCopies_) {
    // If this is a REF data set it also goes in RefList_.
    if (ds->Type() == DataSet::REF_FRAME) { 
      RefList_.push_back( ds );
      if (activeRef_ == 0) SetActiveReference( ds );
    } else if (ds->Type() == DataSet::TOPOLOGY) {
      ((DataSet_Topology*)ds)->SetPindex( TopList_.size() );
      TopList_.push_back( ds );
    }
  }
}

// DataSetList::operator+=()
DataSetList& DataSetList::operator+=(DataSetList const& rhs) {
  // It is OK if rhs does not have copies, but this should have copies.
  // For now just set hasCopies to true.
  hasCopies_ = true; 
  // Append rhs entries to here
  for (DataListType::const_iterator DS = rhs.begin(); DS != rhs.end(); ++DS)
    Push_Back( *DS );
  return *this;
}

// DataSetList::MakeDataSetsEnsemble()
/*
void DataSetList::MakeDataSetsEnsemble(int ensembleNumIn) {
  ensembleNum_ = ensembleNumIn;
  for (DataListType::const_iterator ds = DataList_.begin();
                                    ds != DataList_.end(); ++ds)
    if ( (*ds)->Meta().EnsembleNum() == -1 )
      (*ds)->SetEnsemble( ensembleNum_ );
}
*/

/** Call Allocate for each time series in the list. */
void DataSetList::AllocateSets(long int maxFrames) {
  maxFrames_ = maxFrames;
  if (maxFrames < 1L) return;
  DataSet::SizeArray mfArray(1, maxFrames);
  for (DataListType::iterator ds = DataList_.begin(); ds != DataList_.end(); ++ds)
  {
    if ( (*ds)->Meta().TimeSeries() == MetaData::IS_TS )
      if ( (*ds)->Allocate( mfArray ) )
        mprinterr("Error: Could not allocate time series for '%s'\n", (*ds)->legend());
  }
}

// DataSetList::SetPrecisionOfDataSets()
void DataSetList::SetPrecisionOfDataSets(std::string const& nameIn, int widthIn,
                                         int precisionIn)
{
  if (widthIn < 1)
    mprinterr("Error: Invalid data width (%i)\n", widthIn);
  else {
    DataSetList Sets = SelectSets( nameIn );
    for (DataSetList::const_iterator ds = Sets.begin(); ds != Sets.end(); ++ds)
      (*ds)->SetupFormat().SetFormatWidthPrecision(widthIn, precisionIn);
  }
}

// DataSetList::PendingWarning()
void DataSetList::PendingWarning() const {
  if (dataSetsPending_)
    mprintf("Warning: Some Actions currently in Action list need to be run in order to create\n"
            "Warning:   data sets. Try processing currently loaded trajectories with 'run' or\n"
            "Warning:   'go' to generate these data sets.\n");
}

// -----------------------------------------------------------------------------
/** Remove the specified set from all lists if found.
  * \param freeMemory If true, also free DataSet memory.
  * \return dsIn if found and removed/erased, 0 otherwise.
  */
DataSet* DataSetList::EraseSet( DataSet* dsIn, bool freeMemory ) {
  if (dsIn == 0) return 0;
  for (DataListType::iterator pos = DataList_.begin();
                              pos != DataList_.end(); ++pos)
  {
    if ( *pos == dsIn ) {
      // Also remove from reference list if applicable.
      if ( dsIn->Type() == DataSet::REF_FRAME ) {
        for (DataListType::iterator ref = RefList_.begin(); ref != RefList_.end(); ++ref)
          if ( *ref == dsIn ) {
            RefList_.erase( ref );
            break;
          }
      } else if (dsIn->Type() == DataSet::TOPOLOGY ) {
        for (DataListType::iterator top = TopList_.begin(); top != TopList_.end(); ++top)
          if ( *top == dsIn ) {
            TopList_.erase( top );
            break;
          }
        // Reset P indices so they are unique.
        for (DataListType::iterator top = TopList_.begin(); top != TopList_.end(); ++top)
          ((DataSet_Topology*)*top)->SetPindex( top - TopList_.begin() );
      }
      if (!hasCopies_ && freeMemory) delete *pos;
      DataList_.erase( pos );
      return dsIn; // NOTE: This is invalid if memory was freed
    }
  }
  return 0;
}

/** Remove DataSet from the list and destroy it. */
void DataSetList::RemoveSet( DataSet* dsIn ) {
  EraseSet( dsIn, true );
}

/** \return dsIn if found and removed from list, 0 if not found. */
DataSet* DataSetList::PopSet( DataSet* dsIn ) {
  return EraseSet( dsIn, false );
}

// -----------------------------------------------------------------------------
/** The set argument must match EXACTLY, so Data will not return Data:1 */
DataSet* DataSetList::CheckForSet(MetaData const& md) const
{
  // If incoming MetaData ensemble number is not set but DataSetList is,
  // use DataSetList member number.
  if (md.EnsembleNum() == -1 && ensembleNum_ != -1) {
    MetaData md_ensNum( md );
    md_ensNum.SetEnsembleNum( ensembleNum_ );
    for (DataListType::const_iterator ds = DataList_.begin(); ds != DataList_.end(); ++ds)
      if ( (*ds)->Matches_Exact( md_ensNum ) )
        return *ds;
  } else
    for (DataListType::const_iterator ds = DataList_.begin(); ds != DataList_.end(); ++ds)
      if ( (*ds)->Matches_Exact( md ) )
        return *ds;
  return 0;
}

// DataSetList::GetDataSet()
DataSet* DataSetList::GetDataSet( std::string const& nameIn ) const {
  DataSetList dsetOut = SelectSets( nameIn );
  if (dsetOut.empty()) {
    mprintf("Warning: Data set '%s' not found.\n", nameIn.c_str());
    PendingWarning();
    return 0;
  } else if (dsetOut.size() > 1)
    mprintf("Warning: '%s' selects multiple sets, only using first set.\n", nameIn.c_str());
  return dsetOut[0];
}

// DataSetList::GetMultipleSets()
DataSetList DataSetList::GetMultipleSets( std::string const& dsargIn ) const {
  DataSetList dsetOut = SelectSets(dsargIn);
  if ( dsetOut.empty() ) {
    mprintf("Warning: '%s' selects no data sets.\n", dsargIn.c_str());
    PendingWarning();
  }
  return dsetOut;
}

// DataSetList::GetSetsOfType()
DataSetList DataSetList::GetSetsOfType( std::string const& dsargIn, DataSet::DataType typeIn ) const
{
  DataSetList dsetOut;
  dsetOut.hasCopies_ = true;
  DataSetList selected = SelectSets(dsargIn);
  for (const_iterator ds = selected.begin(); ds != selected.end(); ++ds)
    if ( (*ds)->Type() == typeIn )
      dsetOut.Push_Back( *ds );
  return dsetOut;
}

// DataSetList::SelectSets()
DataSetList DataSetList::SelectSets( std::string const& nameIn ) const {
  return SelectSets( nameIn, DataSet::UNKNOWN_DATA );
}

/** \return a list of all DataSets matching the given argument. */
DataSetList DataSetList::SelectSets( std::string const& dsargIn,
                                     DataSet::DataType typeIn ) const
{
  DataSetList dsetOut;
  dsetOut.hasCopies_ = true;
  // Find matching sets.
  MetaData::SearchString search(dsargIn);
  for (DataListType::const_iterator ds = DataList_.begin(); ds != DataList_.end(); ++ds)
    if ((*ds)->Matches_WC( search, typeIn ))
      dsetOut.Push_Back( *ds );

  return dsetOut;
}

/** \return a list of all DataSets matching given argument and group. */
DataSetList DataSetList::SelectGroupSets( std::string const& dsargIn,
                                          DataSet::DataGroup groupIn ) const
{
  DataSetList dsetOut;
  dsetOut.hasCopies_ = true;
  MetaData::SearchString search(dsargIn);
  for (DataListType::const_iterator ds = DataList_.begin(); ds != DataList_.end(); ++ds)
    if ((*ds)->Group() == groupIn && (*ds)->Matches_WC( search, DataSet::UNKNOWN_DATA ))
      dsetOut.Push_Back( *ds );
  return dsetOut;
}

// DataSetList::FindSetOfType()
DataSet* DataSetList::FindSetOfType(std::string const& nameIn, DataSet::DataType typeIn) const
{
  DataSetList dsetOut = SelectSets( nameIn, typeIn );
  if (dsetOut.empty())
    return 0;
  else if (dsetOut.size() > 1)
    mprintf("Warning: '%s' selects multiple sets. Only using first.\n", nameIn.c_str());
  return dsetOut[0];
}

/** \return single set belonging to specified group. */
DataSet* DataSetList::FindSetOfGroup(std::string const& nameIn, DataSet::DataGroup groupIn) const
{
  DataSetList dsetOut = SelectGroupSets( nameIn, groupIn );
  if (dsetOut.empty())
    return 0;
  else if (dsetOut.size() > 1)
    mprintf("Warning: '%s' selects multiple sets. Only using first.\n", nameIn.c_str());
  return dsetOut[0];
}

/** Search for a COORDS DataSet. If no name specified, create a default
  * COORDS DataSet named _DEFAULTCRD_.
  */
DataSet* DataSetList::FindCoordsSet(std::string const& setname) {
  DataSet* outset = 0;
  if (setname.empty()) {
    // crdset not given, search for the default set
    outset = FindSetOfType("_DEFAULTCRD_", DataSet::COORDS);
    if (outset == 0) {
      // No default set; create one.
      outset = AddSet(DataSet::COORDS, "_DEFAULTCRD_", "CRD");
    }
  } else {
    DataSetList dslist = SelectGroupSets(setname, DataSet::COORDINATES);
    if (!dslist.empty()) outset = dslist[0];
  }
  return outset;
}

// -----------------------------------------------------------------------------
/** Create a name based on the given defaultName and # of DataSets,
  * i.e. defaultName_XXXXX 
  */
std::string DataSetList::GenerateDefaultName(std::string const& defaultName) const {
  // Determine # chars needed to hold text version of set number (min 5).
  size_t extsize = (size_t) DigitWidth( size() );
  if (extsize < 5) extsize = 5;
  if (defaultName.empty())
    return ( "D" + integerToString(size(), extsize) );
  else
    return ( defaultName + "_" + integerToString(size(), extsize) ); 
}

/** Add a DataSet with given MetaData; if no name given create a name based on
  * defaultName and DataSet position.
  */
DataSet* DataSetList::AddSet( DataSet::DataType inType, MetaData const& metaIn,
                              const char* defaultName )
{
  MetaData meta = metaIn;
  if (meta.Name().empty()) {
    if (defaultName != 0)
      meta.SetName( GenerateDefaultName(defaultName) );
  }
  return AddSet( inType, meta );
}

/** Allocate a DataSet of specified type, set it up and return pointer to it.
  * It will not be added to the list. 
  * \param inType type of DataSet to add.
  * \param metaIn DataSet MetaData.
  * \return pointer to successfully set-up DataSet or 0 if error.
  */ 
DataSet* DataSetList::AllocateSet(DataSet::DataType inType, MetaData const& metaIn)
const
{ // TODO Always generate default name if empty?
# ifdef TIMER
  time_total_.Start();
# endif
  // Do not add to a list with copies
  if (hasCopies_) {
    mprinterr("Internal Error: Attempting to add DataSet (%s) to DataSetList with copies.\n",
              metaIn.PrintName().c_str());
    return 0;
  }
  MetaData meta( metaIn );
  meta.SetEnsembleNum( ensembleNum_ );
# ifdef TIMER
  time_check_.Start();
# endif
  // Check if DataSet with same attributes already present.
  DataSet* DS = CheckForSet(meta);
# ifdef TIMER
  time_check_.Stop();
# endif
  if (DS != 0) {
    mprintf("Warning: DataSet '%s' already present.\n", DS->Meta().PrintName().c_str());
    // NOTE: Should return found dataset?
    return 0; 
  }
  DS = Allocate( inType );
  if (DS==0) {
    mprinterr("Internal Error: DataSet %s memory allocation failed.\n", meta.PrintName().c_str());
    return 0;
  }
# ifdef TIMER
  time_setup_.Start();
# endif
  //mprintf("DEBUG: Adding and setting up data set '%s'\n", meta.PrintName().c_str());
  //if (meta.TimeSeries() == MetaData::UNKNOWN_TS)
  //  mprintf("DEBUG: UNKNOWN_TS\n");
  //mprintf("DEBUG: ndim= %i\n", DS->Ndim());
  // If 1 dim set and time series status not set, set to true.
  if (meta.TimeSeries() == MetaData::UNKNOWN_TS && DS->Ndim() == 1) {
    meta.SetTimeSeries( MetaData::IS_TS );
    // Also set dimension default
    DS->SetDim(Dimension::X, Dimension(1.0, 1.0, "Frame") );
  }
  // Set up dataset 
  if ( DS->SetMeta( meta ) ) {
    mprinterr("Error setting up data set %s.\n", meta.PrintName().c_str());
    delete DS;
    return 0;
  }
# ifdef MPI
  if (newSetsNeedSync_) DS->SetNeedsSync( true );
# endif
# ifdef TIMER
  time_setup_.Stop();
# endif
  //fprintf(stderr,"ALLOCATED dataset %s\n",dsetName);
  return DS;
}

/** Add a DataSet of specified type, set it up and return pointer to it. 
  * \param inType type of DataSet to add.
  * \param metaIn DataSet MetaData.
  * \return pointer to successfully set-up DataSet or 0 if error.
  */ 
DataSet* DataSetList::AddSet(DataSet::DataType inType, MetaData const& metaIn)
{
  DataSet* DS = AllocateSet(inType, metaIn);
  if (DS != 0) {
#   ifdef TIMER
    time_push_.Start();
#   endif
    Push_Back(DS);
#   ifdef TIMER
    time_push_.Stop();
#   endif
    //fprintf(stderr,"ADDED dataset %s\n",dsetName);
  }
  return DS;
}

#ifdef TIMER
void DataSetList::Timing() const {
  time_check_.WriteTiming(3,"DSL Check Set", time_total_.Total());
  time_setup_.WriteTiming(3,"DSL Setup Set", time_total_.Total());
  time_push_.WriteTiming(3,"DSL Push Set", time_total_.Total());
  time_total_.WriteTiming(2,"DSL Total");
}
#endif

// DataSetList::Allocate()
DataSet* DataSetList::Allocate(DataSet::DataType inType) {
  return NewSet(inType);
}

// FIXME Should probably just make a more efficient search of DSL
/** Special version of AddSet that does NOT check if set already exists.
  * Intended for use in Action Setup/DoAction where it is assumed that
  * the Action is setting up DataSets in such a way that there will not
  * be name conflicts, i.e. the DataSet name at least is unique.
  * \param inType type of DataSet to add.
  * \param metaIn DataSet MetaData.
  * \return pointer to successfully set-up DataSet or 0 if error.
  */
DataSet* DataSetList::AddSet_NoCheck(DataSet::DataType inType, MetaData const& metaIn)
{ // TODO Pass in Nframes?
  // Assume list does NOT have copies.
  MetaData meta( metaIn );
  meta.SetEnsembleNum( ensembleNum_ );
  // Allocate DataSet
  DataSet* DS = Allocate(inType);
  if (DS==0) {
    mprinterr("Internal Error: DataSet %s memory allocation failed.\n", meta.PrintName().c_str());
    return 0;
  }
  // If 1 dim set and time series status not set, set to true, allocate for frames.
  if (meta.TimeSeries() == MetaData::UNKNOWN_TS && DS->Ndim() == 1) {
    meta.SetTimeSeries( MetaData::IS_TS );
    // Also set dimension default
    DS->SetDim(Dimension::X, Dimension(1.0, 1.0, "Frame") );
    //DS->Allocate( DataSet::SizeArray(1, Nframes) );
  }
  // Set up DataSet MetaData
  if ( DS->SetMeta( meta ) ) {
    mprinterr("Error setting up data set %s.\n", meta.PrintName().c_str());
    delete DS;
    return 0;
  }
# ifdef MPI
  if (newSetsNeedSync_) DS->SetNeedsSync( true );
# endif
  // Add to list
  Push_Back(DS);
  return DS;
}

// DataSetList::AddSet()
int DataSetList::AddSet( DataSet* dsIn ) {
  if (dsIn == 0 ) return 1;
  DataSet* ds = CheckForSet( dsIn->Meta() );
  if (ds != 0) {
    mprintf("Warning: DataSet '%s' already present.\n", ds->Meta().PrintName().c_str());
    return 1;
  }
  Push_Back( dsIn );
  return 0;
}

/** Given an array of already set up data sets (from e.g. input data files)
  * and optional X values, add the sets to this list if they do not exist
  * or append any existing sets.
  */
int DataSetList::AddOrAppendSets(std::string const& XlabelIn, Darray const& Xvals,
                                 DataListType const& Sets)
{
  if (debug_ > 0)
    mprintf("DEBUG: Calling AddOrAppendSets for %zu sets, %zu X values, Xlabel= '%s'.\n",
            Sets.size(), Xvals.size(), XlabelIn.c_str());
  if (Sets.empty()) return 0; // No error for now.
  // If no X label assume 'Frame' for backwards compat.
  std::string Xlabel;
  if (XlabelIn.empty())
    Xlabel.assign("Frame");
  else
    Xlabel = XlabelIn;
  Dimension Xdim;
  // First determine if X values increase monotonically with a regular step
  bool isMonotonic = Xdim.SetDimension(Xvals, Xlabel);
  if (debug_ > 0) {
    mprintf("DEBUG: xstep %g xmin %g\n", Xdim.Step(), Xdim.Min());
    if (isMonotonic) mprintf("DEBUG: Xdim is monotonic.\n");
  }
  for (const_iterator ds = Sets.begin(); ds != Sets.end(); ++ds) {
    if (*ds == 0) continue; // TODO: Warn about blanks?
    if (debug_ > 0) mprintf("DEBUG: AddOrAppend set '%s'", (*ds)->legend());
    if (isMonotonic) (*ds)->SetDim(Dimension::X, Xdim);
    DataSet* existingSet = CheckForSet( (*ds)->Meta() );
    if (existingSet == 0) {
      // New set. If set is scalar 1D but X values are not monotonic,
      // convert to XY mesh if necessary.
      if ( !isMonotonic &&
           (*ds)->Group() == DataSet::SCALAR_1D &&
           (*ds)->Type() != DataSet::XYMESH )
      {
        DataSet* xyptr = AddSet( DataSet::XYMESH, (*ds)->Meta() );
        if (xyptr == 0) {
          mprinterr("Error: Could not convert set %s to XY mesh.\n", (*ds)->legend());
          continue; // TODO exit instead?
        }
        if ( (*ds)->Size() != Xvals.size() ) { // Sanity check
          mprinterr("Error: # of X values does not match set %s size.\n", (*ds)->legend());
          continue;
        }
        DataSet_1D const& set = static_cast<DataSet_1D const&>( *(*ds) );
        DataSet_Mesh& xy = static_cast<DataSet_Mesh&>( *xyptr );
        for (unsigned int i = 0; i != set.Size(); i++)
          xy.AddXY( Xvals[i], set.Dval(i) );
        xy.SetDim(Dimension::X, Xdim);
        if (debug_ > 0) mprintf(", New set, converted to XY-MESH\n");
        // Free memory since set has now been converted.
        delete *ds;
      } else { // No conversion. Just add.
        (*ds)->SetDim(Dimension::X, Xdim);
        AddSet( *ds );
        if (debug_ > 0) mprintf(", New set\n");
      }
    } else {
      if (debug_ > 0) mprintf(", appending to existing set\n");
      // Set exists. Try to append this set to existing set.
      bool canAppend = true;
      if ( (*ds)->Group() == DataSet::GENERIC ) {
        // For GENERIC sets, base Type must match. // TODO: Should this logic be in DataSets?
        if ( (*ds)->Type() != existingSet->Type() )
          canAppend = false;
      } else {
        // For non-GENERIC sets, Group must match
        if ( (*ds)->Group() != existingSet->Group() )
          canAppend = false;
      }
      if (!canAppend)
        mprinterr("Error: Cannot append set of type %s to set of type %s\n",
                  DataSet::description((*ds)->Type()),
                  DataSet::description(existingSet->Type()));
      // If cannot append or attempt to append fails, rename and add as new set.
      if (!canAppend || existingSet->Append( *ds )) { // TODO Dimension check?
        if (canAppend)
          mprintf("Warning: Append currently not supported for type %s\n",
                  DataSet::description(existingSet->Type()));
        MetaData md = (*ds)->Meta();
        md.SetName( GenerateDefaultName("X") );
        mprintf("Warning: Renaming %s to %s\n", (*ds)->Meta().PrintName().c_str(),
                md.PrintName().c_str());
        (*ds)->SetMeta( md );
       
        AddSet( *ds );
      } else // Set was appended to existing set; memory can be freed.
        delete *ds;
    }
  } // END loop over input sets
  return 0;
}

// DataSetList::AddCopyOfSet()
void DataSetList::AddCopyOfSet(DataSet* dsetIn) {
  if (!hasCopies_ && !DataList_.empty()) {
    mprinterr("Internal Error: Attempting to add copy of DataSet (%s) to DataSetList"
              " not set up to hold copies.\n",
              dsetIn->Meta().PrintName().c_str());
    return;
  }
  hasCopies_ = true;
  Push_Back( dsetIn );
}

void DataSetList::PrintList(DataListType const& dlist) {
  size_t memTotal = 0;
  for (DataListType::const_iterator ds = dlist.begin(); ds != dlist.end(); ++ds) {
    DataSet const& dset = static_cast<DataSet const&>( *(*ds) );
    mprintf("\t%s \"%s\" (%s%s), size is %zu", dset.Meta().PrintName().c_str(), 
            dset.legend(), DataSet::description(dset.Type()),
            dset.Meta().ScalarDescription().c_str(),
            dset.Size());
    size_t memUsage = dset.MemUsageInBytes();
    if (memUsage > 0)
      mprintf(" (%s)", ByteString(memUsage, BYTE_DECIMAL).c_str());
    memTotal += memUsage;
    dset.Info();
    mprintf("\n");
  }
# ifdef MPI
  // Print sets from remaining ranks.
  if (Parallel::EnsembleIsSetup()) {
    Parallel::EnsembleComm().Barrier();
    for (int rank = 1; rank < Parallel::EnsembleComm().Size(); rank++) {
      if (rank == Parallel::EnsembleComm().Rank() && Parallel::TrajComm().Master()) {
        for (DataListType::const_iterator ds = dlist.begin(); ds != dlist.end(); ++ds) {
          DataSet const& dset = static_cast<DataSet const&>( *(*ds) );
          rprintf("%s \"%s\" (%s%s), size is %zu\n", dset.Meta().PrintName().c_str(),
                  dset.legend(), DataSet::description(dset.Type()),
                  dset.Meta().ScalarDescription().c_str(), dset.Size());
        }
      }
      mflush();
      Parallel::EnsembleComm().Barrier();
    }
    Parallel::EnsembleComm().Barrier();
  }
# endif
  mprintf("    Total data set memory usage is at least %s\n", ByteString(memTotal, BYTE_DECIMAL).c_str());
}

// DataSetList::List()
void DataSetList::List() const {
  if (!hasCopies_) { // No copies; this is a Master DSL.
    if (DataList_.empty()) return;
    mprintf("\nDATASETS (%zu total):\n", DataList_.size());
  } else if (DataList_.empty()) {
    mprintf("  No data sets.");
    return;
  }
  PrintList( DataList_ );
}

/** List all non-Topology/Reference data sets. */
void DataSetList::ListDataOnly() const {
  DataListType temp;
  for (DataListType::const_iterator ds = DataList_.begin(); ds != DataList_.end(); ++ds)
    if ( (*ds)->Type() != DataSet::REF_FRAME && (*ds)->Type() != DataSet::TOPOLOGY )
      temp.push_back( *ds );
  if (temp.empty()) return;
  mprintf("\nDATASETS (%zu total):\n", temp.size());
  PrintList( temp );
}

void DataSetList::ListStringVar() const {
  for (DataListType::const_iterator ds = DataList_.begin(); ds != DataList_.end(); ++ds)
  {
    if ( (*ds)->Type() == DataSet::STRINGVAR) {
      DataSet_StringVar const& var = static_cast<DataSet_StringVar const&>( *(*ds) );
      mprintf("\t%s = %s\n", var.legend(), var.Value().c_str());
    }
  }
}

#ifdef MPI
// DataSetList::SynchronizeData()
/** Synchronize timeseries data from child ranks to master. */
int DataSetList::SynchronizeData(Parallel::Comm const& commIn) {
  if (commIn.Size() < 2) return 0;
  // Ensure that the number of sets that require sync is same on each rank.
  // FIXME: Make sure this allgather does not end up taking too much time.
  //        Should it be debug only?
  std::vector<int> size_on_rank;
  size_on_rank.reserve( DataList_.size() );
  DataListType SetsToSync;
  for (DataListType::iterator ds = DataList_.begin(); ds != DataList_.end(); ++ds)
    if ( (*ds)->NeedsSync() ) {
      SetsToSync.push_back( *ds );
      size_on_rank.push_back( (*ds)->Size() );
    }
# ifdef PARALLEL_DEBUG_VERBOSE
  // DEBUG
  rprintf("DEBUG: SYNCING SETS\n");
  for (int rank = 0; rank != commIn.Size(); rank++) {
    if (rank == commIn.Rank()) {
      std::vector<int>::const_iterator sz = size_on_rank.begin();
      for (DataListType::const_iterator ds = SetsToSync.begin(); ds != SetsToSync.end(); ++ds, ++sz)
        rprintf("SET '%s' %i\n", (*ds)->legend(), *sz);
    }
    commIn.Barrier();
  }
# endif
  std::vector<int> n_on_rank( commIn.Size(), 0 );
  int nSets = (int)SetsToSync.size();
  commIn.AllGather( &nSets, 1, MPI_INT, &n_on_rank[0] );
  for (int rank = 1; rank < commIn.Size(); rank++)
    if (n_on_rank[rank] != n_on_rank[0]) {
      mprinterr("Internal Error: Number of sets to sync on rank %i (%i) != number on master %i\n",
                rank, n_on_rank[rank], n_on_rank[0]);
      return 1;
    }
  // Send all data set sizes to master.
  std::vector<int> all_rank_sizes;
  if (commIn.Master()) {
    all_rank_sizes.resize( nSets * commIn.Size() );
    commIn.GatherMaster( &size_on_rank[0], nSets, MPI_INT, &all_rank_sizes[0] );
  } else {
    commIn.GatherMaster( &size_on_rank[0], nSets, MPI_INT, 0 );
  }
  size_on_rank.clear();
  // Call Sync only for sets that need it.
  std::vector<int> rank_frames( commIn.Size() );
  int total = 0; //TODO size_t?
  int idx0 = 0;
  for (DataListType::iterator ds = SetsToSync.begin(); ds != SetsToSync.end(); ++ds, ++idx0) {
    if (commIn.Master()) {
      total = all_rank_sizes[idx0];
      rank_frames[0] = all_rank_sizes[idx0];
      int idx1 = idx0 + nSets;
      for (int rank = 1; rank < commIn.Size(); rank++, idx1 += nSets) {
        total += all_rank_sizes[idx1];
        rank_frames[rank] = all_rank_sizes[idx1];
      }
      //mprintf("DEBUG: Syncing '%s' (size=%zu, total=%i)\n", (*ds)->Meta().PrintName().c_str(),
      //        (*ds)->Size(), total);
    }
    if ( (*ds)->Sync(total, rank_frames, commIn) ) {
      rprintf( "Warning: Could not sync dataset '%s'\n",(*ds)->legend());
      //return;
    }
    (*ds)->SetNeedsSync( false );
  }
  return 0;
}
#endif
// -----------------------------------------------------------------------------

int DataSetList::UpdateStringVar(std::string const& varnameIn, std::string const& newValue)
const
{
  DataSet* ds = FindSetOfType( varnameIn, DataSet::STRINGVAR );
  if (ds == 0) {
    mprinterr("Error: No string variable named '%s'\n", varnameIn.c_str());
    return 1;
  }
  DataSet_StringVar& var = static_cast<DataSet_StringVar&>( *ds );
  var.assign( newValue );
  return 0;
}

/** Get a value from specified DataSet as a string. */
std::string DataSetList::GetVariable(std::string const& varnameIn) const {
  DataSet* ds = CheckForSet( varnameIn );
  if (ds == 0) return std::string();
  std::string val;
  if (ds->Size() > 0) {
    if (ds->Type() == DataSet::STRINGVAR) {
      val = ((DataSet_StringVar*)ds)->Value();
    } else if (ds->Type() == DataSet::STRING) {
      val = (*((DataSet_string*)ds))[0];
    } else if (ds->Group() == DataSet::SCALAR_1D) { 
      val = doubleToString( ((DataSet_1D*)ds)->Dval(0) );
    } else {
      mprinterr("Internal Error: DataSetList::GetVariable(): Invalid set type: '%s'\n",
                ds->legend());
      return std::string();
    }
  }
  
  return val;
}

/** Replace all variables (beginning with $) in string with their values.
  * \param varname Final string containing values.
  * \param varnameIn Initial string containing variables to replace.
  */
int DataSetList::ReplaceVariables(std::string& varname, std::string const& varnameIn)
const
{
  int nReplaced = 0;
  varname = varnameIn;
  size_t pos = varname.find("$");
  while (pos != std::string::npos) {
    // Argument is/contains a variable. Find first non-alphanumeric char
    size_t len = 1;
    for (size_t pos1 = pos+1; pos1 < varname.size(); pos1++, len++)
      if (!isalnum(varname[pos1])) break;
    std::string var_in_arg = varname.substr(pos+1, len-1);
    // String variables will never have aspect/index etc. See if string
    // variable with this name exists.
    DataSet* ds = 0;
    //mprintf("DEBUG: Check for string var: '%s'\n", var_in_arg.c_str());
    for (const_iterator it = begin(); it != end(); ++it)
    {
      if ( (*it)->Type() == DataSet::STRINGVAR && (*it)->Matches_Exact(MetaData(var_in_arg)) ) {
        ds = *it;
        break;
      }
    }
    if (ds == 0) {
      // String variable not found; see if this is a DataSet.
      for (size_t pos1 = pos+len; pos1 < varname.size(); pos1++, len++)
      {
        if (!isalnum(varname[pos1]) &&
            varname[pos1] != '[' &&
            varname[pos1] != ':' &&
            varname[pos1] != ']' &&
            varname[pos1] != '_' &&
            varname[pos1] != '-' &&
            varname[pos1] != '%')
          break;
      }
      var_in_arg = varname.substr(pos+1, len-1);
      ds = CheckForSet( var_in_arg );
    }
    if (ds == 0) {
      mprinterr("Error: Unrecognized variable in command: %s\n", var_in_arg.c_str());
      return -1;
    } else {
      if (ds->Type() != DataSet::STRINGVAR &&
          ds->Type() != DataSet::STRING && 
          ds->Group() != DataSet::SCALAR_1D)
      {
        mprinterr("Error: Only strings and 1D data sets supported for variable replacement.\n");
        return -1;
      }
      if (ds->Size() < 1) {
        mprinterr("Error: Set is empty.\n");
        return -1;
      }
      if (ds->Size() > 1)
        mprintf("Warning: Only using first value.\n");
      std::string value;
      if (ds->Type() == DataSet::STRINGVAR)
        value = (*((DataSet_StringVar*)ds)).Value(); 
      else if (ds->Type() == DataSet::STRING)
        value = (*((DataSet_string*)ds))[0];
      else
        value = doubleToString(((DataSet_1D*)ds)->Dval(0));
      if (debug_ > 0)
        mprintf("DEBUG: Replaced variable '$%s' with value '%s' from DataSet '%s'\n",
                var_in_arg.c_str(), value.c_str(), ds->legend());
      varname.replace(pos, var_in_arg.size()+1, value);
      nReplaced++;
    }
    pos = varname.find("$");
  } // END loop over this argument
  return nReplaced;
}

// -----------------------------------------------------------------------------
const char* DataSetList::RefArgs = "reference | ref <name> | refindex <#>";

/** Search for a REF_FRAME DataSet. Provided for backwards compatibility
  * with the FrameList::GetFrameFromArgs() routine.
  * The keywords in order of precedence are:
  *   - 'ref <name>'  : Get reference frame by full/base filename or tag.
  *   - 'reference'   : First reference frame in list.
  *   - 'refindex <#>': Reference frame at position.
  * \param err Set to 1 if keyword present but no reference found, 0 otherwise.
  */
DataSet* DataSetList::GetReferenceSet(ArgList& argIn, int& err) const {
  DataSet* ref = 0;
  err = 0;
  std::string refname = argIn.GetStringKey("ref");
  if (!refname.empty()) {
    ref = FindSetOfType( refname, DataSet::REF_FRAME );
    if (ref == 0) {
      mprinterr("Error: Reference '%s' not found.\n", refname.c_str());
      err = 1;
    }
  } else {
    int refindex = argIn.getKeyInt("refindex", -1);
    if (argIn.hasKey("reference")) refindex = 0;
    if (refindex > -1 && refindex < (int)RefList_.size())
      ref = RefList_[refindex];
    if (refindex != -1 && ref == 0) {
      mprinterr("Error: Reference index %i not found.\n", refindex);
      err = 1;
    }
  }
  return ref;
}

ReferenceFrame DataSetList::GetReferenceFrame(ArgList& argIn) const {
  int err = 0;
  DataSet* ds = GetReferenceSet(argIn, err);
  if (ds == 0) return ReferenceFrame(err);
  return ReferenceFrame((DataSet_Coords_REF*)ds);
}

int DataSetList::SetActiveReference(ArgList &argIn) {
  int err = 0;
  DataSet* ds = GetReferenceSet( argIn, err );
  if (ds == 0) {
    // For backwards compat, see if there is a single integer arg.
    ArgList refArg( "refindex " + argIn.GetStringNext() );
    ds = GetReferenceSet( refArg, err );
  }
  return SetActiveReference( ds );
}

int DataSetList::SetActiveReference(DataSet* dsIn) {
  if (dsIn == 0) return 1;
  activeRef_ = dsIn;
  ReferenceFrame REF((DataSet_Coords_REF*)dsIn);
  mprintf("\tSetting active reference for distance-based masks: '%s'\n", REF.refName());
  // Set in all Topologies and COORDS data sets.
  for (DataListType::const_iterator ds = DataList_.begin(); ds != DataList_.end(); ++ds)
  {
    if ( (*ds)->Type() == DataSet::TOPOLOGY )
      ((DataSet_Topology*)(*ds))->TopPtr()->SetDistMaskRef( REF.Coord() );
    else if ( (*ds)->Group() == DataSet::COORDINATES )
      ((DataSet_Coords*)(*ds))->TopPtr()->SetDistMaskRef( REF.Coord() );
  }
  return 0;
}

// DataSetList::ListReferenceFrames()
void DataSetList::ListReferenceFrames() const {
  if (!RefList_.empty()) {
    mprintf("\nREFERENCE FRAMES (%zu total):\n", RefList_.size());
    for (DataListType::const_iterator ref = RefList_.begin(); ref != RefList_.end(); ++ref)
      mprintf("    %li: %s\n", ref - RefList_.begin(), (*ref)->Meta().PrintName().c_str());
    if (activeRef_ != 0)
      mprintf("\tActive reference frame for distance-based masks is '%s'\n", activeRef_->legend());
  }
}

// DataSetList::ClearRef()
void DataSetList::ClearRef() {
  DataListType setsToKeep;
  setsToKeep.reserve( DataList_.size() - RefList_.size() );
  for (DataListType::const_iterator ds = DataList_.begin(); ds != DataList_.end(); ++ds)
    if ( (*ds)->Type() != DataSet::REF_FRAME )
      setsToKeep.push_back( *ds );
  if (!hasCopies_)
    for (DataListType::const_iterator ds = RefList_.begin(); ds != RefList_.end(); ++ds)
      delete *ds;
  RefList_.clear();
  DataList_ = setsToKeep;
}

// -----------------------------------------------------------------------------
const char* DataSetList::TopArgs = "parm <name> | parmindex <#>";

// DataSetList::GetTopByKeyword()
DataSet* DataSetList::GetTopByKeyword(ArgList& argIn, int& err) const {
  DataSet* top = 0;
  err = 0;
  std::string topname = argIn.GetStringKey("parm");
  if (!topname.empty()) {
    top = FindSetOfType( topname, DataSet::TOPOLOGY );
    if ( top == 0 ) {
      mprinterr("Error: Topology '%s' not found.\n", topname.c_str());
      err = 1;
    }
  } else {
    int topindex = argIn.getKeyInt("parmindex", -1);
    if (topindex > -1 && topindex < (int)TopList_.size())
      top = TopList_[topindex];
    if (topindex != -1 && top == 0) {
      mprinterr("Error: Topology index %i not found.\n", topindex);
      err = 1;
    }
  }
  return top;
}

// DataSetList::GetTopology()
/** \return Topology specified by a keyword. If no keywords are specified, the
  *         first loaded Topology is returned. If no Topology loaded or an
  *         error happens, 0 is returned.
  */
Topology* DataSetList::GetTopology(ArgList& argIn) const {
  if (TopList_.empty()) return 0;
  int err;
  DataSet* top = GetTopByKeyword( argIn, err );
  if (err) return 0;
  if (top == 0) // By default return first parm if nothing else specified.
    top = TopList_.front();
  return ((DataSet_Topology*)top)->TopPtr();
}

const char* DataSetList::TopIdxArgs = "parm <name> | crdset <set> | parmindex <#> | <#>";

// DataSetList::GetTopByIndex()
/** \return Topology specfied by a keyword, or if no keywords specified
  *         the Topology specified by a single integer argument (index).
  *         If no Topology loaded, return 0 and print error message.
  */
Topology* DataSetList::GetTopByIndex(ArgList& argIn) const {
  //if (TopList_.empty()) {
  //  mprinterr("Error: No Topologies are loaded.\n");
  //  return 0;
  //}
  int err;
  DataSet* top = GetTopByKeyword( argIn, err );
  if (err) return 0;
  // Check coords sets
  std::string crdset = argIn.GetStringKey("crdset");
  if (!crdset.empty()) {
    top = FindSetOfGroup(crdset, DataSet::COORDINATES);
    if ( top == 0) return 0;
    return ((DataSet_Coords*)top)->TopPtr();
  }
  if (top == 0) { // For backwards compat., check for single integer arg.
    int topindex = argIn.getNextInteger(-1);
    if (topindex > -1 && topindex < (int)TopList_.size())
      top = TopList_[topindex];
    if (topindex != -1 && top == 0) {
      mprinterr("Error: Topology index %i not found.\n", topindex);
      return 0;
    }
  }
  if (top == 0) {
    // By default return first parm if nothing else specified.
    if (TopList_.empty())
      return 0;
    top = TopList_.front();
  }
  return ((DataSet_Topology*)top)->TopPtr();
}

// DataSetList::ListTopologies()
void DataSetList::ListTopologies() const {
  if (!TopList_.empty()) {
    mprintf("\nPARAMETER FILES (%zu total):\n", TopList_.size());
    for (DataListType::const_iterator ds = TopList_.begin(); ds != TopList_.end(); ++ds)
    {
      Topology const& top = ((DataSet_Topology*)*ds)->Top();
      mprintf(" %i:", top.Pindex());
      if ( (*ds)->Meta().Name() != (*ds)->Meta().Fname().Base() )
        mprintf(" %s", (*ds)->Meta().Name().c_str());
      top.Brief(0);
      mprintf("\n");
    }
  }
}

// DataSetList::ClearTop()
void DataSetList::ClearTop() {
  DataListType setsToKeep;
  setsToKeep.reserve( DataList_.size() - TopList_.size() );
  for (DataListType::const_iterator ds = DataList_.begin(); ds != DataList_.end(); ++ds)
    if ( (*ds)->Type() != DataSet::TOPOLOGY )
      setsToKeep.push_back( *ds );
  if (!hasCopies_)
    for (DataListType::const_iterator ds = TopList_.begin(); ds != TopList_.end(); ++ds)
      delete *ds;
  TopList_.clear();
  DataList_ = setsToKeep;
}
