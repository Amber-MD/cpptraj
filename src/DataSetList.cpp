// DataSetList
#include <algorithm> // sort
// This also includes basic DataSet class and dataType
#include "DataSetList.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // convertToInteger and DigitWidth
#include "ArgList.h"
#include "Range.h"
// Data types go here
#include "DataSet_double.h"
#include "DataSet_string.h"
#include "DataSet_integer.h"
#include "DataSet_float.h"
#include "DataSet_Vector.h"
#include "DataSet_Matrix.h"
#include "Histogram.h"
#include "TriangleMatrix.h"
#include "Matrix_2D.h"
#include "DataSet_Modes.h"
#include "DataSet_Coords.h"

// CONSTRUCTOR
DataSetList::DataSetList() :
  debug_(0),
  hasCopies_(false), 
  maxFrames_(0)
{}

// DESTRUCTOR
DataSetList::~DataSetList() {
  Clear();
}

void DataSetList::Clear() {
  if (!hasCopies_)
    for (DataListType::iterator ds = DataList_.begin(); ds != DataList_.end(); ++ds) 
      delete *ds;
  DataList_.clear();
  hasCopies_ = false;
  maxFrames_ = 0;
} 

DataSetList& DataSetList::operator+=(DataSetList const& rhs) {
  // It is OK if rhs does not have copies, but this should have copies.
  // For now just set hasCopies to true.
  hasCopies_ = true; 
  // Append rhs entries to here
  for (DataListType::const_iterator DS = rhs.begin(); DS != rhs.end(); ++DS)
    DataList_.push_back( *DS );
  // Update maxframes
  if (rhs.maxFrames_ > maxFrames_)
    maxFrames_ = rhs.maxFrames_;
  return *this;
}

// DataSetList::erase()
// NOTE: In order to call erase, must use iterator and not const_iterator.
//       Hence, the conversion. The new standard *should* allow const_iterator
//       to be passed to erase(), but this is currently not portable.
/** Erase element pointed to by posIn from the list. */
void DataSetList::erase( const_iterator posIn ) {
  std::vector<DataSet*>::iterator pos = DataList_.begin() + (posIn - DataList_.begin());  
  DataList_.erase( pos ); 
} 

// DataSetList::sort()
void DataSetList::sort() {
  std::sort( DataList_.begin(), DataList_.end(), dsl_cmp() );
}

// DataSetList::SetDebug()
void DataSetList::SetDebug(int debugIn) {
  debug_ = debugIn;
  if (debug_>0) 
    mprintf("DataSetList Debug Level set to %i\n",debug_);
}

// DataSetList::SetMax()
void DataSetList::SetMax(int expectedMax) {
  maxFrames_ = expectedMax;
  if (maxFrames_<0) maxFrames_ = 0;
}

/** Call Allocate for each DataSet in the list. */
void DataSetList::AllocateSets() {
  if (maxFrames_ == 0) return;
  for (DataListType::iterator ds = DataList_.begin(); ds != DataList_.end(); ++ds)
    (*ds)->Allocate( maxFrames_ );
}

/* DataSetList::SetPrecisionOfDatasets()
 * Set the width and precision for all datasets in the list.
 */
void DataSetList::SetPrecisionOfDatasets(int widthIn, int precisionIn) {
  for (DataListType::iterator ds = DataList_.begin(); ds != DataList_.end(); ++ds) 
    (*ds)->SetPrecision(widthIn,precisionIn);
}

// DataSetList::ParseArgString()
/** Separate argument nameIn specifying DataSet into name, index, and 
  * attribute parts.
  * Possible formats:
  *  - "<name>"         : Plain dataset name.
  *  - "<name>:<index>" : Dataset within larger overall set (e.g. perres:1)
  *  - "<name>[<attr>]" : Dataset with name and given attribute (e.g. rog[max])
  *  - "<name>[<attr>]:<index>" : 
  *       Dataset with name, given attribute, and index (e.g. NA[shear]:1)
  */
std::string DataSetList::ParseArgString(std::string const& nameIn, std::string& idx_arg,
                                        std::string& attr_arg)
{
  std::string dsname( nameIn );
  attr_arg.clear();
  //mprinterr("DBG: ParseArgString called with %s\n", nameIn.c_str());
  // Separate out index arg if present
  size_t idx_pos = dsname.find( ':' );
  if ( idx_pos != std::string::npos ) {
    // Advance to after the ':'
    idx_arg = dsname.substr( idx_pos + 1 );
    //mprinterr("DBG:\t\tIndex Arg [%s]\n", idx_arg.c_str());
    // Drop the index arg
    dsname.resize( idx_pos );
  }

  // Separate out aspect arg if present
  size_t attr_pos0 = dsname.find_first_of( '[' );
  size_t attr_pos1 = dsname.find_last_of( ']' );
  if ( attr_pos0 != std::string::npos && attr_pos1 != std::string::npos ) {
    if ( (attr_pos0 != std::string::npos && attr_pos1 == std::string::npos) ||
         (attr_pos0 == std::string::npos && attr_pos1 != std::string::npos) )
    {
      mprinterr("Error: Malformed attribute ([<attr>]) in dataset name %s\n", nameIn.c_str());
      return 0;
    }
    // Advance to after '[', length is position of ']' minus '[' minus 1 
    attr_arg = dsname.substr( attr_pos0 + 1, attr_pos1 - attr_pos0 - 1 );
    //mprinterr("DBG:\t\tAttr Arg [%s]\n", attr_arg.c_str());
    // Drop the attribute arg
    dsname.resize( attr_pos0 );
  }
  //mprinterr("DBG:\t\tName Arg [%s]\n", dsname.c_str());

  return dsname;
}

// DataSetList::GetDataSet()
DataSet* DataSetList::GetDataSet( std::string const& nameIn ) const {
  std::string attr_arg;
  std::string idx_arg;
  std::string dsname = ParseArgString( nameIn, idx_arg, attr_arg );
  int idx = -1;
  if (!idx_arg.empty()) idx = convertToInteger(idx_arg); // TODO: Set idx_arg to -1
  return GetSet( dsname, idx, attr_arg );
}

// DataSetList::GetMultipleSets()
/** \return a list of all DataSets matching the given argument. */
DataSetList DataSetList::GetMultipleSets( std::string const& nameIn ) const {
  DataSetList dsetOut;
  dsetOut.hasCopies_ = true;
  Range idxrange;

  std::string attr_arg;
  std::string idx_arg;
  std::string dsname = ParseArgString( nameIn, idx_arg, attr_arg );
  //mprinterr("DBG: GetMultipleSets \"%s\": Looking for %s[%s]:%s\n",nameIn.c_str(), dsname.c_str(), attr_arg.c_str(), idx_arg.c_str());
  if (idx_arg.empty())
    idxrange.SetRange( -1, 0 ); 
  else
    idxrange.SetRange( idx_arg );

  // All start selected
  std::vector<char> SelectedSets(DataList_.size(), 'T');
  // First check name
  std::vector<char>::iterator selected = SelectedSets.begin();
  if ( dsname != "*" ) {
    for (DataListType::const_iterator ds = DataList_.begin(); ds != DataList_.end(); ++ds) {
      if ( (*ds)->Name() != dsname ) *selected = 'F';
      ++selected;
    }
  }
  // Second check aspect
  if ( attr_arg != "*" ) {
    selected = SelectedSets.begin();
    for (DataListType::const_iterator ds = DataList_.begin(); ds != DataList_.end(); ++ds) {
      if ( *selected == 'T' && (*ds)->Aspect() != attr_arg ) *selected = 'F';
      ++selected;
    }
  }
  // Last check index
  if ( !idx_arg.empty() && idx_arg != "*" ) {
    selected = SelectedSets.begin();
    for (DataListType::const_iterator ds = DataList_.begin(); ds != DataList_.end(); ++ds) {
      if ( *selected == 'T' && !idxrange.InRange( (*ds)->Idx() ) ) *selected = 'F';
      ++selected;
    }
  }
  // Add selected DataSets to dsetOut
  selected = SelectedSets.begin();
  for (DataListType::const_iterator ds = DataList_.begin(); ds != DataList_.end(); ++ds)
    if ( *(selected++) == 'T' ) dsetOut.DataList_.push_back( *ds );
  return dsetOut;
}

// DataSetList::GetSet()
/** \return Specified Dataset or null if not found.
  */
DataSet* DataSetList::GetSet(std::string const& dsname, int idx, std::string const& aspect) const 
{
  for (DataListType::const_iterator ds = DataList_.begin(); ds != DataList_.end(); ++ds) 
    if ( (*ds)->Matches( dsname, idx, aspect ) ) return *ds;
  return 0;
}

// DataSetList::AddSet()
/** Add a DataSet with given name, or if no name given create a name based on 
  * defaultName and DataSet position.
  */
DataSet* DataSetList::AddSet( DataSet::DataType inType, std::string const& nameIn,
                              const char* defaultName )
{
  if (nameIn.empty()) {
    if (defaultName == 0) {
      mprinterr("Internal Error: DataSetList::AddSet() called without default name.\n");
      return 0;
    }
    return AddSetIdxAspect( inType, GenerateDefaultName(defaultName), -1, std::string() ); 
  } else
    return AddSetIdxAspect( inType, nameIn, -1, std::string() );
}

// DataSetList::GenerateDefaultName()
/** Create a name based on the given defaultName and # of DataSets,
  * i.e. defaultName_XXXXX 
  */
std::string DataSetList::GenerateDefaultName(const char* defaultName) const {
  // Determine # chars needed to hold text version of set number (min 5).
  size_t extsize = (size_t) DigitWidth( size() );
  if (extsize < 5) extsize = 5;
  return ( std::string(defaultName) + "_" + integerToString(size(), extsize) ); 
}

// DataSetList::AddSetIdx()
/** Add DataSet of specified type with given name and index to list. */
DataSet* DataSetList::AddSetIdx(DataSet::DataType inType,
                                std::string const& nameIn, int idxIn)
{
  return AddSetIdxAspect( inType, nameIn, idxIn, std::string() );
}

// DataSetList::AddSetAspect()
/** Add DataSet of specified type with given name and aspect to list. */
DataSet* DataSetList::AddSetAspect(DataSet::DataType inType,
                                   std::string const& nameIn,
                                   std::string const& aspectIn)
{
  return AddSetIdxAspect( inType, nameIn, -1, aspectIn );
}

// DataSetList::AddSetIdxAspect()
DataSet* DataSetList::AddSetIdxAspect(DataSet::DataType inType,
                                      std::string const& nameIn,
                                      int idxIn, std::string const& aspectIn,
                                      std::string const& legendIn)
{
  DataSet* ds = AddSetIdxAspect( inType, nameIn, idxIn, aspectIn );
  if (ds != 0)
    ds->SetLegend( legendIn );
  return ds;
}

// DataSetList::AddSetIdxAspect()
/** Add a DataSet of specified type, set it up and return pointer to it. 
  * \param inType type of DataSet to add.
  * \param nameIn DataSet name.
  * \param idxIn DataSet index, -1 if not specified.
  * \param aspectIn DataSet aspect, empty if not specified.
  * \param MAXin Size to set dataset to; DataSet will be set to maxFrames if < 1.
  * \return pointer to successfully set-up dataset.
  */ 
DataSet* DataSetList::AddSetIdxAspect(DataSet::DataType inType, 
                                     std::string const& nameIn, int idxIn,
                                     std::string const& aspectIn) 
{
  // Do not add to a list with copies
  if (hasCopies_) {
    mprinterr("Internal Error: Adding DataSet %s copy to invalid list.\n", nameIn.c_str());
    return 0;
  }

  // Check if DataSet with same attributes already present.
  DataSet* DS = GetSet(nameIn, idxIn, aspectIn);
  if (DS != 0) {
    mprintf("Warning: DataSet %s:%i already present.\n", nameIn.c_str(), idxIn);
    // NOTE: Should return found dataset?
    return 0; 
  }

  switch (inType) {
    case DataSet::DOUBLE       : DS = new DataSet_double(); break;
    case DataSet::FLOAT        : DS = new DataSet_float(); break;
    case DataSet::STRING       : DS = new DataSet_string(); break;
    case DataSet::INT          : DS = new DataSet_integer(); break;
    case DataSet::HIST         : DS = new Histogram(); break;
    case DataSet::TRIMATRIX    : DS = new TriangleMatrix(); break;
    case DataSet::MATRIX2D     : DS = new Matrix_2D(); break;
    case DataSet::VECTOR       : DS = new DataSet_Vector(); break;
    case DataSet::MATRIX       : DS = new DataSet_Matrix(); break;
    case DataSet::MODES        : DS = new DataSet_Modes(); break;
    case DataSet::COORDS       : DS = new DataSet_Coords(); break;
    case DataSet::UNKNOWN_DATA :
    default:
      mprinterr("Error: DataSetList::Add: Unknown set type.\n");
      return 0;
  }
  if (DS==0) {
    mprinterr("Internal Error: DataSet %s memory allocation failed.\n", nameIn.c_str());
    return 0;
  }

  // Set up dataset 
  if ( DS->SetupSet(nameIn, idxIn, aspectIn) ) {
    mprinterr("Error setting up data set %s.\n",nameIn.c_str());
    delete DS;
    return 0;
  }

  DataList_.push_back(DS); 
  //fprintf(stderr,"ADDED dataset %s\n",dsetName);
  return DS;
}

// DataSetList::AddCopyOfSet()
void DataSetList::AddCopyOfSet(DataSet* dsetIn) {
  if (!hasCopies_ && !DataList_.empty()) {
    mprinterr("Internal Error: Adding DataSet (%s) copy to invalid list\n", 
    dsetIn->Legend().c_str());
    return;
  }
  hasCopies_ = true;
  DataList_.push_back( dsetIn );
}

// DataSetList::List()
/** Print information on all data sets in the list, as well as any datafiles
  * that will be written to.
  */
void DataSetList::List() const {
  if (DataList_.empty())
    mprintf("  No data sets.");
  else if (DataList_.size()==1)
    mprintf("  1 data set: ");
  else
    mprintf("  %zu data sets: ", DataList_.size());

  mprintf("\n");
  for (unsigned int ds=0; ds<DataList_.size(); ds++) {
    mprintf("\t%s", DataList_[ds]->Name().c_str());
    if (!DataList_[ds]->Aspect().empty())
      mprintf("[%s]", DataList_[ds]->Aspect().c_str());
    if (DataList_[ds]->Idx() != -1)
      mprintf(":%i", DataList_[ds]->Idx());
    mprintf(" \"%s\"", DataList_[ds]->Legend().c_str());
    mprintf(" (%s)", DataList_[ds]->TypeName());
    mprintf(", size is %i", DataList_[ds]->Size());
    DataList_[ds]->Info();
    mprintf("\n");
  }
}

// DataSetList::Sync()
void DataSetList::Sync() {
  // Sync datasets - does nothing if worldsize is 1
  for (DataListType::iterator ds = DataList_.begin(); ds != DataList_.end(); ++ds) {
    if ( (*ds)->Sync() ) {
      rprintf( "Error syncing dataset %s\n",(*ds)->Legend().c_str());
      //return;
    }
  }
}

// DataSetList::FindSetOfType()
DataSet* DataSetList::FindSetOfType(std::string const& nameIn, DataSet::DataType typeIn) const
{
  for (DataListType::const_iterator ds = DataList_.begin(); ds != DataList_.end(); ++ds) {
    if ( (*ds)->Type() == typeIn ) {
      if ( (*ds)->Name() == nameIn )
        return (*ds);
    }
  }
  return 0;
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
    // crdset specified
    outset = FindSetOfType(setname, DataSet::COORDS);
  }
  return outset;
}
