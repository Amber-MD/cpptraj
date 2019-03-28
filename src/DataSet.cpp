#include "DataSet.h"
#include "CpptrajStdio.h"

// NOTE: IT IS IMPORTANT THAT THIS ARRAY CORRESPOND TO DataSet::DataType
/** Description of each DataType. */
const char* DataSet::Descriptions_[] = {
  "unknown",        // UNKNOWN_DATA
  "double",         // DOUBLE
  "float",          // FLOAT
  "integer",        // INTEGER
  "string",         // STRING
  "double matrix",  // MATRIX_DBL
  "float matrix",   // MATRIX_FLT
  "coordinates",    // COORDS
  "vector",         // VECTOR
  "eigenmodes",     // MODES
  "float grid",     // GRID_FLT
  "double grid",    // GRID_DBL
  "remlog",         // REMLOG
  "X-Y mesh",       // XYMESH
  "trajectories",   // TRAJ
  "reference",      // REF_FRAME
  "3x3 matrices",   // MAT3X3
  "topology",       // TOPOLOGY
  "pH",                         // PH
  "pH REMD (explicit)",         // PH_EXPL
  "pH REMD (implicit)",         // PH_IMPL
  "parameters",                 // PARAMETERS
  "pairwise matrix (mem)",      // PMATRIX_MEM
  "pairwise matrix (NetCDF)"    // PMATRIX_NC
};

// CONSTRUCTOR
DataSet::DataSet() : dType_(UNKNOWN_DATA), dGroup_(GENERIC)
# ifdef MPI
, needsSync_( false )
# endif
{}

/// CONSTRUCTOR - Take type, group, width, precision, and dimension
DataSet::DataSet(DataType typeIn, DataGroup groupIn, TextFormat const& fmtIn, int dimIn) :
  format_(fmtIn),
  dim_(dimIn, Dimension(1.0, 1.0)), // default min=1.0, step=1.0
  dType_(typeIn),
  dGroup_(groupIn)
# ifdef MPI
, needsSync_( false )
# endif
{ }

// COPY CONSTRUCTOR
DataSet::DataSet(const DataSet& rhs) :
  format_(rhs.format_),
  dim_(rhs.dim_),
  dType_(rhs.dType_),
  dGroup_(rhs.dGroup_),
  meta_(rhs.meta_)
# ifdef MPI
, needsSync_(rhs.needsSync_)
# endif
{
  for (AdataArray::const_iterator a = rhs.associatedData_.begin();
                                  a != rhs.associatedData_.end(); ++a)
    associatedData_.push_back( (*a)->Copy() );
}

/** Clear all associated data. */
void DataSet::ClearAssociatedData() {
  for (AdataArray::iterator ad = associatedData_.begin(); ad != associatedData_.end(); ++ad)
    delete *ad;
  associatedData_.clear();
}

// DESTRUCTOR
DataSet::~DataSet() { ClearAssociatedData(); }

// ASSIGNMENT
DataSet& DataSet::operator=(const DataSet& rhs) {
  if (this != &rhs) {
    format_ = rhs.format_;
    dim_ = rhs.dim_;
    dType_ = rhs.dType_;
    dGroup_ = rhs.dGroup_;
    meta_ = rhs.meta_;
    ClearAssociatedData(); 
#   ifdef MPI
    needsSync_ = rhs.needsSync_;
#   endif
    for (AdataArray::const_iterator a = rhs.associatedData_.begin();
                                    a != rhs.associatedData_.end(); ++a)
      associatedData_.push_back( (*a)->Copy() );
  }
  return *this;
}

/** Set up DataSet name, index and/or aspect etc. Also create
  * default legend if not already set.
  * \param In The DataSet meta data 
  */
int DataSet::SetMeta(MetaData const& In) {
  // Dataset name
  if (In.Name().empty()) {
    mprinterr("Internal Error: DataSet has no name.\n"); //FIXME allow?
    return 1;
  }
  meta_ = In;
  if (meta_.Legend().empty())
    meta_.SetDefaultLegend();
  return 0;
}

/** This version allows wildcards and ranges. */
bool DataSet::Matches_WC(MetaData::SearchString const& search, DataType typeIn) const
{
  // Match type if specified
  if ( typeIn != UNKNOWN_DATA && typeIn != dType_ ) return false;
  return (meta_.Match_WildCard( search ));
}

AssociatedData* DataSet::GetAssociatedData(AssociatedData::AssociatedType typeIn) const {
  for (AdataArray::const_iterator ad = associatedData_.begin();
                                  ad != associatedData_.end(); ++ad)
    if ((*ad)->Type() == typeIn) return *ad;
  return 0;
}
