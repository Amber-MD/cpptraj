// DataSet
#include "DataSet.h"
#include "StringRoutines.h" // WildcardMatch
#include "CpptrajStdio.h"

// CONSTRUCTOR
DataSet::DataSet() : dType_(UNKNOWN_DATA), dGroup_(GENERIC) {}

/// CONSTRUCTOR - Take type, group, width, precision, and dimension
DataSet::DataSet(DataType typeIn, DataGroup groupIn, TextFormat const& fmtIn, int dimIn) :
  format_(fmtIn),
  dim_(dimIn),
  dType_(typeIn),
  dGroup_(groupIn)
{ }

// COPY CONSTRUCTOR
DataSet::DataSet(const DataSet& rhs) :
  format_(rhs.format_),
  dim_(rhs.dim_),
  dType_(rhs.dType_),
  dGroup_(rhs.dGroup_),
  meta_(rhs.meta_)
{
  for (AdataArray::const_iterator a = rhs.associatedData_.begin();
                                  a != rhs.associatedData_.end(); ++a)
    associatedData_.push_back( (*a)->Copy() );
}

// ASSIGNMENT
DataSet& DataSet::operator=(const DataSet& rhs) {
  if (this != &rhs) {
    format_ = rhs.format_;
    dim_ = rhs.dim_;
    dType_ = rhs.dType_;
    dGroup_ = rhs.dGroup_;
    meta_ = rhs.meta_;
    associatedData_.clear();
    for (AdataArray::const_iterator a = rhs.associatedData_.begin();
                                    a != rhs.associatedData_.end(); ++a)
      associatedData_.push_back( (*a)->Copy() );
  }
  return *this;
}

void DataSet::WriteCoord(CpptrajFile& file, const char* fmt, unsigned int dim, size_t pos) const {
  file.Printf( fmt, dim_[dim].Coord(pos) );
}

/** Set up DataSet name, index and/or aspect etc. Also create
  * default legend if not already set.
  * \param In The DataSet meta data 
  */
int DataSet::SetMetaData(MetaData const& In) {
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

// DataSet::Matches()
/** This version allows wildcards and ranges. */
bool DataSet::Matches_WC(std::string const& dsname, Range const& idxRange,
                         std::string const& aspect,
                         Range const& memberRange, DataType typeIn) const
{
  //mprintf("DEBUG: Input: %s[%s]:%s  This Set: %s[%s]:%i\n",
  //        dsname.c_str(), aspect.c_str(), idxRange.RangeArg(), 
  //        name_.c_str(), aspect_.c_str(), idx_);
  // Match type if specified
  if ( typeIn != UNKNOWN_DATA && typeIn != dType_ ) return false;
  // Match name
  if ( WildcardMatch(dsname, meta_.Name()) == 0 ) return false;
  // If aspect specified make sure it matches.
  if ( WildcardMatch(aspect, meta_.Aspect()) == 0 ) return false;
  // Currently match any index if not specified.
  if (idxRange.Front() != -1 && !idxRange.InRange(meta_.Idx())) return false;
  // Match any ensemble if not specified
  if (memberRange.Front() != -1 && !memberRange.InRange(meta_.EnsembleNum())) return false;
  // If no aspect specified but dataset has aspect do not match.
  //if (aspect.empty() && !aspect_.empty()) return false;
  //mprintf("\tMATCH\n");
  return true;
}
