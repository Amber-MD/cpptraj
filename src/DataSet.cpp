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

// DataSet::Matches()
/** This version allows wildcards and ranges. */
bool DataSet::Matches_WC(SearchString const& search, DataType typeIn) const
{
  //mprintf("DEBUG: Input: %s[%s]:%s  This Set: %s[%s]:%i\n",
  //        dsname.c_str(), aspect.c_str(), idxRange.RangeArg(), 
  //        name_.c_str(), aspect_.c_str(), idx_);
  // Match type if specified
  if ( typeIn != UNKNOWN_DATA && typeIn != dType_ ) return false;
  if (meta_.Fname().empty()) {
    // No filename. Match name if specified.
    if ( WildcardMatch(search.NameArg(), meta_.Name()) == 0 ) return false;
  } else {
    // Has associated file name. Exit if name specified no match with
    // name/base name/full path.
    if ( WildcardMatch(search.NameArg(), meta_.Name()) == 0 &&
         !meta_.Fname().MatchFullOrBase( search.NameArg() ) ) return false;
  }
  // If aspect specified make sure it matches.
  if ( WildcardMatch(search.AspectArg(), meta_.Aspect()) == 0 ) return false;
  // Currently match any index if not specified.
  if (search.IdxRange().Front() != -1 && !search.IdxRange().InRange(meta_.Idx()))
    return false;
  // Match any ensemble if not specified
  if (search.MemberRange().Front() != -1 && !search.MemberRange().InRange(meta_.EnsembleNum()))
    return false;
  // If no aspect specified but dataset has aspect do not match.
  //if (aspect.empty() && !aspect_.empty()) return false;
  //mprintf("\tMATCH\n");
  return true;
}

AssociatedData* DataSet::GetAssociatedData(AssociatedData::AssociatedType typeIn) const {
  for (AdataArray::const_iterator ad = associatedData_.begin();
                                  ad != associatedData_.end(); ++ad)
    if ((*ad)->Type() == typeIn) return *ad;
  return 0;
}

// -----------------------------------------------------------------------------
/** Separate argument nameIn specifying DataSet into name, index, and 
  * attribute parts.
  * Possible formats:
  *  - "<name>"         : Plain dataset name.
  *  - "<name>:<index>" : Dataset within larger overall set (e.g. perres:1)
  *  - "<name>[<attr>]" : Dataset with name and given attribute (e.g. rog[max])
  *  - "<name>[<attr>]:<index>" : 
  *       Dataset with name, given attribute, and index (e.g. NA[shear]:1)
  */
int DataSet::SearchString::ParseArgString(std::string const& argIn)
{
  name_arg_.assign( argIn );
  aspect_arg_.clear();
  idx_range_.Clear();
  member_range_.Clear();
  std::string idx_arg, member_arg;
  //mprinterr("DBG: ParseArgString called with %s\n", nameIn.c_str());
  // Separate out ensemble arg if present
  size_t idx_pos = name_arg_.find( '%' );
  if ( idx_pos != std::string::npos ) {
    // Advance to after the '.'
    member_arg = name_arg_.substr( idx_pos + 1 );
    // Drop the ensemble arg
    name_arg_.resize( idx_pos );
  }
  // Separate out index arg if present
  idx_pos = name_arg_.find( ':' );
  if ( idx_pos != std::string::npos ) {
    // Advance to after the ':'
    idx_arg = name_arg_.substr( idx_pos + 1 );
    //mprinterr("DBG:\t\tIndex Arg [%s]\n", idx_arg.c_str());
    // Drop the index arg
    name_arg_.resize( idx_pos );
  }
  // Separate out aspect arg if present
  size_t attr_pos0 = name_arg_.find_first_of( '[' );
  size_t attr_pos1 = name_arg_.find_last_of( ']' );
  if ( attr_pos0 != std::string::npos && attr_pos1 != std::string::npos ) {
    if ( (attr_pos0 != std::string::npos && attr_pos1 == std::string::npos) ||
         (attr_pos0 == std::string::npos && attr_pos1 != std::string::npos) )
    {
      mprinterr("Error: Malformed attribute ([<attr>]) in dataset arg %s\n", argIn.c_str());
      return 1;
    }
    // If '[' is at very beginning, this is a tag. Otherwise attribute.
    if (attr_pos0 > 0) {
      // Advance to after '[', length is position of ']' minus '[' minus 1 
      aspect_arg_ = name_arg_.substr( attr_pos0 + 1, attr_pos1 - attr_pos0 - 1 );
      //mprinterr("DBG:\t\tAttr Arg [%s]\n", aspect_arg_.c_str());
      // Drop the attribute arg
      name_arg_.resize( attr_pos0 );
    }
  }
  //mprinterr("DBG:\t\tName Arg [%s]\n", name_arg_.c_str());
  // If index arg is empty make wildcard (-1)
  if (idx_arg.empty() || idx_arg == "*")
    idx_range_.SetRange( -1, 0 );
  else
    idx_range_.SetRange( idx_arg );
  // If member arg is empty make none (-1)
  if (member_arg.empty() || member_arg == "*")
    member_range_.SetRange( -1, 0 );
  else
    member_range_.SetRange( member_arg );
  // If aspect arg not set and name is wildcard, make attribute wildcard.
  if (aspect_arg_.empty() && name_arg_ == "*")
    aspect_arg_.assign("*");

  return 0;
}
