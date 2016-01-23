#include "DataSet.h"
#include "CpptrajStdio.h"

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

// ASSIGNMENT
DataSet& DataSet::operator=(const DataSet& rhs) {
  if (this != &rhs) {
    format_ = rhs.format_;
    dim_ = rhs.dim_;
    dType_ = rhs.dType_;
    dGroup_ = rhs.dGroup_;
    meta_ = rhs.meta_;
    associatedData_.clear();
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
