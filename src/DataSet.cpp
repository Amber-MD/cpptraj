// DataSet
#include "DataSet.h"
#include "StringRoutines.h" // SetStringFormatString etc
#include "CpptrajStdio.h"

// CONSTRUCTOR
DataSet::DataSet() :
  data_format_(0),
  dType_(UNKNOWN_DATA),
  dGroup_(GENERIC),
  colwidth_(0),
  width_(0),
  precision_(0),
  leftAlign_(false)
{ }

/// CONSTRUCTOR - Take type, group, width, precision, and dimension
DataSet::DataSet(DataType typeIn, DataGroup groupIn, int widthIn, int precisionIn, int dimIn) :
  data_format_(0),
  dim_(dimIn),
  dType_(typeIn),
  dGroup_(groupIn),
  colwidth_(widthIn),
  width_(widthIn),
  precision_(precisionIn),
  leftAlign_(false)
{
  SetDataSetFormat(leftAlign_);
}

// COPY CONSTRUCTOR
DataSet::DataSet(const DataSet& rhs) :
  data_format_(0),
  dim_(rhs.dim_),
  dType_(rhs.dType_),
  dGroup_(rhs.dGroup_),
  colwidth_(rhs.colwidth_),
  width_(rhs.width_),
  precision_(rhs.precision_),
  leftAlign_(rhs.leftAlign_),
  format_(rhs.format_),
  meta_(rhs.meta_)
{
  if (!format_.empty())
    data_format_ = format_.c_str();
}

// ASSIGNMENT
DataSet& DataSet::operator=(const DataSet& rhs) {
  if (this != &rhs) {
    dim_ = rhs.dim_;
    dType_ = rhs.dType_;
    dGroup_ = rhs.dGroup_;
    colwidth_ = rhs.colwidth_;
    width_ = rhs.width_;
    precision_ = rhs.precision_;
    leftAlign_ = rhs.leftAlign_;
    format_ = rhs.format_;
    if (!format_.empty()) 
      data_format_ = format_.c_str();
  }
  return *this;
}

void DataSet::WriteCoord(CpptrajFile& file, const char* fmt, unsigned int dim, size_t pos) const {
  file.Printf( fmt, dim_[dim].Coord(pos) );
}

// DataSet::SetWidth()
/** Set only DataSet width */
void DataSet::SetWidth(int widthIn) {
  width_ = widthIn;
  SetDataSetFormat( leftAlign_ );
}

// DataSet::SetPrecision()
/** Set DataSet width and precision; recalc. output format string.
  */
void DataSet::SetPrecision(int widthIn, int precisionIn) {
  width_ = widthIn;
  precision_ = precisionIn;
  SetDataSetFormat( leftAlign_ );
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

// DataSet::SetDataSetFormat()
/** Sets the output format strings for DataSet data and name.
  * \param leftAlign if true the data and header will be left-aligned,
  *        otherwise they will be preceded by a space.
  * \return 0 on success, 1 on error.
  */
// TODO: Each data set has their own 
int DataSet::SetDataSetFormat(bool leftAlignIn) {
  leftAlign_ = leftAlignIn;
  // Set data format string.
  // NOTE: According to C++ std 4.7/4 (int)true == 1
  colwidth_ = width_ + (int)(!leftAlign_);
  switch (dType_) {
    case MODES :
    case REMLOG:
    case MATRIX_DBL:
    case XYMESH   :
    case TRAJ     :
    case REF_FRAME:
    case DOUBLE   : format_ = SetDoubleFormatString(width_, precision_, 0); break;
    case MATRIX_FLT:
    case GRID_FLT  :
    case COORDS : 
    case FLOAT  : format_ = SetDoubleFormatString(width_, precision_, 1); break;
    case INTEGER: format_ = SetIntegerFormatString(width_); break;
    case STRING : format_ = SetStringFormatString(width_, leftAlign_); break;
    case VECTOR:
      format_ = SetDoubleFormatString(width_, precision_, 0); 
      colwidth_ = (width_ + 1) * 6; // Vx Vy Vz Ox Oy Oz
      break;
    case MAT3X3:
      format_ = SetDoubleFormatString(width_, precision_, 0);
      colwidth_ = (width_ + 1) * 9;
      break;
    default:
      mprinterr("Error: No format string defined for this data type (%s).\n", 
                meta_.PrintName().c_str());
      return 1;
  }
  // If we are not left-aligning prepend a space to the format string.
  if (!leftAlign_) format_ = " " + format_;
  // Assign format to a constant ptr to avoid continuous calls to c_str
  data_format_ = format_.c_str();
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
