// DataSet
#include "DataSet.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // SetStringFormatString etc

// CONSTRUCTOR
DataSet::DataSet() :
  data_format_(0),
  idx_(-1),
  dType_(UNKNOWN_DATA),
//  dim_(0),
  colwidth_(0),
  width_(0),
  precision_(0),
  scalarmode_(UNKNOWN_MODE),
  scalartype_(UNDEFINED)
{ }

// CONSTRUCTOR
DataSet::DataSet(DataType typeIn, int widthIn, int precisionIn, int dimIn) :
  data_format_(0),
  idx_(-1),
  dType_(typeIn),
  dim_(dimIn),
  colwidth_(widthIn),
  width_(widthIn),
  precision_(precisionIn),
  scalarmode_(UNKNOWN_MODE),
  scalartype_(UNDEFINED)
{
  SetDataSetFormat(false);
  // Allocate default values for dimensions
/*
  for (unsigned int d = 0; d < dim_.size(); ++d) {
    switch (d) {
      case 0: dim_[d].SetLabel("X"); break;
      case 1: dim_[d].SetLabel("Y"); break;
      case 2: dim_[d].SetLabel("Z"); break;
      default: dim_[d].SetLabel("D" + integerToString(d));
    }
    dim_[d].SetStep(1.0);
  }
*/   
}  

// COPY CONSTRUCTOR
DataSet::DataSet(const DataSet& rhs) :
  data_format_(0),
  name_(rhs.name_),
  idx_(rhs.idx_),
  aspect_(rhs.aspect_),
  legend_(rhs.legend_),
  dType_(rhs.dType_),
//  dim_(rhs.dim_),
  colwidth_(rhs.colwidth_),
  width_(rhs.width_),
  precision_(rhs.precision_),
  format_(rhs.format_),
  scalarmode_(rhs.scalarmode_),
  scalartype_(rhs.scalartype_)
{
  if (!format_.empty())
    data_format_ = format_.c_str();
}

// ASSIGNMENT
DataSet& DataSet::operator=(const DataSet& rhs) {
  if (this == &rhs) return *this;
  name_ = rhs.name_;
  idx_ = rhs.idx_;
  aspect_ = rhs.aspect_;
  legend_ = rhs.legend_;
  dType_ = rhs.dType_;
//  dim_ = rhs.dim_;
  colwidth_ = rhs.colwidth_;
  width_ = rhs.width_;
  precision_ = rhs.precision_;
  format_ = rhs.format_;
  if (!format_.empty()) 
    data_format_ = format_.c_str();
  scalarmode_ = rhs.scalarmode_;
  scalartype_ = rhs.scalartype_;
  return *this;
}

// DataSet::SetPrecision()
/** Set dataset width and precision and recalc output format string.
  */
void DataSet::SetPrecision(int widthIn, int precisionIn) {
  width_ = widthIn;
  precision_ = precisionIn;
  SetDataSetFormat(false);
}


// DataSet::SetupSet()
/** Set up common to all data sets. The dataset name should be unique and is
  * checked for in DataSetList prior to this call. Nin is the expected size 
  * of the dataset. 
  */
int DataSet::SetupSet(std::string const& nameIn, int idxIn, std::string const& aspectIn)
{
  // Dataset name
  if (nameIn.empty()) {
    mprintf("Dataset has no name.\n");
    return 1;
  }
  name_ = nameIn;
  // Set index and aspect if given
  if (idxIn != -1) idx_ = idxIn;
  if (!aspectIn.empty()) aspect_ = aspectIn;
  // If no legend set yet create a default one. Possible formats are:
  //  - Name[Aspect]
  //  - Name:idx
  //  - Aspect:Idx
  //  - Name
  if (legend_.empty()) {
    if (!aspect_.empty() && idx_ == -1)
      legend_ = name_ + "[" + aspect_ + "]";
    else if (aspect_.empty() && idx_ != -1)
      legend_ = name_ + ":" + integerToString( idx_ );
    else if (!aspect_.empty() && idx_ != -1)
      legend_ = aspect_ + ":" + integerToString( idx_ );
    else
      legend_ = name_;
  }
  return 0;
}

// DataSet::SetDataSetFormat()
/** Sets the output format strings for DataSet data and name.
  * \param leftAlign if true the data and header will be left-aligned,
  *        otherwise they will be preceded by a space.
  * \return 0 on success, 1 on error.
  */
int DataSet::SetDataSetFormat(bool leftAlign) {
  // Set data format string.
  // NOTE: According to C++ std 4.7/4 (int)true == 1
  colwidth_ = width_ + (int)leftAlign;
  switch (dType_) {
//    case HIST  :
//    case MATRIX2D:
//    case MATRIX_VEC3:
    case MATRIX_DBL:
    case DOUBLE : format_ = SetDoubleFormatString(width_, precision_, 0, leftAlign); break;
//    case TRIMATRIX:
    case MATRIX_FLT:
    case COORDS : 
    case FLOAT  : format_ = SetDoubleFormatString(width_, precision_, 1, leftAlign); break;
    case INTEGER: format_ = SetIntegerFormatString(width_, leftAlign); break;
    case STRING : format_ = SetStringFormatString(width_, leftAlign); break;
    case MODES :
//    case MATRIX:
    case VECTOR: // No left-align allowed for now with VECTOR.
      format_ = SetDoubleFormatString(width_, precision_, 0, false); 
      colwidth_ = (width_ + 1) * 6; // Vx Vy Vz Ox Oy Oz
      break;
    default:
      mprinterr("Error: No format string defined for this data type (%s).\n", 
                Legend().c_str());
      return 1;
  }
  // Assign format to a constant ptr to avoid continuous calls to c_str
  data_format_ = format_.c_str();
  return 0;
}

// DataSet::Matches()
bool DataSet::Matches( std::string const& dsname, int idxnum, std::string const& aspect )
{
  /*mprintf("DEBUG: Input: %s[%s]:%i  This Set: %s[%s]:%i\n",
          dsname.c_str(), aspect.c_str(), idxnum, 
          name_.c_str(), aspect_.c_str(), idx_);*/
  if ( dsname != name_ && dsname != "*") return false;
  // Currently match any index if not specified.
  if (idxnum != -1 && idxnum != idx_) return false;
  // If aspect specified make sure it matches. 
  if (!aspect.empty() && (aspect != aspect_ && aspect != "*")) return false;
  // If no aspect specified but dataset has aspect do not match.
  if (aspect.empty() && !aspect_.empty()) return false;
  //mprintf("\tMATCH\n");
  return true;
}
