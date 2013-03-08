// DataSet
#include "DataSet.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // SetStringFormatString etc

// CONSTRUCTOR
DataSet::DataSet() :
  idx_(-1),
  dType_(UNKNOWN_DATA),
  dim_(1),
  width_(0),
  precision_(0),
  data_format_(0),
  scalarmode_(UNKNOWN_MODE),
  scalartype_(UNDEFINED)
{ }

const char* DataSet::SetStrings[] = {
  "unknown",   "double", "string",    "integer", "float", "vector",
  "matrix",    "modes",  "histogram", "upper-triangle matrix",
  "2D matrix", "coords"
};

/// CONSTRUCTOR - Take type, width, precision, and dimension
DataSet::DataSet(DataType typeIn, int widthIn, int precisionIn, int dimIn) :
  idx_(-1),
  dType_(typeIn),
  dim_(dimIn),
  width_(widthIn),
  precision_(precisionIn),
  data_format_(0),
  scalarmode_(UNKNOWN_MODE),
  scalartype_(UNDEFINED)
{
  SetDataSetFormat(false);
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
 
  return 0;
}

// DataSet::Empty()
/** \return true if size==0, which indicates set has not been written to. 
  * \return false otherwise.
  */
bool DataSet::Empty() {
  if (Size()==0) return 1;
  return 0;
}

// DataSet::SetDataSetFormat()
/** Sets the output format strings for DataSet data and name.
  * \param leftAlign if true the data and header will be left-aligned,
  *        otherwise they will be preceded by a space.
  * \return 0 on success, 1 on error.
  */
int DataSet::SetDataSetFormat(bool leftAlign) {
  // Set data format string
  switch (dType_) {
    case HIST  :
    case MATRIX2D:
    case DOUBLE: SetDoubleFormatString(format_, width_, precision_, 0, leftAlign); break;
    case TRIMATRIX:
    case COORDS:
    case FLOAT : SetDoubleFormatString(format_, width_, precision_, 1, leftAlign); break;
    case INT   : SetIntegerFormatString(format_, width_, leftAlign); break;
    case STRING: SetStringFormatString(format_, width_, leftAlign); break;
    case MODES :
    case MATRIX:
    case VECTOR: SetDoubleFormatString(format_, width_, precision_, 0, false); break;
    default:
      mprinterr("Error: No format string defined for this data type (%s).\n", 
                Legend().c_str());
      return 1;
  }
  // Assign format to a constant ptr to avoid continuous calls to c_str
  data_format_ = format_.c_str();
  return 0;
}

// DataSet::Legend()
/** Return DataSet legend. If the legend is empty create one based on 
  * DataSet name (and aspect/index if present). Possible formats are:
  * - Name[Aspect]
  * - Name:idx
  * - Aspect:Idx
  * - Name
  */
std::string const& DataSet::Legend() {
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
  return legend_;
}

// DataSet::SetScalar()
void DataSet::SetScalar( scalarMode modeIn, scalarType typeIn ) {
  scalarmode_ = modeIn;
  scalartype_ = typeIn;
}

// DataSet::Matches()
bool DataSet::Matches( std::string const& dsname, int idxnum, std::string const& attr_arg )
{
  /*mprintf("DEBUG: Input: %s[%s]:%i  This Set: %s[%s]:%i\n",
          dsname.c_str(), attr_arg.c_str(), idxnum, 
          name_.c_str(), aspect_.c_str(), idx_);*/
  if ( dsname != name_ && dsname != "*") return false;
  // Currently match any index if not specified.
  if (idxnum != -1 && idxnum != idx_) return false;
  // If aspect specified make sure it matches. 
  if (!attr_arg.empty() && (attr_arg != aspect_ && attr_arg != "*")) return false;
  // If no aspect specified but dataset has aspect do not match.
  if (attr_arg.empty() && !aspect_.empty()) return false;
  //mprintf("\tMATCH\n");
  return true;
}
