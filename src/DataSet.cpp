// DataSet
#include <cmath> // sqrt
#include "DataSet.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
DataSet::DataSet() :
  idx_(-1),
  dType_(UNKNOWN_DATA),
  dim_(1),
  width_(0),
  precision_(0),
  leadingSpace_(1),
  data_format_(NULL),
  scalarmode_(UNKNOWN_MODE),
  scalartype_(UNDEFINED)
{
  //fprintf(stderr,"DataSet Constructor.\n");
}

// DESTRUCTOR
DataSet::~DataSet() {
  //fprintf(stderr,"DataSet Destructor\n");
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
  * of the dataset. If Nin<=0 the dataset will be allocated dynamically.
  */
int DataSet::SetupSet(std::string const& nameIn, int Nin, int idxIn,
                      std::string const& aspectIn)
{
  // Dataset name
  if (nameIn.empty()) {
    mprintf("Dataset has no name.\n");
    return 1;
  }
  name_ = nameIn;
 
  // Attempt to allocate DataSet if necessary
  if (Nin > 0) {
    if ( Allocate( Nin ) ) return 1;
  }

  // Set index and aspect if given
  if (idxIn != -1) idx_ = idxIn;
  if (!aspectIn.empty()) aspect_ = aspectIn;
 
  return 0;
}

// DataSet::Info()
void DataSet::Info() {
  mprintf("    Data set %s",name_.c_str());
  mprintf(", size is %i\n", Size());
  //mprintf(", current is %i\n",current_);
}


// DataSet::CheckSet()
/** Return 1 if size==0, which indicates set has not been written to.
  * Otherwise return 0.
  */
int DataSet::CheckSet() {
  //if (current_==0) return 1;
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
    case FLOAT : SetDoubleFormatString(format_, width_, precision_, 1, leftAlign); break;
    case INT   : SetIntegerFormatString(format_, width_, leftAlign); break;
    case STRING: SetStringFormatString(format_, width_, leftAlign); break;
    case VECTOR: SetDoubleFormatString(format_, width_, precision_, 0, false); break;
    default:
      mprinterr("Error: No format string defined for this data type (%s).\n", c_str());
      return 1;
  }
  // Assign format to a constant ptr to avoid continuous calls to c_str
  data_format_ = format_.c_str();
  // If left aligning, no leading space. 
  if (leftAlign)  
    leadingSpace_ = 0;
  else
    leadingSpace_ = 1;
  // Set header format string
  SetStringFormatString(header_format_, width_, leftAlign);
  return 0;
}

// DataSet::WriteNameToBuffer()
/** Write the dataset name to the given character buffer. Use Sprintf
  * so that allocation happens automatically.
  */
void DataSet::WriteNameToBuffer(CharBuffer &cbuffer) {
  std::string temp_name;
  if (!aspect_.empty() && idx_ == -1)
    temp_name = name_ + aspect_;
  else if (!aspect_.empty() && idx_ != -1)
    //temp_name = aspect_ + integerToString( idx_ );
    temp_name = aspect_;
  else
    temp_name = name_;
  // If left aligning, add '#' to name; ensure that name will not be
  // larger than column width.
  if (leadingSpace_ == 0) {
    if (temp_name[0]!='#')
      temp_name.insert(0,"#");
  }
  if ((int)temp_name.size() > width_)
    temp_name.resize( width_ );
  cbuffer.Sprintf(header_format_.c_str(), temp_name.c_str());
}

// DataSet::Legend()
std::string DataSet::Legend() {
  std::string temp_name;
  if (!aspect_.empty() && idx_ == -1)
    temp_name = name_ + aspect_;
  else if (!aspect_.empty() && idx_ != -1)
    //temp_name = aspect_ + integerToString( idx_ );
    temp_name = aspect_;
  else
    temp_name = name_;
  return temp_name;
}

// DataSet::SetScalar()
void DataSet::SetScalar( scalarMode modeIn, scalarType typeIn ) {
  scalarmode_ = modeIn;
  scalartype_ = typeIn;
}

// DataSet::Matches()
bool DataSet::Matches( std::string const& dsname, int idxnum, std::string const& attr_arg )
{
  if ( dsname != name_ ) return false;
  if (idxnum != -1 && idxnum != idx_) return false;
  if (!attr_arg.empty() && attr_arg != aspect_) return false;
  return true;
}

// DataSet::Avg()
/** Calculate the average over values in this set if this set
  * is an atomic type (i.e. int, double, float).
  */
double DataSet::Avg(double *stdev) {
  double sum;

  // Check # values
  if ( Size() == 0 ) return 0;
  int numvalues = Size();
  double avg = 0;
  // Check if this set is a good type
  if (dType_==DOUBLE || 
      dType_==FLOAT ||
      dType_==INT)
  {
    sum = 0;
    for ( int i = 0; i < numvalues; ++i )  
      sum += Dval( i );
    avg = sum / (double)numvalues;
    if (stdev==NULL) return avg;

    // Stdev
    sum = 0;
    for ( int i = 0; i < numvalues; ++i ) {
      double diff = avg - Dval( i );
      diff *= diff;
      sum += diff;
    } 
    sum /= (double)numvalues;
    *stdev = sqrt(sum);
  }
  return avg;
}

// DataSet_double::Max()
/** Return the maximum value in the dataset.  */
double DataSet::Max() {
  // Check # values
  if ( Size() == 0 ) return 0;
  double max = 0;
  // Check if this set is a good type
  if (dType_==DOUBLE || 
      dType_==FLOAT ||
      dType_==INT)
  {
    max = Dval( 0 );
    for (int i = 1; i < Size(); ++i) {
      double val = Dval( i );
      if (val > max) max = val;
    } 
  }
  return max;
}

// DataSet::Min()
/** Return the minimum value in the dataset.  */
double DataSet::Min() {
  // Check # values
  if (Size()==0) return 0;
  double min = 0;
  // Check if this set is a good type
  if (dType_==DOUBLE ||
      dType_==FLOAT ||
      dType_==INT)
  { 
    min = Dval( 0 );
    for (int i = 1; i < Size(); ++i) {
      double val = Dval( i );
      if (val < min) min = val;
    } 
  }
  return min;
}

