// DataSet
#include <cmath> // sqrt
#include "DataSet.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
DataSet::DataSet() {
  //fprintf(stderr,"DataSet Constructor.\n");
  idx_ = -1;
  dType_ = UNKNOWN_DATA;
//  N=0;
  current_ = 0;
  width_ = 0;
  precision_ = 0;
  leadingSpace_ = 1;
  data_format_ = NULL;
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

// DataSet::Setup()
/** Set up common to all data sets. The dataset name should be unique and is
  * checked for in DataSetList prior to this call. Nin is the expected size 
  * of the dataset. If Nin<=0 the dataset will be allocated dynamically.
  */
int DataSet::Setup(char *nameIn, int Nin) {
  // Dataset name
  if (nameIn==NULL) {
    mprintf("Dataset has no name.\n");
    return 1;
  }
  name_.assign( nameIn );
  // Dataset memory
  //N=Nin;
  //if (N<=0) N=0;
  
  return 0;
}

// DataSet::Info()
void DataSet::Info() {
  mprintf("    Data set %s",name_.c_str());
  //mprintf(", size is %i",N);
  mprintf(", current is %i\n",current_);
}


// DataSet::CheckSet()
/** Return 1 if current==0, which indicates set has not been written to.
  * Otherwise return 0.
  */
int DataSet::CheckSet() {
  if (current_==0) return 1;
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
    case DOUBLE:
    case FLOAT : SetDoubleFormatString(format_, width_, precision_, leftAlign); break;
    case INT   : SetIntegerFormatString(format_, width_, leftAlign); break;
    case STRING: SetStringFormatString(format_, width_, leftAlign); break;
    default:
      mprinterr("Error: No format string defined for this data type (%s).\n",Name());
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
  std::string temp_name = name_;
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

// DataSet::Avg()
/** Calculate the average over values in this set if this set
  * is an atomic type (i.e. int, double, float).
  */
double DataSet::Avg(double *stdev) {
  double sum, numvalues, avg, diff;
  // Check # values
  if (current_==0) return 0;
  avg = 0;
  // Check if this set is a good type
  if (dType_==DOUBLE || 
      dType_==FLOAT ||
      dType_==INT)
  {
    sum = 0;
    numvalues = 0;
    Begin();
    do {
      sum += CurrentValue();
      ++numvalues;
    } while (NextValue());
    // NOTE: Manually calc # values since current is increment
    //       whenever Add is called even if a value is already present.
    //numvalues = (double) current;
    //mprintf("DEBUG: AVG: %.0lf values.\n",numvalues);
    avg = sum / numvalues;
    if (stdev==NULL) return avg;

    // Stdev
    sum = 0;
    Begin();
    do {
      diff = avg - CurrentValue();
      diff *= diff;
      sum += diff;
    } while(NextValue());
    sum /= numvalues;
    *stdev = sqrt(sum);
  }
  return avg;
}

// DataSet_double::Max()
/** Return the maximum value in the dataset.  */
double DataSet::Max() {
  double max;
  // Check # values
  if (current_==0) return 0;
  max = 0;
  // Check if this set is a good type
  if (dType_==DOUBLE || 
      dType_==FLOAT ||
      dType_==INT)
  {
    Begin();
    max = CurrentValue();
    do {
      double val = CurrentValue();
      if (val > max) max = val;
    } while (NextValue());
  }
  return max;
}

// DataSet::Min()
/** Return the minimum value in the dataset.  */
double DataSet::Min() {
  double min;
  // Check # values
  if (current_==0) return 0;
  min = 0;
  // Check if this set is a good type
  if (dType_==DOUBLE ||
      dType_==FLOAT ||
      dType_==INT)
  { 
    Begin();
    min = CurrentValue();
    do {
      double val = CurrentValue();
      if (val < min) min = val;
    } while (NextValue());
  }
  return min;
}

// DataSet::Capacity()
//int DataSet::Capacity() {
//  return N;
//}

// DataSet::Name()
char *DataSet::Name() { 
  return (char*)name_.c_str();
}

// DataSet::SetIdx()
void DataSet::SetIdx(int idxIn) { 
  idx_ = idxIn;
}

// DataSet::Idx()
int DataSet::Idx() { 
  return idx_;
}

// DataSet::Type()
DataSet::DataType DataSet::Type() { 
  return dType_; 
}

