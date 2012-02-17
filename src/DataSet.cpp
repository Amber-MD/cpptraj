// DataSet
#include <cstdio>
#include <cmath> // sqrt
#include "DataSet.h"
#include "CpptrajStdio.h"

// SetFormatString()
/** Set up a printf-style format string with a leading space.
  * \param formatString string to be set
  * \param dType the dataType of the format string.
  * \param width the width of the data
  * \param precision the precision of the data (floating point only)
  * \param leftAlign if true do not put a leading space in format string.
  * \return 0 on success, 1 on error.
  */
int SetFormatString(std::string &formatString, dataType dType, int width, int precision, 
                    bool leftAlign) 
{
  char *format = NULL;
  int wWidth;
  int pWidth;
  size_t stringWidth;
  char leftSpace[2];
  char alignChar[2];

  // If not left-aligned, set leading space
  if (!leftAlign)
    leftSpace[0]=' ';
  else
    leftSpace[0]='\0';
  leftSpace[1]='\0';

  // Calc num of chars necessary to hold width
  wWidth = (width / 10) + 1;

  switch (dType) {
    case DOUBLE :
    case FLOAT  :
      // Calc num of chars necessary to hold precision
      pWidth = (precision / 10) + 1;
      // String fmt: " %w.plf\0"
      stringWidth = pWidth + wWidth + 6;
      format = new char [ stringWidth ];
      sprintf(format, "%s%%%i.%ilf", leftSpace, width, precision);
      break;
    case STRING :
      if (!leftAlign)
        alignChar[0]='\0';
      else
        alignChar[0]='-';
      alignChar[1]='\0';
      // String fmt: " %-ws"
      stringWidth = wWidth + 5;
      format = new char[ stringWidth ];
      sprintf(format, "%s%%%s%is", leftSpace, alignChar, width);
      break;
    case INT    :
      // String fmt: " %wi"
      stringWidth = wWidth + 4;
      format = new char[ stringWidth ];
      sprintf(format, "%s%%%ii", leftSpace, width);
      break;
    case XYZ :
      // Calc num of chars necessary to hold precision
      pWidth = (precision / 10) + 1;
      // String fmt: "%w.plf %w.plf %w.plf\0"
      stringWidth = pWidth + wWidth + 5;
      stringWidth *= 3;
      ++stringWidth;
      format = new char[ stringWidth ];
      sprintf(format, "%s%%%i.%ilf %%%i.%ilf %%%i.%ilf",leftSpace, 
              width,precision,width,precision, width,precision);
      break;
    case UNKNOWN_DATA :
      mprintf("Internal Error: SetFormatString called with unknown data type.\n");
  }

  if (format==NULL) { 
    mprintf("Error: SetFormatString: Could not allocate memory for string.\n");
    return 1;
  // DEBUG
  } else {
    formatString.assign( format );
    //mprintf("DEBUG: Format string: [%s]\n",format);
    delete[] format;
  }
  return 0;
} 
// -----------------------------------------------------------------------------

// CONSTRUCTOR
DataSet::DataSet() {
  //fprintf(stderr,"DataSet Constructor.\n");
  idx=-1;
  N=0;
  current=0;
  width = 0;
  precision = 0;
  dType = UNKNOWN_DATA;
  data_format=NULL;
  leadingSpace=1;
}

// DESTRUCTOR
DataSet::~DataSet() {
  //fprintf(stderr,"DataSet Destructor\n");
}

// DataSet::SetPrecision()
/** Set dataset width and precision and recalc output format string.
  */
void DataSet::SetPrecision(int widthIn, int precisionIn) {
  width=widthIn;
  precision=precisionIn;
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
  name.assign( nameIn );
  // Dataset memory
  N=Nin;
  if (N<=0) N=0;
  
  return 0;
}

// DataSet::Info()
void DataSet::Info() {
  mprintf("    Data set %s",name.c_str());
  mprintf(", size is %i",N);
  mprintf(", current is %i\n",current);
}

// DataSet::WriteNameToBuffer()
/** Write the dataset name to the given character buffer.
  */
void DataSet::WriteNameToBuffer(CharBuffer &cbuffer) {
  cbuffer.WriteString(header_format.c_str(), name.c_str());
}

// DataSet::CheckSet()
/** Return 1 if current==0, which indicates set has not been written to.
  * Otherwise return 0.
  * Call setFormatString; mostly just needed for string data sets which
  * have variable width based on the size of the strings that have been
  * stored.
  */
int DataSet::CheckSet() {
  if (current==0) return 1;
  if (SetDataSetFormat(false)) return 1;
  //mprinterr("Dataset %s has format [%s]\n",name,format);
  return 0;
}

// DataSet::SetDataSetFormat()
/** Sets the output format strings for DataSet data and name.
  * \param leftAlign if true the data and header will be left-aligned,
  *        otherwise they will be preceded by a space.
  * \return 0 on success, 1 on error.
  */
int DataSet::SetDataSetFormat(bool leftAlign) {
  if (SetFormatString(format, dType, width, precision, leftAlign)) return 1;
  data_format = format.c_str();
  // If left aligning, add '#' to name. Ensure that name will not overflow
  if (leftAlign) { 
    if (name[0]!='#') name.insert(0, "#");
    leadingSpace = 0;
  } else
    leadingSpace = 1;
  if ((int)name.size() > width) name.resize( width );
  if (SetFormatString(header_format, STRING, width, 0, leftAlign)) return 1;
  return 0;
}

// DataSet::Avg()
/** Calculate the average over values in this set if this set
  * is an atomic type (i.e. int, double, float).
  */
double DataSet::Avg(double *stdev) {
  double sum, numvalues, avg, diff;
  // Check # values
  if (current==0) return 0;
  avg = 0;
  // Check if this set is a good type
  if (dType==DOUBLE || 
      dType==FLOAT ||
      dType==INT)
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
  if (current==0) return 0;
  max = 0;
  // Check if this set is a good type
  if (dType==DOUBLE || 
      dType==FLOAT ||
      dType==INT)
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
  if (current==0) return 0;
  min = 0;
  // Check if this set is a good type
  if (dType==DOUBLE ||
      dType==FLOAT ||
      dType==INT)
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
 
