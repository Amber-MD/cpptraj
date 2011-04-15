// DataSet
#include <cstdlib>
#include <cstring>
#include <cstdio> // sprintf
#include "DataSet.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
DataSet::DataSet() {
  //fprintf(stderr,"DataSet Constructor.\n");
  name=NULL;
  N=0;
  isDynamic=false;
  current=0;
  width = 0;
  precision = 0;
  format = NULL;
  dType = UNKNOWN_DATA;
}

// DESTRUCTOR
DataSet::~DataSet() {
  //fprintf(stderr,"DataSet Destructor\n");
  if (name!=NULL) free(name);
  if (format!=NULL) free(format);
}

/*
 * DataSet::setFormatString()
 * Set up the output format string for each data element based on the given 
 * dataType and the current width, and precision.
 */
void DataSet::setFormatString() {
  size_t stringWidth = 0;
  int wWidth = 0;
  int pWidth = 0;

  if (format!=NULL) {free(format); format=NULL;}

  // Calc num of chars necessary to hold width
  wWidth = (width / 10) + 1;

  switch (dType) {
    case DOUBLE :
      // Calc num of chars necessary to hold precision
      pWidth = (precision / 10) + 1;
      // String fmt: " %w.plf\0"
      stringWidth = pWidth + wWidth + 6;
      format = (char*) malloc( stringWidth * sizeof(char) );
      sprintf(format, " %%%i.%ilf", width, precision);
      break;
    case STRING :
      // String fmt: " %s"
      format = (char*) malloc( 4 * sizeof(char) );
      strcpy(format, " %s");
      break;
    case INT :
      // String fmt: " %wi"
      stringWidth = wWidth + 4;
      format = (char*) malloc( stringWidth * sizeof(char) );
      sprintf(format, " %%%ii", width);
      break;
    case UNKNOWN_DATA :
      mprintf("Internal Error: setFormatString called with unknown data type.\n");
  }

  if (format==NULL) 
    mprintf("Error: setFormatString: Could not allocate memory for string.\n");
  // DEBUG
  //else
  //  mprintf("DEBUG: Format string: [%s]\n",format);
}    

/*
 * DataSet::SetPrecision()
 * Set dataset width and precision and recalc output format string.
 */
void DataSet::SetPrecision(int widthIn, int precisionIn) {
  width=widthIn;
  precision=precisionIn;
  setFormatString();
}

/* 
 * DataSet::Setup()
 * Set up common to all data sets. The dataset name should be unique and is
 * checked for in DataSetList prior to this call. Nin is the expected size 
 * of the dataset. If Nin<=0 the dataset will be allocated dynamically.
 */
int DataSet::Setup(char *nameIn, int Nin) {
  // Dataset name
  if (nameIn==NULL) {
    mprintf("Dataset has no name.\n");
    return 1;
  }
  name=(char*) malloc( (strlen(nameIn)+1) * sizeof(char));
  strcpy(name, nameIn);
  // Dataset memory
  N=Nin;
  if (N<=0) {
    isDynamic=true;
    N=0;
  }
  if ( this->Allocate() ) return 1;
  return 0;
}

/*
 * DataSet::Info()
 * Print dataset information.
 */
void DataSet::Info() {
  mprintf("    Data set %s",name);
  mprintf(", size is ");
  if (isDynamic)
    mprintf("dynamic");
  else
    mprintf("%i",N);
  mprintf(", current is %i",current);
  mprintf(".\n");
}

/*
 * DataSet::CheckSet()
 * Return 1 if current==0, which indicates set has not been written to.
 * Otherwise return 0;
 */
int DataSet::CheckSet() {
  if (current==0) return 1;
  return 0;
}

