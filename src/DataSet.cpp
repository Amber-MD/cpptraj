// DataSet
#include <cmath> // sqrt
#include <cstring> // memset
#include "DataSet.h"
#include "CpptrajStdio.h"
#include "PubFFT.h"

// CONSTRUCTOR
DataSet::DataSet() :
  idx_(-1),
  dType_(UNKNOWN_DATA),
  dim_(1),
  width_(0),
  precision_(0),
  data_format_(NULL),
  scalarmode_(UNKNOWN_MODE),
  scalartype_(UNDEFINED)
{
  //fprintf(stderr,"DataSet Constructor.\n");
}

/// CONSTRUCTOR - Take type, width, precision, and dimension
DataSet::DataSet(DataType typeIn, int widthIn, int precisionIn, int dimIn) :
  idx_(-1),
  dType_(typeIn),
  dim_(dimIn),
  width_(widthIn),
  precision_(precisionIn),
  data_format_(NULL),
  scalarmode_(UNKNOWN_MODE),
  scalartype_(UNDEFINED)
{
  SetDataSetFormat(false);
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

// DataSet::Info()
void DataSet::Info() {
  mprintf("\tData set %s",name_.c_str());
  mprintf(", size is %i", Size());
  mprintf(" [%s]", aspect_.c_str());
  mprintf(":%i \"%s\"\n", idx_,Legend().c_str());
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
  * - Aspect:Idx
  * - Name
  */
std::string const& DataSet::Legend() {
  if (legend_.empty()) {
    if (!aspect_.empty() && idx_ == -1)
      legend_ = name_ + "[" + aspect_ + "]";
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
  if ( dsname != name_ ) return false;
  // Currently match any index if not specified.
  if (idxnum != -1 && idxnum != idx_) return false;
  // If aspect specified make sure it matches. 
  if (!attr_arg.empty() && attr_arg != aspect_ ) return false;
  // If no aspect specified but dataset has aspect do not match.
  if (attr_arg.empty() && !aspect_.empty()) return false;
  //mprintf("\tMATCH\n");
  return true;
}

// DataSet::GoodCalcType()
bool DataSet::GoodCalcType() {
  if (dType_==DOUBLE || dType_==FLOAT || dType_==INT)
    return true;
  mprinterr("Error: DataSet %s is not a valid type for this calc.\n", 
            Legend().c_str());
  return false;
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
  if ( GoodCalcType( ) ) {
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
  if ( GoodCalcType( ) ) {
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
  if ( GoodCalcType( ) ) {
    min = Dval( 0 );
    for (int i = 1; i < Size(); ++i) {
      double val = Dval( i );
      if (val < min) min = val;
    } 
  }
  return min;
}

// DataSet::CrossCorr()
/** Calculate time correlation between two DataSets.
  * \D2 DataSet to calculate correlation to.
  * \Ct DataSet to store time correlation fn, must be DOUBLE.
  * \lagmaxIn Max lag to calculate corr. -1 means use size of dataset.
  * \calccovar If true calculate covariance (devation from avg).
  * \return 0 on success, 1 on error.
  */
int DataSet::CrossCorr( DataSet& D2, DataSet& Ct, int lagmaxIn, bool calccovar, bool usefft ) 
{
  int lagmax;
  double ct;
  // Check if D1 and D2 are valid types
  if ( !GoodCalcType( )    ) return 1;
  if ( !D2.GoodCalcType( ) ) return 1;
  // Check that D1 and D2 have same # data points.
  int Nelements = Size();
  if (Nelements != D2.Size()) {
    mprinterr("Error: CrossCorr: # elements in dataset %s (%i) not equal to\n", Legend().c_str(), 
              Nelements);
    mprinterr("                  # elements in dataset %s (%i)\n", D2.Legend().c_str(), 
              D2.Size());
    return 1;
  }
  if (Nelements < 2) {
    mprinterr("Error: CrossCorr: # elements is less than 2 (%i)\n", Nelements);
    return 1;
  }
  // Check return dataset type
  if ( Ct.Type() != DOUBLE ) {
    mprinterr("Internal Error: CrossCorr: Ct must be of type DataSet::DOUBLE.\n");
    return 1;
  }
  // Check if lagmaxIn makes sense. Set default lag to be Nelements 
  // if not specified.
  if (lagmaxIn == -1)
    lagmax = Nelements;
  else if (lagmaxIn > Nelements) {
    mprintf("Warning: CrossCorr [%s][%s]: max lag (%i) > Nelements (%i), setting to Nelements.\n",
            Legend().c_str(), D2.Legend().c_str(), lagmaxIn, Nelements);
    lagmax = Nelements;
  } else 
    lagmax = lagmaxIn;
  // If calculating covariance calculate averages
  double avg1 = 0;
  double avg2 = 0;
  if ( calccovar ) {
    avg1 = Avg(NULL);
    avg2 = D2.Avg(NULL);
  }
  // Calculate correlation
  double norm = 1.0;

  if ( usefft ) {
    // Calc using FFT
    PubFFT pubfft1(Nelements);
    int ndata = pubfft1.size() * 2; // Allocate space for real + img component
    double* data1 = new double[ ndata ];
    memset( data1, 0, ndata * sizeof(double) );
    for (int i = 0; i < Nelements; ++i) 
      data1[i*2] = Dval(i) - avg1;
    if (&D2 == this) 
      pubfft1.CorF_FFT(ndata, data1, NULL);
    else {
      // Populate second dataset if different
      double* data2 = new double[ ndata ];
      memset( data2, 0, ndata * sizeof(double) );
      for (int i = 0; i < Nelements; ++i)
        data2[i*2] = D2.Dval(i) - avg2;
      pubfft1.CorF_FFT(ndata, data1, data2);
      delete[] data2;
    }
    // Put real components of data1 in output DataSet
    norm = 1.0 / fabs( data1[0] );
    for (int i = 0; i < lagmax; ++i) {
      ct = data1[i*2] * norm;
      Ct.Add(i, &ct);
    }
    delete[] data1;
  } else {
    // Direct calc
    for (int lag = 0; lag < lagmax; ++lag) {
      ct = 0;
      int jmax = Nelements - lag;
      for (int j = 0; j < jmax; ++j) 
        ct += ((Dval(j) - avg1) * (D2.Dval(j+lag) - avg2));
      if (lag == 0) {
        if (ct != 0)
          norm = fabs( ct );
      }
      ct /= norm;
      Ct.Add(lag, &ct);
    }
  }
  return 0;
}

// DataSet::Corr()
/** Calculate Pearson product-moment correlation between DataSets.
  * \D2 DataSet to caclulate correlation to.
  * \return Pearson product-moment correlation coefficient.
  */
double DataSet::Corr( DataSet& D2 ) {
  double d1, d2;
  // Check if D1 and D2 are valid types
  if ( !GoodCalcType( )    ) return 0; 
  if ( !D2.GoodCalcType( ) ) return 0;
  // Check that D1 and D2 have same # data points.
  int Nelements = Size();
  if (Nelements != D2.Size()) {
    mprinterr("Error: Corr: # elements in dataset %s (%i) not equal to\n", 
              Legend().c_str(), Nelements);
    mprinterr("             # elements in dataset %s (%i)\n", 
              D2.Legend().c_str(), D2.Size());
    return 0;
  }
  // Calculate averages
  double avg1 = Avg(NULL);
  double avg2 = D2.Avg(NULL);
  // Calculate average deviations. 
  double sumdiff1_2 = 0.0;
  double sumdiff2_2 = 0.0;
  double corr_coeff = 0.0;
  //mprinterr("DATASETS %s and %s\n", c_str(), D2.c_str());
  for (int i = 0; i < Nelements; i++) {
    d1 = Dval(i);
    d2 = D2.Dval(i);
    //mprinterr("\t\t\t%f\t%f\n",d1,d2);
    double diff1 = d1 - avg1;
    double diff2 = d2 - avg2;
    sumdiff1_2 += (diff1 * diff1);
    sumdiff2_2 += (diff2 * diff2);
    corr_coeff += (diff1 * diff2);
  }
  if (sumdiff1_2 == 0.0 || sumdiff2_2 == 0.0) {
    mprintf("Warning: Corr: %s to %s, Normalization is 0\n",
            Legend().c_str(),  D2.Legend().c_str());
    return 0;
  }
  // Correlation coefficient
  corr_coeff /= ( sqrt( sumdiff1_2 ) * sqrt( sumdiff2_2 ) );
  //mprintf("    CORRELATION COEFFICIENT %6s to %6s IS %10.4f\n",
  //        D1_->c_str(), D2_->c_str(), corr_coeff );
  return corr_coeff;
}
