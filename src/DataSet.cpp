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
    case VECTOR: SetDoubleFormatString(format_, width_, precision_, 0, false); break;
    default:
      mprinterr("Error: No format string defined for this data type (%s).\n", c_str());
      return 1;
  }
  // Assign format to a constant ptr to avoid continuous calls to c_str
  data_format_ = format_.c_str();
  return 0;
}

// DataSet::Legend()
std::string DataSet::Legend() {
  if (!legend_.empty())
   return legend_;
  else {
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

// DataSet::GoodCalcType()
bool DataSet::GoodCalcType() {
  if (dType_==DOUBLE || dType_==FLOAT || dType_==INT)
    return true;
  mprinterr("Error: DataSet %s is not a valid type for this calc.\n", c_str());
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

// DataSet::Corr()
/** Calculate correlation between DataSets.
  * \D2 DataSet to caclulate correlation to.
  * \Ct If not NULL, DataSet to store correlation function (must be DOUBLE).
  * \lagmax If Ct is not NULL, maximum lag to calculate correlation fn for.
  * \return Pearson product-moment correlation coefficient.
  */
double DataSet::Corr( DataSet& D2, DataSet* Ct, int lagmax ) {
  double d1, d2, ct;
  // Check if D1 and D2 are valid types
  if ( !GoodCalcType( )    ) return 0; 
  if ( !D2.GoodCalcType( ) ) return 0;
  // Check that D1 and D2 have same # data points.
  int Nelements = Size();
  if (Nelements != D2.Size()) {
    mprinterr("Error: Corr: # elements in dataset %s (%i) not equal to\n", c_str(), Nelements);
    mprinterr("             # elements in dataset %s (%i)\n", D2.c_str(), D2.Size());
    return 0;
  }
  // Calculate averages
  double avg1 = Avg(NULL);
  double avg2 = D2.Avg(NULL);
  // Compute normalization
  double sumdiff1_2 = 0;
  double sumdiff2_2 = 0;
  double corr_coeff = 0;
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
  double norm = sumdiff1_2 * sumdiff2_2;
  if (norm <= 0) {
    mprintf("Warning: Corr: %s to %s, Normalization sqrt <= 0.\n",Legend().c_str(),
            D2.Legend().c_str());
    return 0;
  }
  norm = sqrt( norm );
  // Correlation coefficient
  corr_coeff /= ( sqrt( sumdiff1_2 ) * sqrt( sumdiff2_2 ) );
  //mprintf("    CORRELATION COEFFICIENT %6s to %6s IS %10.4f\n",
  //        D1_->c_str(), D2_->c_str(), corr_coeff );

  if (Ct != NULL) {
    if ( Ct->Type() != DOUBLE ) {
      mprinterr("Internal Error: DataSet Corr Ct calc must be called with DataSet::DOUBLE.\n");
      return corr_coeff;
    }
    // Calculate correlation
    for (int lag = 0; lag < lagmax; lag++) {
      ct = 0;
      for (int j = 0; j < Nelements - lag; j++) {
        d1 = Dval(j);
        d2 = D2.Dval(j+lag);
        ct += ((d1 - avg1) * (d2 - avg2));
      }
      ct /= norm;
      Ct->Add(lag, &ct);
    }
  }
  return corr_coeff;
}
