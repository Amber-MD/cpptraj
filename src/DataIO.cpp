#include "DataIO.h"
#include "StringRoutines.h"
#include "CpptrajStdio.h"
#include "Constants.h"
#include "DataSet_Mesh.h"

// DataIO::CheckValidFor()
bool DataIO::CheckValidFor( DataSet const& dataIn ) const {
  if (valid1d_ && dataIn.Ndim() == 1) return true; 
  if (valid2d_ && dataIn.Ndim() == 2) return true; 
  if (valid3d_ && dataIn.Ndim() == 3) return true;
  for (std::vector<DataSet::DataType>::const_iterator dt = valid_.begin(); 
                                                      dt != valid_.end(); ++dt)
    if (dataIn.Type() == *dt) return true;
  return false;
}

// DataIO::SetupCoordFormat()
std::string DataIO::SetupCoordFormat(size_t maxFrames, Dimension const& dim, 
                                     int default_width, int default_precision)
{
  int col_precision = default_precision;
  // Determine maximum coordinate.
  double maxCoord = (dim.Step() * (double)maxFrames) + dim.Min();
  // Determine character width necessary to hold largest coordinate.
  int col_width = DigitWidth( (long int)maxCoord );
  // Check if the precision is enough to support the step size.
  if (dim.Step() < 1.0) {
    int prec_exp_width = FloatWidth( dim.Step() );
    if (prec_exp_width > col_precision)
      col_precision = prec_exp_width;
  }
  // If the width for the column plus the characters needed for precision
  // (plus 1 for decimal point) would be greated than default_width, increment 
  // the column width by (precision+1).
  if (col_precision != 0) {
    int precision_width = col_width + col_precision + 1;
    if (precision_width > default_width) col_width = precision_width;
  }
  // Default width for column is at least default_width.
  if (col_width < default_width) col_width = default_width;
  // Set column data format string, left-aligned (no leading space).
  return SetDoubleFormatString( col_width, col_precision, 0 );
}

/// For backwards compat. FIXME
Dimension DataIO::DetermineXdim( std::vector<double> const& Xvals ) {
  int nerr;
  return DetermineXdim( Xvals, nerr );
}

/** Given X values, try to determine step size etc */
Dimension DataIO::DetermineXdim( std::vector<double> const& Xvals, int& nerr ) {
  nerr = 0;
  if ( Xvals.empty() ) {
    mprinterr("Error: Cannot determine X dimension - no X values.\n");
    nerr = 1;
    return Dimension();
  }
  if ( Xvals.size() == 1)
    return Dimension(Xvals.front(), 1.0, 1);
  // Determine from max/min
  double xstep = (Xvals.back() - Xvals.front()) / (double)(Xvals.size() - 1);
  // Check if xstep is reasonable. Check only the first few values.
  double xval = Xvals.front();
  for (unsigned int i = 0; i < std::min(10U, (unsigned int)Xvals.size()); i++) {
    if (xval - Xvals[i] > Constants::SMALL) {
      //mprintf("Warning: X dimension step may be off. Xval[%u] is %f, expected %f\n",
      //        i, Xvals[i], xval);
      ++nerr;
    }
    xval += xstep;
  }
  if (nerr > 0) mprintf("Warning: Could not accurately determine X dimension step size.\n");
  return Dimension(Xvals.front(), xstep, Xvals.size());
}

// TODO: DataSetList should have a function that allows you to add/append sets
int DataIO::AddSetsToList(DataSetList& datasetlist, Xarray const& TimeVals,
                          ArrayDD const& Sets, std::string const& dsname)
{
  int dim_err;
  DataSet::DataType Dtype;
  Dimension Xdim = DetermineXdim( TimeVals, dim_err );
  if (dim_err == 0)
    Dtype = DataSet::DOUBLE;
  else {
    mprintf("Warning: %s data sets will be X-Y mesh.\n",dsname.c_str());
    Dtype = DataSet::XYMESH;
  }
  // ----- ADD NON-EMPTY DATA SETS -----
  for (ArrayDD::const_iterator set = Sets.begin(); set != Sets.end(); ++set)
  {
    if (set->Size() > 0) {
      DataSet* ds = datasetlist.AddSetIdxAspect( Dtype, set->Name(), set->Idx(),
                                                 set->Aspect(), set->Legend() );
      if (ds == 0)
        mprinterr("Error: Could not create set for %s \"%s\"\n", set->Name().c_str(), set->Legend().c_str());
      else {
        ds->SetDim(Dimension::X, Xdim);
        if (Dtype == DataSet::DOUBLE) {
          DataSet_double& dsD = static_cast<DataSet_double&>( *ds );
          dsD = set->Data();
        } else { // DataSet::XYMESH
          DataSet_Mesh& dsM = static_cast<DataSet_Mesh&>( *ds );
          dsM.SetMeshXY( TimeVals, set->Data() );
        }
      }
    }
  }
  return 0;
}
