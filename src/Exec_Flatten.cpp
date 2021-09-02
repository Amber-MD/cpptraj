#include "Exec_Flatten.h"
#include "CpptrajStdio.h"
#include "DataSet_2D.h"
#include "DataSet_1D.h"
#include "OnlineVarT.h"
#include <algorithm> // std::max

// Exec_Flatten::Help()
void Exec_Flatten::Help() const
{
  mprintf("\tname <output set name> [mode {sum|avg}] <input set args>\n"
          "  Flatten 1 or more matrices into 1D array(s) by summing or\n"
          "  averaging elements.\n"
         );
}

// Exec_Flatten::Execute()
Exec::RetType Exec_Flatten::Execute(CpptrajState& State, ArgList& argIn)
{
  std::vector<DataSet_2D*> inpSets;
  std::vector<DataSet_1D*> outSets;

  std::string dsname = argIn.GetStringKey("name");
  if (dsname.empty()) {
    mprinterr("Error: Must specify output set name with 'name'\n");
    return CpptrajState::ERR;
  }
  mprintf("\tOutput set name: %s\n", dsname.c_str());

  const char* modeStr[] = {"sum", "avg"};
  enum ModeType { SUM = 0, AVG };
  ModeType mode;
  std::string modeArg = argIn.GetStringKey("mode");
  if (!modeArg.empty()) {
    if (modeArg == "sum")
      mode = SUM;
    else if (modeArg == "avg")
      mode = AVG;
    else {
      mprinterr("Error: Unrecognized keyword for 'mode': %s\n", modeArg.c_str());
      return CpptrajState::ERR;
    }
  } else
    mode = SUM;
  mprintf("\tFlatten mode: %s\n", modeStr[mode]);

  // Get input sets
  std::string setarg = argIn.GetStringNext();
  while (!setarg.empty()) {
    DataSetList dsl = State.DSL().GetMultipleSets( setarg );
    for (DataSetList::const_iterator ds = dsl.begin(); ds != dsl.end(); ++ds)
    {
      if ( (*ds)->Group() != DataSet::MATRIX_2D ) {
        mprintf("Warning: Set '%s' is not a matrix, skipping.\n", (*ds)->legend());
      } else {
        inpSets.push_back( (DataSet_2D*)(*ds) );
      }
    }
    setarg = argIn.GetStringNext();
  }
  mprintf("\t%zu matrices to flatten.\n", inpSets.size());
  if (inpSets.empty()) return CpptrajState::OK;

  // Set up output sets.
  if (inpSets.size() == 1) {
    DataSet* ds = State.DSL().AddSet( DataSet::DOUBLE, dsname );
    if (ds == 0) return CpptrajState::ERR;
    outSets.push_back( (DataSet_1D*)ds );
  } else {
    int idx = 1;
    for (std::vector<DataSet_2D*>::const_iterator it = inpSets.begin(); it != inpSets.end(); ++it)
    {
      MetaData md(dsname, idx);
      DataSet* ds = State.DSL().AddSet( DataSet::DOUBLE, md );
      if (ds == 0) return CpptrajState::ERR;
      outSets.push_back( (DataSet_1D*)ds );
      idx++;
    }
  }

  // Flatten
  std::vector<double> sumArray;
  std::vector< Stats<double> > avgArray;
  
  for (unsigned int idx = 0; idx != inpSets.size(); idx++) {
    DataSet_2D const& Mat = static_cast<DataSet_2D const&>( *(inpSets[idx]) );
    DataSet_1D&       Out = static_cast<DataSet_1D&      >( *(outSets[idx]) );
    // Determine max size for Out and allocate
    unsigned int maxSize = std::max(Mat.Ncols(), Mat.Nrows());
    mprintf("\tMatrix: %s, %zu columns, %zu rows. Output: %u elements\n", Mat.legend(), Mat.Ncols(), Mat.Nrows(), maxSize);
    Out.Allocate( DataSet::SizeArray(1, maxSize) );
    // Allocate temp array
    if (mode == SUM) {
      sumArray.clear();
      sumArray.resize( maxSize, 0.0 );
    } else {
      avgArray.clear();
      avgArray.resize( maxSize, Stats<double>() );
    }
    // Loop over array.
    for (unsigned int row = 0; row != Mat.Nrows(); row++) {
      for (unsigned int col = 0; col != Mat.Ncols(); col++) {
        double dval = Mat.GetElement(col, row);
        //mprintf("DBG: %8u %8u %12.4f\n", col, row, dval);
        // Diagonal element gets full interaction
        if (row == col) {
          //if (col == 0) mprintf("DBG: Diag Elt 0 %12.4f\n", dval);
          if (mode == SUM)
            sumArray[row] += dval;
          else
            avgArray[row].accumulate( dval );
        } else {
          // Off-diagonal element: divide element between row and col
          double halfVal = dval / 2.0;
          //mprintf("DBG0: %8u %12.4f\n", col, halfVal);
          //mprintf("DBG1: %8u %12.4f\n", row, halfVal);
          //if (col == 0 || row == 0) mprintf("DBG: Elt 0 %12.4f\n", halfVal);
          // Operation
          if (mode == SUM) {
            sumArray[row] += halfVal;
            sumArray[col] += halfVal;
          } else {
            avgArray[row].accumulate( halfVal );
            avgArray[col].accumulate( halfVal );
          }
        }
      } // END loop over columns
    } // END loop over rows
    // Fill in final values TODO make Resize a DataSet_1D function
    if (mode == SUM) {
      int jdx = 0;
      for (std::vector<double>::const_iterator it = sumArray.begin(); it != sumArray.end(); ++it) {
        double val = *it;
        Out.Add(jdx++, &val);
      }
    } else {
      int jdx = 0;
      for (std::vector< Stats<double> >::const_iterator it = avgArray.begin(); it != avgArray.end(); ++it) {
        double val = it->mean();
        Out.Add(jdx++, &val );
      }
    }
  } // END loop over input sets
  return CpptrajState::OK;  
}
