#include "Analysis_Regression.h"
#include "CpptrajStdio.h"
#include "DataSet_Mesh.h"

Analysis_Regression::Analysis_Regression() {}

void Analysis_Regression::Help() {
  mprintf("\t<dset0> [<dset1> ...]\n"
          "  Calculate linear regression lines for given data sets.\n");
}

// Analysis_Regression::Setup()
Analysis::RetType Analysis_Regression::Setup(ArgList& analyzeArgs, DataSetList* datasetlist,
                            TopologyList* PFLin, DataFileList* DFLin, int debugIn)
{
  //std::string outname = analyzeArgs.GetStringKey("out");
  // Select datasets from remaining args
  if (input_dsets_.AddSetsFromArgs( analyzeArgs.RemainingArgs(), *datasetlist )) {
    mprinterr("Error: Could not add data sets.\n");
    return Analysis::ERR;
  }
  if (input_dsets_.empty()) {
    mprinterr("Error: No input data sets.\n");
    return Analysis::ERR;
  }

  mprintf("    REGRESSION: Calculating linear regression of %i data sets.\n",
          input_dsets_.size());
  //if (!outname.empty())
  //  mprintf("\tWriting results to %s\n", outname.c_str());
  //for (Array1D::const_iterator set = input_dsets_.begin(); set != input_dsets_.end(); ++set)
  //  mprintf("\t%s\n", (*set)->Legend().c_str());
  //if (outfile_.OpenWrite( outname )) return Analysis::ERR;

  return Analysis::OK;
}

// Analysis_Regression::Analyze()
Analysis::RetType Analysis_Regression::Analyze() {
  int nerr = 0;
  //outfile_.Printf("#SetNum\tAverage\tStdev\tMin\tMax\tName\n");
  for (Array1D::const_iterator DS = input_dsets_.begin();
                               DS != input_dsets_.end(); ++DS)
  {
    if ( (*DS)->Size() < 2)
      mprintf("Warning: Set \"%s\" does not have enough data for regression (%zu points).\n", 
              (*DS)->Legend().c_str(), (*DS)->Size());
    else {
      DataSet_Mesh mesh, output;
      // Set XY mesh
      mesh.SetMeshXY( *(*DS) );
      mprintf("  %zu: %s\n", DS - input_dsets_.begin(), (*DS)->Legend().c_str());
      nerr += mesh.LinearRegression( output );
    }
  }
  if (nerr > 0) return Analysis::ERR;
  return Analysis::OK;
}
