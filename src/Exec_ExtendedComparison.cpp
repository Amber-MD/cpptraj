#include "Exec_ExtendedComparison.h"
#include "CpptrajStdio.h"
#include "ExtendedSimilarity.h"

// Exec_ExtendedComparison::Help()
void Exec_ExtendedComparison::Help() const
{

}

// Exec_ExtendedComparison::Execute()
Exec::RetType Exec_ExtendedComparison::Execute(CpptrajState& State, ArgList& argIn)
{
  // Metric
  std::string mstr = argIn.GetStringKey("metric");
  ExtendedSimilarity::MetricType metric = ExtendedSimilarity::NO_METRIC;
  if (!mstr.empty()) {
    metric = ExtendedSimilarity::TypeFromKeyword( mstr );
    if (metric == ExtendedSimilarity::NO_METRIC) {
      mprinterr("Error: Metric '%s' not recognized.\n", mstr.c_str());
      return CpptrajState::ERR;
    }
  } else {
    metric = ExtendedSimilarity::MSD;
  }
  // Get output set
  std::string outname = argIn.GetStringKey("name");
  DataSet* out = State.DSL().AddSet( DataSet::DOUBLE, outname, "EXTCOMP" );
  if (out == 0) {
    mprinterr("Error: Could not set up output data set for extended comparison.\n");
    return CpptrajState::ERR;
  }
  // Get COORDS set
  std::string setname = argIn.GetStringNext();
  if (setname.empty()) {
    mprinterr("Error: %s: Specify input COORDS dataset name.\n", argIn.Command());
    return CpptrajState::ERR;
  }
  DataSet_Coords* CRD = (DataSet_Coords*)State.DSL().FindSetOfGroup( setname, DataSet::COORDINATES );
  if (CRD == 0) {
    mprinterr("Error: %s: No COORDS set with name %s found.\n", argIn.Command(), setname.c_str());
    return CpptrajState::ERR;
  }
  mprintf("\tUsing set '%s'\n", CRD->legend());
  if (CRD->Size() < 1) {
    mprinterr("Error: Set '%s' has no frames.\n", CRD->legend());
    return CpptrajState::ERR;
  }
  out->SetDim(Dimension::X, Dimension(1.0, 1.0, "Frame") ); // TODO make time an option?
  out->Allocate( DataSet::SizeArray(1, CRD->Size()) );

  mprintf("\tCalculating extended comparison similarity values.\n");
  mprintf("\tInput coords: %s\n", CRD->legend());
  mprintf("\tOutput set: %s\n", out->legend());
  mprintf("\tUsing metric: %s\n", ExtendedSimilarity::metricStr(metric));

  // Do extended similarity calculation for each frame
  // TODO have ExtendedSimilarity return a DataSet_double?
  ExtendedSimilarity ExtSim;
  if (ExtSim.SetOpts( metric, CRD->Size(), CRD->Top().Natom()*3 )) {
    mprinterr("Error: Extended similarity setup failed.\n");
    return CpptrajState::ERR;
  }
  ExtendedSimilarity::Darray csimvals = ExtSim.CalculateCompSim( *CRD );
  if (csimvals.empty()) {
    mprinterr("Error: No comparitive similarity values calculated.\n");
    return CpptrajState::ERR;
  }
  const double* dptr = &csimvals[0];
  for (unsigned int idx = 0; idx != csimvals.size(); idx++, dptr++)
    out->Add(idx, dptr);

  return CpptrajState::OK;
}
