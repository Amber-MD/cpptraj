#include "Exec_CrdTransform.h"
#include "CpptrajStdio.h"
#include "Constants.h"
#include <algorithm> // std::max,min
#include <cmath> // floor

/// Get the minimum and maximum coordinates in a given frame, store in min and max
static inline void getMaxMin(Frame const& frmIn, Vec3& max, Vec3& min) {
  for (int at = 0; at != frmIn.Natom(); at++) {
    const double* xyz = frmIn.XYZ(at);
    for (int i = 0; i != 3; i++) {
      max[i] = std::max(max[i], xyz[i]);
      min[i] = std::min(min[i], xyz[i]);
    }
  }
}

/** Normalize coordinates between 0 and 1. */
int Exec_CrdTransform::normalizeCoords(DataSet_Coords* crdIn,
                                       DataSet_Coords* crdOut)
const
{
  mprintf("\tNormalize coordinates between 0 and 1.\n");
  mprintf("\tInput coords: %s\n", crdIn->legend());
  mprintf("\tOutput coords: %s\n", crdIn->legend());
  // Get max and min X Y and Z
  Frame frmIn = crdIn->AllocateFrame();
  crdIn->GetFrame(0, frmIn);
  Vec3 xyzmax( frmIn.XYZ(0) );
  Vec3 xyzmin( frmIn.XYZ(0) );
  getMaxMin(frmIn, xyzmax, xyzmin);
  for (unsigned int idx = 1; idx < crdIn->Size(); idx++) {
    crdIn->GetFrame(idx, frmIn);
    getMaxMin(frmIn, xyzmax, xyzmin);
  }
  mprintf("\tMax: %g %g %g\n", xyzmax[0], xyzmax[1], xyzmax[2]);
  mprintf("\tMin: %g %g %g\n", xyzmin[0], xyzmin[0], xyzmin[0]);
  // Choose the overall max and min
  double max = xyzmax[0];
  max = std::max(max, xyzmax[1]);
  max = std::max(max, xyzmax[2]);
  double min = xyzmin[0];
  min = std::min(min, xyzmin[1]);
  min = std::min(min, xyzmin[2]);

  Vec3 norm( max - min );
  // Protect against bad values
  bool hasBadValues = false;
  static const char dirStr[] = { 'X', 'Y', 'Z' };
  for (int ii = 0; ii != 3; ii++) {
    if (norm[ii] < 0) {
      mprinterr("Error: Min value > max value for %c coordinate.\n", dirStr[ii]);
      hasBadValues = true;
    }
    if (norm[ii] < Constants::SMALL) {
      mprinterr("Error: Min value == max value for %c coordinate.\n", dirStr[ii]);
      hasBadValues = true;
    }
  }
  if (hasBadValues) return 1;
  // Transform coords between 0 and 1
  for (unsigned int idx = 0; idx < crdIn->Size(); idx++) {
    crdIn->GetFrame(idx, frmIn);
    for (int crdidx = 0; crdidx < frmIn.size(); crdidx+=3) {
      frmIn[crdidx  ] = (frmIn[crdidx  ] - min) / norm[0];
      frmIn[crdidx+1] = (frmIn[crdidx+1] - min) / norm[1];
      frmIn[crdidx+2] = (frmIn[crdidx+2] - min) / norm[2];
    }
    crdOut->SetCRD(idx, frmIn);
  }

  return 0;
}

/** Transform coordinates by RMS-fitting to an average structure, calculating
  * a new average, then RMS-fitting to that average and so on until a
  * tolerance is reached. Essentially the procedure described by 
  * Klem et al. J. Chem. Theory Comput. 2022, 18, 3218−3230.
  */
int Exec_CrdTransform::iterativeRmsRefinement(AtomMask const& maskIn,
                                              bool useMass,
                                              double tolIn,
                                              DataSet_Coords* crdIn,
                                              DataSet_Coords* crdOut)
const
{
  mprintf("\tRMS iterative refinement.\n");
  mprintf("\tInput coords: %s\n", crdIn->legend());
  mprintf("\tOutput coords: %s\n", crdIn->legend());
  mprintf("\tAtom mask: %s\n", maskIn.MaskString());
  mprintf("\tRMS Tolerance: %g Ang.\n", tolIn);
  if (useMass)
    mprintf("\tMass-weighting on.\n");
  else
    mprintf("\tMass-weighting off.\n");
  // Do the initial fit to the first frame.
  Frame frmIn = crdIn->AllocateFrame();
  crdIn->GetFrame(0, frmIn);
  Frame selectedRef;
  selectedRef.SetupFrameFromMask( maskIn, crdIn->Top().Atoms() );
  selectedRef.SetCoordinates( frmIn, maskIn );
  // Ensure reference is centered on the origin
  Vec3 refTrans = selectedRef.CenterOnOrigin( useMass );
  // Set up frame for selected incoming atoms
  Frame selectedTgt = selectedRef;
  // Set up frame to hold average 
  Frame avgFrm = selectedTgt;

  double currentTol = tolIn + 9999.0;

  unsigned int iteration = 0;
  while (currentTol > tolIn) {
    avgFrm.ZeroCoords();
    Vec3 tgtTrans(0.0);
    Matrix_3x3 rot(0.0);
    for (unsigned int idx = 0; idx != crdIn->Size(); idx++) {
      crdIn->GetFrame(idx, frmIn);
      selectedTgt.SetCoordinates( frmIn, maskIn );
      selectedTgt.RMSD_CenteredRef( selectedRef, rot, tgtTrans, useMass );
      frmIn.Trans_Rot_Trans(tgtTrans, rot, refTrans);
      crdOut->SetCRD(idx, frmIn);
      avgFrm.AddByMask( frmIn, maskIn );
    }
    avgFrm.Divide( (double)crdIn->Size() );
    // Calc RMS of current average to current reference
    currentTol = avgFrm.RMSD_CenteredRef( selectedRef, rot, tgtTrans, useMass );
    mprintf("\t%8u %12.4f\n", iteration+1, currentTol);
    // Fit the current average TODO is this necessary?
    avgFrm.Trans_Rot_Trans(tgtTrans, rot, refTrans);
    // Set current average to be new reference
    selectedRef = avgFrm;
    iteration++;
  }

  return 0;
}

//const char* Exec_CrdTransform::TrimMetricStr_[] = {
//  "MSD", "RR", "JT", "SM", "No metric"
//};

const char* Exec_CrdTransform::CriterionStr_[] = {
  "comp_sim", "sim_to_medioid", "No criterion"
};

static inline void printDarray(std::vector<double> const& array) {
  mprintf("[");
  for (std::vector<double>::const_iterator it = array.begin(); it != array.end(); ++it)
    mprintf(" %f", *it);
  mprintf("]\n");
}

/** Trim a desired percentage of outliers (most dissimilar) from the COORDS
  * data set by calculating the largest complement similarity.
  */
int Exec_CrdTransform::trimOutliers(int n_trimmed, double cutoffIn,
                                    ExtendedSimilarity::MetricType metric,
                                    CriterionType criterion,
                                    DataSet_Coords* crdIn,
                                    DataSet_Coords* crdOut)
const
{
  mprintf("\tTrimming outliers.\n");
  mprintf("\tUsing metric: %s\n", ExtendedSimilarity::metricStr(metric));
  mprintf("\tCriterion: %s\n", CriterionStr_[criterion]);
  unsigned int Ncoords = crdIn->Top().Natom() * 3;
  unsigned int Nframes = crdIn->Size();
  mprintf("\t'%s' has %u coordinates, %u frames.\n", crdIn->legend(), Ncoords, Nframes);
  // Specify n_trimmed or cutoff, but not both.
  if (n_trimmed < 0 && cutoffIn < 0) {
    mprinterr("Internal Error: Must specify either number to trim or cutoff.\n");
    return 1;
  }
  if (n_trimmed >= 0 && cutoffIn > 0) {
    mprinterr("Error: Must specify either number to trim or cutoff, but not both.\n");
    return 1;
  }
  int cutoff;
  if (n_trimmed >= 0) {
    cutoff = n_trimmed;
    mprintf("\t# to trim: %i\n", n_trimmed);
  } else {
    cutoff = (int)(floor(Nframes * cutoffIn));
    mprintf("\tFraction of outliers to remove: %f\n", cutoffIn);
  }
  mprintf("\tUsing cutoff value: %i\n", cutoff);

  if (criterion == COMP_SIM) {
    // Comp sim
    std::vector<double> c_sum( Ncoords, 0.0 );
    std::vector<double> sq_sum_total( Ncoords, 0.0 );
    Frame frmIn = crdIn->AllocateFrame();
    // Get sum and sum squares for each coordinate
    for (unsigned int idx = 0; idx < crdIn->Size(); idx++) {
      crdIn->GetFrame(idx, frmIn);
      for (unsigned int icrd = 0; icrd < Ncoords; icrd++) {
        c_sum[icrd] += frmIn[icrd];
        sq_sum_total[icrd] += frmIn[icrd] * frmIn[icrd];
      }
    }
    //printDarray(sq_sum_total);
    // For each frame, get the comp. similarity
    std::vector<double> c_arr(Ncoords, 0.0);
    std::vector<double> sq_arr(Ncoords, 0.0);
    ExtendedSimilarity ExtSim;
    for (unsigned int idx = 0; idx < crdIn->Size(); idx++) {
      crdIn->GetFrame(idx, frmIn);
      for (unsigned int icrd = 0; icrd < Ncoords; icrd++) {
        c_arr[icrd]  = c_sum[icrd] - frmIn[icrd];
        sq_arr[icrd] = sq_sum_total[icrd] - (frmIn[icrd]*frmIn[icrd]);
      }
      //mprintf("%u\n", Nframes-1);
      //printDarray(sq_arr);
      //printDarray(c_sum);
      //mprintf("%i\n", crdIn->Top().Natom());
      double val = ExtSim.Comparison(c_arr, sq_arr, metric, Nframes-1, crdIn->Top().Natom());
      mprintf("DBG:\t%u %.10f\n", idx, val);
    }
    
  }

  return 0;
}


// Exec_CrdTransform::Help()
void Exec_CrdTransform::Help() const
{
  mprintf("\t<crd set>\n"
          "\t{ rmsrefine [mask <mask>] [mass] [rmstol <tolerance>] |\n"
          "\t  normcoords |\n"
          "\t  trim\n"
          "\t}\n");
}

// Exec_CrdTransform::Execute()
Exec::RetType Exec_CrdTransform::Execute(CpptrajState& State, ArgList& argIn)
{
  AtomMask mask;
  bool useMass = false;
  double rmsTol = -1.0;
  int n_trimmed = -1;
  double cutoff = -1.0;
  // Determine mode
  enum ModeType { RMSREFINE = 0, NORMCOORDS, TRIM, UNSPECIFIED };
  ModeType mode = UNSPECIFIED;
  if (argIn.hasKey("rmsrefine")) {
    mode = RMSREFINE;
    mask.SetMaskString( argIn.GetStringKey("mask") );
    useMass = argIn.hasKey("mass");
    rmsTol = argIn.getKeyDouble("rmstol", 0.0001);
  } else if (argIn.hasKey("normcoords")) {
    mode = NORMCOORDS;
  } else if (argIn.hasKey("trim")) {
    mode = TRIM;
    n_trimmed = argIn.getKeyInt("ntrimmed", -1);
    cutoff = argIn.getKeyDouble("cutoff", -1.0);
  } else {
    mprinterr("Error: Expected 'trim', 'rmsrefine', or 'normcoords'\n");
    return CpptrajState::ERR;
  }
 
  // Get COORDS set
  std::string setname = argIn.GetStringNext();
  if (setname.empty()) {
    mprinterr("Error: %s: Specify COORDS dataset name.\n", argIn.Command());
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
  if (CRD->Type() == DataSet::TRAJ) {
    mprinterr("Error: TRAJ sets not yet supported.\n"); // FIXME
    return CpptrajState::ERR;
  }

  // Set up mask
  if (mask.MaskStringSet()) {
    if (CRD->Top().SetupIntegerMask( mask )) {
      mprinterr("Error: Could not set up mask.\n");
      return CpptrajState::ERR;
    }
    mask.MaskInfo();
  }

  int err = 0;
  switch (mode) {
    case RMSREFINE  : err = iterativeRmsRefinement(mask, useMass, rmsTol, CRD, CRD); break;
    case NORMCOORDS : err = normalizeCoords(CRD, CRD); break;
    // TODO pass in criterion and metric
    case TRIM       : err = trimOutliers(n_trimmed, cutoff, ExtendedSimilarity::MSD, COMP_SIM, CRD, CRD); break;
    default         : err = 1; break;
  }
  if (err != 0) return CpptrajState::ERR;

  return CpptrajState::OK;
}
