#include "Minimize_SteepestDescent.h"
#include <algorithm> // std::fill
#include <cmath> // sqrt
#include "PotentialFunction.h"
#include "Frame.h"
#include "Trajout_Single.h"
#include "CpptrajStdio.h"
#include "ArgList.h"
#include "DataSetList.h"
#include "CpptrajFile.h"

/** CONSTRUCTOR */
Minimize_SteepestDescent::Minimize_SteepestDescent() :
  min_tol_(1.0E-5),
  dx0_(0.01),
  nMinSteps_(1)
{}

/** Set up minimization. */
int Minimize_SteepestDescent::SetupMin(std::string const& nameIn, double tolIn, double dx0In,
                                       int stepsIn)
{
  trajoutName_ = nameIn;
  min_tol_ = tolIn;
  dx0_ = dx0In;
  nMinSteps_ = stepsIn;
  return 0;
}

/** Run minimization. */
int Minimize_SteepestDescent::RunMin(PotentialFunction& potential, Frame& frameIn,
                                     CpptrajFile& outfile)
const
{
  // Output trajectory
  int iteration = 0;
  Trajout_Single trajOut; // TODO change type
  if (!trajoutName_.empty()) {
    if (trajOut.InitTrajWrite(trajoutName_, ArgList(), DataSetList(), TrajectoryFile::UNKNOWN_TRAJ))
      return 1;
    if (trajOut.SetupTrajWrite(potential.CurrentTop(), frameIn.CoordsInfo(), 0))
      return 1;
    if (trajOut.WriteSingle(iteration, frameIn)) return 1;
  }
  
  // Zero forces
  if (!frameIn.HasForce()) {
    mprinterr("Internal Error: Frame not set up with forces.\n");
    return 1;
  }
  frameIn.ZeroForces();

  // Degrees of freedom
  double deg_of_freedom = potential.DegreesOfFreedom(); 
  double fnq = sqrt(deg_of_freedom);
  // Main loop for steepest descent
  const double dxstm = 1.0E-5;
  const double crits = 1.0E-6;
  double rms = 1.0;
  double dxst = dx0_;
  double last_e = 0.0;
  outfile.Printf("%-8s %12s %12s", "#Iter.", "ENE", "RMS");
  // Do not print labels when only a single term.
  potential.Energy().PrintActiveLabels(outfile, false);
  outfile.Printf("\n");
  // Start min 
  while (rms > min_tol_ && iteration < nMinSteps_) {
    // Calculate forces.
    if (potential.CalculateForce( frameIn )) {
      mprinterr("Error: Could not calculate force.\n");
      return 1;
    }
    double e_total = potential.Energy().Total();

    // Calculate the magnitude of the force vector.
    double sum = 0.0;
    const double* fxyz = frameIn.fAddress();
    for (int idx = 0; idx < frameIn.Natom(); idx++, fxyz += 3) {
      //mprintf("FDEBUG %8li%20.10f\n", fxyz-frameIn.fAddress()+1, fxyz[0]);
      //mprintf("FDEBUG %8li%20.10f\n", fxyz-frameIn.fAddress()+2, fxyz[1]);
      //mprintf("FDEBUG %8li%20.10f\n", fxyz-frameIn.fAddress()+3, fxyz[2]);
      sum += (fxyz[0]*fxyz[0] + fxyz[1]*fxyz[1] + fxyz[2]*fxyz[2]);
    }
    rms = sqrt( sum ) / fnq;
    //mprintf("DBG sum rms %20.10f%20.10f\n", sum, rms);
    // Adjust search step size
    if (dxst < crits) dxst = dxstm;
    dxst = dxst / 2.0;
    if (e_total < last_e) dxst = dxst * 2.4;
    double dxsth = dxst / sqrt( sum );
    last_e = e_total;
    // Update positions
    double* Xptr = frameIn.xAddress();
    fxyz = frameIn.fAddress();
    for (int idx = 0; idx != frameIn.Natom(); idx++, Xptr += 3, fxyz += 3)
    {
      //mprintf("xyz0= %g %g %g  Fxyz= %g %g %g\n", Xptr[0], Xptr[1], Xptr[2], fxyz[0], fxyz[1], fxyz[2]);
      Xptr[0] += fxyz[0] * dxsth;
      Xptr[1] += fxyz[1] * dxsth;
      Xptr[2] += fxyz[2] * dxsth;
      //fxyz[0] = 0.0;
      //fxyz[1] = 0.0;
      //fxyz[2] = 0.0;
      //mprintf("xyz1= %g %g %g\n", Xptr[0], Xptr[1], Xptr[2]);
      //*XV += (*FV * dxsth);
      //*FV = 0.0;
    }
    // Write out current E.
    outfile.Printf("%-8i %12.4E %12.4E", iteration+1, e_total, rms);
    // Do not print terms when only a single term.
    potential.Energy().PrintActiveTerms(outfile, false);
    outfile.Printf("\n");
    //outfile.Printf("%-8i %12.4E %12.4E EB=%12.4E EV=%12.4E EC=%12.4E\n",
    //               iteration+1, e_total, rms,
    //               potential.Energy().Ene(EnergyArray::E_BOND),
    //               potential.Energy().Ene(EnergyArray::E_VDW),
    //               potential.Energy().Ene(EnergyArray::E_COULOMB));
    //mprintf("Iteration:\t%8i %12.4E %12.4E EB=%12.4E EV=%12.4E EC=%12.4E\n",
    //        iteration, e_total, rms, E_bond, E_vdw, E_elec);
    iteration++;
    // Write out current coords
    if (trajOut.IsInitialized()) {
      if (trajOut.WriteSingle(iteration, frameIn)) return 1;
    }
  } // END minimization loop

  return 0;
}
