#include "Exec_Emin.h"
#include "CpptrajStdio.h"
#include "PotentialFunction.h"
#include "Minimize_SteepestDescent.h"
#include "MdOpts.h"

// Exec_Emin::Help()
void Exec_Emin::Help() const
{
  mprintf("\tcrdset <name> [trajoutname <name>] [rmstol <tol>] [nsteps <#>]\n"
          "\t[<mask>] [frame <#>] [dx0 <step0>] [out <file>]\n");
}

// Exec_Emin::Execute()
Exec::RetType Exec_Emin::Execute(CpptrajState& State, ArgList& argIn)
{
  mprintf("Warning: THIS COMMAND IS STILL UNDER DEVELOPMENT.\n");

  std::string setname = argIn.GetStringKey("crdset");
  if (setname.empty()) {
    mprinterr("Error: Specify COORDS set to minimize with 'crdset'\n");
    return CpptrajState::ERR;
  }
  DataSet_Coords* crdset = (DataSet_Coords*)State.DSL().FindSetOfGroup(setname, DataSet::COORDINATES);
  if (crdset == 0) {
    mprinterr("Error: No COORDS set found with name '%s'\n", setname.c_str());
    return CpptrajState::ERR;
  }

  // Get the frame # to be minimized.
  int framenum = argIn.getKeyInt("frame", 0);
  mprintf("\tMinimizing COORDS set '%s' frame %i\n", crdset->legend(), framenum+1);

  if (framenum < 0 || framenum >= (int)crdset->Size()) {
    mprinterr("Error: Frame %i is out of range.\n", framenum+1);
    return CpptrajState::ERR;
  }

  // Get the frame. Instead of using AllocateFrame, allocate manually because
  // we need to ensure space for forces.
  //Frame frameIn = crdset->AllocateFrame();
  Frame frameIn;
  CoordinateInfo cinfo = crdset->CoordsInfo();
  cinfo.SetForce( true );
  frameIn.SetupFrameV(crdset->Top().Atoms(), cinfo);
  crdset->GetFrame(framenum, frameIn);

  std::string trajoutname = argIn.GetStringKey("trajoutname");
  if (!trajoutname.empty())
    mprintf("\tOutput trajectory: %s\n", trajoutname.c_str());

  CpptrajFile* outfile = State.DFL().AddCpptrajFile(argIn.GetStringKey("out"),
                                                    "Min. Out", DataFileList::TEXT, true);
  if (outfile == 0) {
    mprinterr("Internal Error: Could not allocate output file for minimization.\n");
    return CpptrajState::ERR;
  }
  mprintf("\tOutput to %s\n", outfile->Filename().full());

  double min_tol = argIn.getKeyDouble("rmstol", 1E-4);
  mprintf("\tMin RMS tolerance: %g\n", min_tol);

  double dx0 = argIn.getKeyDouble("dx0", 0.01);
  mprintf("\tInitial step size: %g\n", dx0);

  int nMinSteps = argIn.getKeyInt("nsteps", 1);
  mprintf("\t%i minimization steps.\n", nMinSteps);

  std::string maskexpr = argIn.GetMaskNext();
  if (!maskexpr.empty())
    mprintf("\tMask expression: %s\n", maskexpr.c_str());

  bool use_openmm = argIn.hasKey("openmm");
  // Create the potential function. This is done last so potential term
  // arguments are parsed last.
  PotentialFunction potential;
  MdOpts opts;
  if (opts.GetOptsFromArgs(argIn)) return CpptrajState::ERR;
  opts.PrintOpts();
  if (use_openmm)
    potential.AddTerm( PotentialTerm::OPENMM, opts );
  else
    potential.AddTerm( PotentialTerm::BOND, opts );
  Minimize_SteepestDescent SD;

  // Set up the potential function
  if (potential.SetupPotential( crdset->Top(), frameIn.BoxCrd(), maskexpr )) {
    mprinterr("Error: Could not set up potential.\n");
    return CpptrajState::ERR;
  }

  // Set up and run minimization
  if (SD.SetupMin(trajoutname, min_tol, dx0, nMinSteps)) {
    mprinterr("Error: Could not set up minimizer.\n");
    return CpptrajState::ERR;
  }

  if (SD.RunMin(potential, frameIn, *outfile)) {
    mprinterr("Error: Minimization failed.\n");
    return CpptrajState::ERR;
  }

  return CpptrajState::OK;
}
