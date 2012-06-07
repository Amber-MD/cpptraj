#include <cmath> // sqrt
#include "Action_AtomicFluct.h"
#include "CpptrajStdio.h"
#include "Constants.h" // PI

// CONSTRUCTOR
AtomicFluct::AtomicFluct() {
  sets_ = 0;
  start_ = 0;
  stop_ = -1;
  offset_ = 1;
  bfactor_ = false;
  targetSet_ = 0;
  outfilename_ = NULL;
  fluctParm_ = NULL;
  outtype_ = BYATOM;
}

// AtomicFluct::init()
/** Usage: atomicfluct [out <filename>] [<mask>] [byres | byatom | bymask] [bfactor]
  *                    [start <start>] [stop <stop>] [offset <offset>]
  */               
int AtomicFluct::init() {
  // Get frame # keywords
  int startArg = actionArgs.getKeyInt("start",1);
  stop_ = actionArgs.getKeyInt("stop",-1);
  if (stop_==-1) 
    stop_ = actionArgs.getKeyInt("end",-1);
  // Cpptraj frame #s start from 0
  if (startArg<1) {
    mprinterr("Error: AtomicFluct: start arg must be >= 1 (%i)\n",startArg);
    return 1;
  }
  if (stop_!=-1 && startArg>stop_) {
    mprinterr("Error: AtomicFluct: start arg (%i) > stop arg (%i)\n",startArg,stop_);
    return 1;
  }
  start_ = startArg - 1;
  targetSet_ = start_;
  offset_ = actionArgs.getKeyInt("offset",1);
  if (offset_<1) {
    mprinterr("Error: AtomicFluct: offset arg must be >= 1 (%i)\n",offset_);
    return 1;
  }
  // Get other keywords
  bfactor_ = actionArgs.hasKey("bfactor");
  outfilename_ = actionArgs.getKeyString("out",NULL);
  if (actionArgs.hasKey("byres"))
    outtype_ = BYRES;
  else if (actionArgs.hasKey("bymask"))
    outtype_ = BYMASK;
  else if (actionArgs.hasKey("byatom") || actionArgs.hasKey("byatm"))
    outtype_ = BYATOM;
  // Get Mask
  char *maskstring = actionArgs.getNextMask();
  Mask.SetMaskString(maskstring);

  mprintf("    ATOMICFLUCT: calculating");
  if (bfactor_)
    mprintf(" B factors");
  else
    mprintf(" atomic positional fluctuations");
  if (outfilename_==NULL)
    mprintf(", output to STDOUT");
  else
    mprintf(", output to file %s",outfilename_);
  mprintf("\n                 Atom mask: [%s]\n",Mask.MaskString());
  if (start_>0 || stop_!=-1 || offset_!=1) {
    mprintf("                 Processing frames %i to",start_+1);
    if (stop_!=-1)
      mprintf(" %i",stop_);
    else
      mprintf(" end");
    if (offset_!=1)
      mprintf(", offset %i\n",offset_);
  }

  return 0;
}

// AtomicFluct::setup()
int AtomicFluct::setup() {

  if (SumCoords_.Natom()==0) {
    // Set up frames if not already set
    SumCoords_.SetupFrame(currentParm->Natom());
    SumCoords2_.SetupFrame(currentParm->Natom());
    SumCoords_.ZeroCoords();
    SumCoords2_.ZeroCoords();
    // This is the parm that will be used for this calc
    fluctParm_ = currentParm;
    // Set up atom mask
    if (currentParm->SetupCharMask( Mask )) {
      mprinterr("Error: AtomicFluct: Could not set up mask [%s]\n",Mask.MaskString());
      return 1;
    }
    if (Mask.None()) {
      mprinterr("Error: AtomicFluct: No atoms selected [%s]\n",Mask.MaskString());
      return 1;
    }
  } else if (currentParm->Natom() != SumCoords_.Natom()) {
    // Check that current #atoms matches
    mprinterr("Error: AtomicFluct not yet supported for mulitple topologies with different\n");
    mprinterr("       #s of atoms (set up for %i, this topology has %i\n",
              SumCoords_.Natom(), currentParm->Natom());
    return 1;
  } 
  // NOTE: Print warning here when setting up multiple topologies?
  return 0;
}

// AtomicFluct::action()
int AtomicFluct::action() {
  if (frameNum == targetSet_) {
    SumCoords_ += *currentFrame;
    SumCoords2_ += ( (*currentFrame) * (*currentFrame) ) ;
    ++sets_;
    targetSet_ += offset_;
    if (targetSet_ >= stop_ && stop_!=-1)
      targetSet_ = -1;
  }
  return 0;
}

// AtomicFluct::print() 
void AtomicFluct::print() {
  CpptrajFile outfile;
  int atom, res, atomnum, resstart, resstop;
  double xi, fluct, mass;

  mprintf("    ATOMICFLUCT: Calculating fluctuations for %i sets.\n",sets_);

  if (outfile.SetupWrite(outfilename_, debug)!=0) {
    mprinterr("Error: AtomicFluct: Could not set up output file %s\n",outfilename_);
    return;
  }
  if (outfile.OpenFile()) return;

  double Nsets = (double)sets_;
  SumCoords_.Divide(Nsets);
  SumCoords2_.Divide(Nsets);
  //SumCoords2_ = SumCoords2_ - (SumCoords_ * SumCoords);
  SumCoords_ *= SumCoords_;
  SumCoords2_ -= SumCoords_;

  std::vector<double> XYZ = SumCoords2_.ConvertToDouble();
  // DEBUG
  //mprintf("DEBUG: Converting to double: Original first coord:\n");
  //SumCoords2_.printAtomCoord(0);
  //mprintf("       Converted first coord: %lf %lf %lf\n",XYZ[0],XYZ[1],XYZ[2]);

  // Hold fluctuation results - initialize to 0
  std::vector<double> Results( SumCoords2_.Natom(), 0 );
  std::vector<double>::iterator result = Results.begin();

  if (bfactor_) {
    // Set up b factor normalization
    // B-factors are (8/3)*PI*PI * <r>**2 hence we do not sqrt the fluctuations
    double bfac = (8.0/3.0)*PI*PI;
    for (unsigned int i = 0; i < XYZ.size(); i+=3) {
      double fluct = XYZ[i] + XYZ[i+1] + XYZ[i+2];
      if (fluct > 0) 
        *result = bfac * fluct;
      ++result;
    }
  } else {
    // Atomic fluctuations
    for (unsigned int i = 0; i < XYZ.size(); i+=3) {
      double fluct = XYZ[i] + XYZ[i+1] + XYZ[i+2];
      if (fluct > 0)
        *result = sqrt(fluct);
      ++result;
    }
  }

  switch (outtype_) {
    // By atom output
    case BYATOM:
      atomnum = 1;
      for (atom = 0; atom < (int)Results.size(); atom++ ) {
        if (Mask.AtomInCharMask(atom)) 
          outfile.Printf(" %6i  %lf\n",atomnum++,Results[atom]);
      }
      break;
    // By residue output
    case BYRES:
      for (res = 0; res < fluctParm_->Nres(); res++) {
        xi = 0;
        fluct = 0;
        fluctParm_->ResAtomRange(res, &resstart, &resstop);
        for (atom = resstart; atom < resstop; atom++) {
          if ( Mask.AtomInCharMask(atom) ) {
            mass = (*fluctParm_)[atom].Mass(); 
            xi += mass;
            fluct += Results[atom] * mass;
          }
        }
        if (xi > SMALL)
          outfile.Printf(" %6i  %lf\n",res+1,fluct/xi);
      }
      break;
    // By mask output
    case BYMASK:
      xi = 0;
      fluct = 0;
      for (atom = 0; atom < (int)Results.size(); atom++) {
        if (Mask.AtomInCharMask(atom)) {
          mass = (*fluctParm_)[atom].Mass();
          xi += mass;
          fluct += Results[atom] * mass;
        }
      }
      if (xi > SMALL)
        outfile.Printf(" %6i  %lf\n",1, fluct/xi);
      break;
  }
  
  outfile.CloseFile();        
}

