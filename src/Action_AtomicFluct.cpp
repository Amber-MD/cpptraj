#include <cmath> // sqrt
#include "Action_AtomicFluct.h"
#include "CpptrajStdio.h"
#include "Constants.h" // PI

// CONSTRUCTOR
Action_AtomicFluct::Action_AtomicFluct() {
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

// Action_AtomicFluct::init()
/** Usage: atomicfluct [out <filename>] [<mask>] [byres | byatom | bymask] [bfactor]
  *                    [start <start>] [stop <stop>] [offset <offset>]
  */               
int Action_AtomicFluct::init() {
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
  // Get DataSet name
  setname_ = actionArgs.GetStringNext();

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
  if (!setname_.empty())
    mprintf("\tData will be saved to set named %s\n", setname_.c_str());

  return 0;
}

// Action_AtomicFluct::setup()
int Action_AtomicFluct::setup() {

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

// Action_AtomicFluct::action()
int Action_AtomicFluct::action() {
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

// Action_AtomicFluct::print() 
void Action_AtomicFluct::print() {
  int atom, res, resstart, resstop; 
  double xi, fluct, mass;

  mprintf("    ATOMICFLUCT: Calculating fluctuations for %i sets.\n",sets_);
  DataSet* dataout = DSL->AddSet( DataSet::DOUBLE, setname_, "Fluct" );
  if (dataout == NULL) {
    mprinterr("Error: AtomicFluct: Could not allocate dataset for output.\n");
    return; 
  }
  DataFile* outfile = DFL->Add( outfilename_, dataout );
  outfile->ProcessArgs("noemptyframes");

  double Nsets = (double)sets_;
  SumCoords_.Divide(Nsets);
  SumCoords2_.Divide(Nsets);
  //SumCoords2_ = SumCoords2_ - (SumCoords_ * SumCoords);
  SumCoords_ *= SumCoords_;
  SumCoords2_ -= SumCoords_;

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
    outfile->ProcessArgs( "ylabel B-factors" );
    double bfac = (8.0/3.0)*PI*PI;
    for (int i = 0; i < SumCoords2_.size(); i+=3) {
      double fluct = SumCoords2_[i] + SumCoords2_[i+1] + SumCoords2_[i+2];
      if (fluct > 0) 
        *result = bfac * fluct;
      ++result;
    }
  } else {
    // Atomic fluctuations
    for (int i = 0; i < SumCoords2_.size(); i+=3) {
      double fluct = SumCoords2_[i] + SumCoords2_[i+1] + SumCoords2_[i+2];
      if (fluct > 0)
        *result = sqrt(fluct);
      ++result;
    }
  }

  switch (outtype_) {
    // By atom output
    case BYATOM:
      outfile->ProcessArgs("xlabel Atom");
      for (atom = 0; atom < (int)Results.size(); atom++ ) {
        if (Mask.AtomInCharMask(atom)) 
          dataout->Add( atom, &(Results[atom]) );
      }
      break;
    // By residue output
    case BYRES:
      outfile->ProcessArgs("xlabel Res");
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
        if (xi > SMALL) {
          mass = fluct / xi;
          dataout->Add( res, &mass );
        }
      }
      break;
    // By mask output
    case BYMASK:
      outfile->ProcessArgs( "xlabel " + Mask.MaskExpression() );
      xi = 0;
      fluct = 0;
      for (atom = 0; atom < (int)Results.size(); atom++) {
        if (Mask.AtomInCharMask(atom)) {
          mass = (*fluctParm_)[atom].Mass();
          xi += mass;
          fluct += Results[atom] * mass;
        }
      }
      if (xi > SMALL) {
        mass = fluct / xi;
        dataout->Add( 0, &mass );
      }
      break;
  }
}

