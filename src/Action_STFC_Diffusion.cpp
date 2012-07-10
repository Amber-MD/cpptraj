#include "Action_STFC_Diffusion.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_STFC_Diffusion::Action_STFC_Diffusion() :
  printDistances_(false),
  calcType_(DEFAULT),
  direction_(DX),
  //outputAverDist_(NULL),
  //outputNumWat_(NULL),
  time_(1.0),
  lowerCutoff_(0),
  upperCutoff_(0),
  hasBox_(false),
  Ninitial_atm_(0),
  Ninitial_crd_(0),
  initialxyz_(0),
  distancexyz_(0),
  deltaxyz_(0),
  previousxyz_(0),
  elapsedFrames_(0)
{}

Action_STFC_Diffusion::~Action_STFC_Diffusion() {
  if (initialxyz_ != 0) delete[] initialxyz_;
  if (distancexyz_ != 0) delete[] distancexyz_;
  if (deltaxyz_ != 0) delete[] deltaxyz_;
  if (previousxyz_ != 0) delete[] previousxyz_;
}

/*const char* Action_STFC_Diffusion::DirectionString[7][4] = {
  "x", "y", "z", "xy", "xz", "yz", "xyz"
};*/

// Action_STFC_Diffusion::init()
/** Usage: diffusion mask <mask> [out <file>] [time <time per frame>]
 *                   ([mask2 <mask>] [lower <distance>] [upper <distance>]
 *                   [nwout <file>]) [avout <file>] [distances] [com]
 *                   [x|y|z|xy|xz|yz|xyz]
 */
int Action_STFC_Diffusion::init() {
  char* maskarg = actionArgs.getKeyString("mask", NULL);
  if (maskarg == NULL) {
    mprinterr("Error: diffusion: No mask specified.\n");
    return 1;
  }
  mask_.SetMaskString( maskarg );

  std::string outfileName = actionArgs.GetStringKey("out");
  if (outfileName.empty())
    outfileName = "diffusion.dat";
  if ( output_.OpenWrite( outfileName ) ) {
    mprinterr("Error: diffusion: Could not open output file %s\n", outfileName.c_str());
    return 1;
  }

  outputAverDist_ = actionArgs.GetStringKey("avout");
  if (!outputAverDist_.empty()) {
    if ( outputad_.OpenWrite( outputAverDist_ ) ) {
      mprinterr("Error: diffusion: Could not open average distance file %s\n", 
                outputAverDist_.c_str()); 
      return 1;
    }
  }

  time_ = actionArgs.getKeyDouble("time", 1.0);
  if (time_ < 0) {
    mprinterr("Error: diffusion: time argument cannot be < 0 (%lf)\n", time_);
    return 1;
  }

  printDistances_ = actionArgs.hasKey("distances");
  if ( actionArgs.hasKey("com") )
    calcType_ = COM;

  // Directions considered for diffusion
  direction_ = DXYZ;
  if ( actionArgs.hasKey("x") )   direction_ = DX;
  if ( actionArgs.hasKey("y") )   direction_ = DY;
  if ( actionArgs.hasKey("z") )   direction_ = DZ;
  if ( actionArgs.hasKey("xy") )  direction_ = DXY;
  if ( actionArgs.hasKey("xz") )  direction_ = DXZ;
  if ( actionArgs.hasKey("yz") )  direction_ = DYZ;
  if ( actionArgs.hasKey("xyz") ) direction_ = DXYZ;
    
  // Process second mask
  maskarg = actionArgs.getKeyString("mask2", NULL);
  double lcut = 0;
  double ucut = 0;
  if (maskarg != NULL) {
    mask2_.SetMaskString( maskarg );
    lcut = actionArgs.getKeyDouble("lower", 0.01);
    ucut = actionArgs.getKeyDouble("upper", 3.5 );
    lowerCutoff_ = lcut * lcut;
    upperCutoff_ = ucut * ucut;
    outputNumWat_ = actionArgs.GetStringKey("nwout");
    if (outputNumWat_.empty())
      outputNumWat_ = "nw.dat";
    if ( outputnw_.OpenWrite( outputNumWat_ ) ) {
      mprinterr("Error: diffusion: Could not open nwout file %s\n", 
                outputNumWat_.c_str());
      return 1;
    }
    calcType_ = DIST;
  }

  if (calcType_ != DEFAULT)
    printDistances_ = false;

  // Write info.
  mprintf("    DIFFUSION (STFC): Calculating diffusion in the");
  switch (direction_) {
    case DX:   mprintf(" x direction"); break;
    case DY:   mprintf(" y direction"); break;
    case DZ:   mprintf(" z direction"); break;
    case DXY:  mprintf(" xy plane"); break;
    case DXZ:  mprintf(" xz plane"); break;
    case DYZ:  mprintf(" yz plane"); break;
    case DXYZ: mprintf(" xyz directions"); break;
  }
  mprintf("\n\t\tMask 1 expression: %s\n", mask_.MaskString());
  if (calcType_ == COM)
    mprintf("\t\tCenter of mass diffusion of atoms in mask1 will be computed\n");
  else if (calcType_ == DIST)
    mprintf("\t\tAtoms in mask 2 (%s) in the range %.3f to %.3f Angstrom will be used\n",
            mask2_.MaskString(), lcut, ucut);

  if (!printDistances_)
    mprintf("\t\tOnly the average");
  else
    mprintf("\t\tThe average and individual");
  mprintf(" results will be written to %s\n", output_.Name());

  if (calcType_ == DIST)
    mprintf("\t\tThe number of atoms in the shell will be written to %s\n",
            outputnw_.Name());

  if (!outputAverDist_.empty())
    mprintf("\t\t<dr^2> will be written to %s\n", outputad_.Name());

  mprintf("\t\tThe time step between frames is %.3f ps.\n", time_);

  return 0;
}

// Action_STFC_Diffusion::setup()
int Action_STFC_Diffusion::setup() {
  // Setup atom mask
  if (currentParm->SetupIntegerMask( mask_ )) return 1;
  if (mask_.None()) {
    mprinterr("Error: diffusion: No atoms selected.\n");
    return 1;
  }
  // Setup second mask if necessary
  if ( calcType_ == DIST ) {
    if (currentParm->SetupIntegerMask( mask2_ )) return 1;
    if (mask2_.None()) {
      mprinterr("Error: diffusion: No atoms selected in second mask.\n");
      return 1;
    }
  }

  // Check for box
  if ( currentParm->BoxType()!=Box::NOBOX )
    hasBox_ = true;
  else
    hasBox_ = false;

  // If initial frame already set and current # atoms > # atoms in initial
  // frame this will probably cause an error.
  // NOTE: Shouldnt matter for COM.
  //if (!initialxyz_.empty() && currentParm->Natom() > (int)initialxyz_.size()) {
  if ( calcType_ != COM ) {
    if ( initialxyz_ != 0 && currentParm->Natom() > Ninitial_atm_ ) {
      mprintf("Warning: # atoms in current parm (%s, %i) > # atoms in initial frame (%i)\n",
               currentParm->c_str(), currentParm->Natom(), Ninitial_atm_);
      mprintf("Warning: This may lead to segmentation faults.\n");
    }
  }

  // Based on the calculation type:
  //   1- Allocate the distance arrays
  //   2- Allocate the delta arrays
  //   3- Reserve space for the initial and previous frame arrays
  // NOTE: Currently only performing allocation on first setup. This means
  //       the code is currently not safe for changing # of atoms.
  if (initialxyz_ == 0) {
    if ( calcType_ == DEFAULT ) { // All
        Ninitial_atm_ = currentParm->Natom();
        Ninitial_crd_ = Ninitial_atm_ * 3;
        initialxyz_ = new double[ Ninitial_crd_ ];
        int selectedcrd = mask_.Nselected() * 3;
        previousxyz_ = new double[ selectedcrd ];
        distancexyz_ = new double[ selectedcrd ];
        distance_.resize( mask_.Nselected() );
        deltaxyz_ = new double[ selectedcrd ]; 
/*
    initialx_.reserve( currentParm->Natom() );
    initialy_.reserve( currentParm->Natom() );
    initialz_.reserve( currentParm->Natom() );
    previousx_.reserve( mask_.Nselected() );
    previousy_.reserve( mask_.Nselected() );
    previousz_.reserve( mask_.Nselected() );
    distancex_.resize( mask_.Nselected() );
    distancey_.resize( mask_.Nselected() );
    distancez_.resize( mask_.Nselected() );
    distance_.resize(  mask_.Nselected() );
    deltax_.assign( mask_.Nselected(), 0 );
    deltay_.assign( mask_.Nselected(), 0 );
    deltaz_.assign( mask_.Nselected(), 0 );
*/
    } else if ( calcType_ == COM ) { // Center of mass
      Ninitial_atm_ = 1;
      Ninitial_crd_ = 3;
      initialxyz_ = new double[ 3 ];
      previousxyz_ = new double[ 3 ];
      distancexyz_ = new double[ 3 ];
      distance_.resize( 1 );
      deltaxyz_ = new double[ 3 ];
/*
    initialx_.resize( 1 );
    initialy_.resize( 1 );
    initialz_.resize( 1 );
    previousx_.resize( 1 );
    previousy_.resize( 1 );
    previousz_.resize( 1 );
    distancex_.resize( 1 );
    distancey_.resize( 1 );
    distancez_.resize( 1 );
    distance_.resize( 1 );
    deltax_.assign( 1, 0 );
    deltay_.assign( 1, 0 );
    deltaz_.assign( 1, 0 );
*/
    } else if ( calcType_ == DIST ) { // Region based
      Ninitial_atm_ = currentParm->Natom();
      Ninitial_crd_ = Ninitial_atm_ * 3;
      initialxyz_ = new double[ Ninitial_crd_ ];
      previousxyz_ = new double[ Ninitial_crd_ ];
      distancexyz_ = new double[ Ninitial_crd_ ];
      distance_.resize( Ninitial_atm_ );
      deltaxyz_ = new double[ Ninitial_crd_ ];
/*
    initialx_.reserve( currentParm->Natom() );
    initialy_.reserve( currentParm->Natom() );
    initialz_.reserve( currentParm->Natom() );
    previousx_.reserve( currentParm->Natom() );
    previousy_.reserve( currentParm->Natom() );
    previousz_.reserve( currentParm->Natom() );
    distancex_.resize( currentParm->Natom() );
    distancey_.resize( currentParm->Natom() );
    distancez_.resize( currentParm->Natom() );
    distance_.resize(  currentParm->Natom() );
    deltax_.assign( currentParm->Natom(), 0 );
    deltay_.assign( currentParm->Natom(), 0 );
    deltaz_.assign( currentParm->Natom(), 0 );
*/
    }
  }

  // Allocate the dSum arrays
  // NOTE: Should this be resize?
  dSum1_.assign( currentParm->Natom(), 0);
  dSum2_.assign( currentParm->Natom(), 0);

  return 0;
}

void Action_STFC_Diffusion::calculateMSD(double* XYZ1, double* XYZ2) {

}

// Action_STFC_Diffusion::action()
int Action_STFC_Diffusion::action() {
  double XYZ[3];
  // Load initial frame if necessary
  if (elapsedFrames_ == 0) {
    if ( calcType_ == DEFAULT ) { // All
      // Save all initial coords, selected previous coords.
      AtomMask::const_iterator selected_atom = mask_.begin();
      double* initial = initialxyz_;
      double* previous = previousxyz_;
      for (int atnum = 0; atnum < Ninitial_atm_; ++atnum) {
        currentFrame->GetAtomXYZ(initial, atnum);
        if ( atnum == *selected_atom ) {
          previous[0] = initial[0];
          previous[1] = initial[1];
          previous[2] = initial[2];
          ++selected_atom;
        }
        initial += 3;
        previous += 3;
      }
    } else if ( calcType_ == COM ) { // Center of Mass
      // Save initial COM
      currentFrame->CenterOfMass( &mask_, initialxyz_ );
      previousxyz_[0] = initialxyz_[0];
      previousxyz_[1] = initialxyz_[1];
      previousxyz_[2] = initialxyz_[2];
    } else if (calcType_ == DIST ) { // Region Based
      // Save all coords 
      double* initial = initialxyz_;
      double* previous = previousxyz_;
      for (int atnum = 0; atnum < currentParm->Natom(); ++atnum) {
        currentFrame->GetAtomXYZ(initial, atnum);
        previous[0] = initial[0];
        previous[1] = initial[1];
        previous[2] = initial[2];
        initial += 3;
        previous += 3;
      }
    }

  // Initial frame is loaded
  } else {
    if ( calcType_ == DEFAULT ) {
      double* previous = previousxyz_;
      for ( AtomMask::const_iterator atom = mask_.begin(); atom != mask_.end(); ++atom) 
      {
        currentFrame->GetAtomXYZ(XYZ, *atom);
        calculateMSD( XYZ, previous );
        previous += 3; 
      }
    } else if (calcType_ == COM) {
      currentFrame->CenterOfMass( &mask_, XYZ );
      calculateMSD( XYZ, previousxyz_ );
    } else if (calcType_ == DIST) {
      // Check which atoms from mask1 are inside or outside the region defined
      // by lower and upper.
    }
  }

  return 0;
}
  
