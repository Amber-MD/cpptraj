#include <cmath>
#include "Action_STFC_Diffusion.h"
#include "CpptrajStdio.h"
#include "DistRoutines.h"

// CONSTRUCTOR
Action_STFC_Diffusion::Action_STFC_Diffusion() :
  printDistances_(false),
  calcType_(DEFAULT),
  direction_(DX),
  output_(0), outputnw_(0), outputad_(0),
  time_(1.0),
  lowerCutoff_(0),
  upperCutoff_(0),
  hasBox_(false),
  n_atom_(-1),
  elapsedFrames_(0)
{}

void Action_STFC_Diffusion::Help() const {
  mprintf("\tmask <mask> [out <file>] [time <time per frame>]\n"
          "\t[mask2 <mask> [lower <distance>] [upper <distance>]]\n"
          "\t[nwout <file>]) [avout <file>] [distances] [com]\n"
          "\t[x|y|z|xy|xz|yz|xyz]\n"
          "  Calculate diffusion of atoms in <mask>\n");
}

/// Must correspond to DirectionType
static const char* DirectionString[] = {
  "x",  "y",  "z",  "xy", "xz", "yz", "xyz" 
};

// Action_STFC_Diffusion::Init()
Action::RetType Action_STFC_Diffusion::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
# ifdef MPI
  if (init.TrajComm().Size() > 1) {
    mprinterr("Error: 'stfcdiffusion' action does not work with > 1 process"
              " (%i processes currently).\n", init.TrajComm().Size());
    return Action::ERR;
  }
# endif
  n_atom_ = -1;
  // Get keywords
  std::string maskarg = actionArgs.GetStringKey("mask");
  if (maskarg.empty()) {
    mprinterr("Error: No mask specified.\n");
    return Action::ERR;
  }
  mask_.SetMaskString( maskarg );

  std::string outfileName = actionArgs.GetStringKey("out");
  if (outfileName.empty())
    outfileName = "diffusion.dat";
  output_ = init.DFL().AddCpptrajFile(outfileName, "Diffusion");
  if ( output_ == 0 ) {
    mprinterr("Error: Could not open output file '%s'\n", outfileName.c_str());
    return Action::ERR;
  }

  outputad_ = init.DFL().AddCpptrajFile(actionArgs.GetStringKey("avout"), "Diffusion Avg Dist");
  time_ = actionArgs.getKeyDouble("time", 1.0);
  if (time_ < 0) {
    mprinterr("Error: time argument cannot be < 0 (%f)\n", time_);
    return Action::ERR;
  }
  printDistances_ = actionArgs.hasKey("distances");
  if ( actionArgs.hasKey("com") )
    calcType_ = COM;

  // Directions considered for diffusion
  direction_ = DXYZ;
  for (int i = 0; i <= (int)DXYZ; i++)
    if (actionArgs.hasKey( DirectionString[i] ))
      direction_ = (DirectionType)i;
  // Process second mask
  maskarg = actionArgs.GetStringKey("mask2");
  double lcut = 0;
  double ucut = 0;
  if (!maskarg.empty()) {
    mask2_.SetMaskString( maskarg );
    lcut = actionArgs.getKeyDouble("lower", 0.01);
    ucut = actionArgs.getKeyDouble("upper", 3.5 );
    lowerCutoff_ = lcut * lcut;
    upperCutoff_ = ucut * ucut;
    std::string outputNumWat = actionArgs.GetStringKey("nwout");
    if (outputNumWat.empty())
      outputNumWat = "nw.dat";
    outputnw_ = init.DFL().AddCpptrajFile(outputNumWat, "Diffusion # waters");
    if ( outputnw_ == 0 ) {
      mprinterr("Error: Could not open diffusion number of waters output file '%s'\n", 
                outputNumWat.c_str());
      return Action::ERR;
    }
    calcType_ = DIST;
    // See if imaging is to be performed.
    image_.InitImaging( !(actionArgs.hasKey("noimage")) );
  } else if (actionArgs.Contains("lower")) {
    mprinterr("Error: 'lower' requires 'mask2'\n");
    return Action::ERR;
  } else if (actionArgs.Contains("upper")) {
    mprinterr("Error: 'upper' requires 'mask2'\n");
    return Action::ERR;
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
  else if (calcType_ == DIST) {
    mprintf("\t\tAtoms in mask 2 (%s) in the range %.3f to %.3f Angstrom will be used\n",
            mask2_.MaskString(), lcut, ucut);
    if (image_.UseImage())
      mprintf("\t\tDistances will be imaged.\n");
    else
      mprintf("\t\tDistances will not be imaged.\n");
  }
  if (!printDistances_)
    mprintf("\t\tOnly the average");
  else
    mprintf("\t\tThe average and individual");
  mprintf(" results will be written to %s\n", output_->Filename().full());

  if (calcType_ == DIST)
    mprintf("\t\tThe number of atoms in the shell will be written to %s\n",
            outputnw_->Filename().full());

  if (outputad_ != 0)
    mprintf("\t\t<dr^2> will be written to %s\n", outputad_->Filename().full());

  mprintf("\t\tThe time step between frames is %.3f ps.\n", time_);

  return Action::OK;
}

// Action_STFC_Diffusion::Setup()
Action::RetType Action_STFC_Diffusion::Setup(ActionSetup& setup) {
  // Setup atom mask
  if (setup.Top().SetupIntegerMask( mask_ )) return Action::ERR;
  mask_.MaskInfo();
  if (mask_.None()) {
    mprintf("Warning: No atoms selected.\n");
    return Action::SKIP;
  }
  if (n_atom_ == -1) { // first time through, write header
    output_->Printf("%-10s %10s %10s %10s %10s","#time","x","y","z",DirectionString[direction_]);
    if (printDistances_) {
      for (AtomMask::const_iterator atom = mask_.begin(); atom != mask_.end(); ++atom) {
        int a1 = *atom + 1;
        output_->Printf(" x%-8i y%-8i z%-8i r%-8i", a1, a1, a1, a1);
      }
    }
    output_->Printf("\n");
  }
  n_atom_ = setup.Top().Natom();
  // Setup second mask if necessary
  if ( calcType_ == DIST ) {
    if (setup.Top().SetupIntegerMask( mask2_ )) return Action::ERR;
    mask2_.MaskInfo();
    if (mask2_.None()) {
      mprinterr("Error: No atoms selected by second mask.\n");
      return Action::ERR;
    }
    // Set up imaging info
    image_.SetupImaging( setup.CoordInfo().TrajBox().Type() );
    if (image_.ImagingEnabled())
      mprintf("\tImaging distances.\n");
    else
      mprintf("\tImaging off.\n");
  }

  // Check for box
  if ( setup.CoordInfo().TrajBox().Type()!=Box::NOBOX )
    hasBox_ = true;
  else
    hasBox_ = false;

  // If initial frame already set and current # atoms > # atoms in initial
  // frame this will probably cause an error.
  // NOTE: Shouldnt matter for COM.
  if ( calcType_ != COM ) {
    int initial_natom = (int)initialxyz_.size() / 3;
    if ( !initialxyz_.empty() && n_atom_ > initial_natom ) {
      mprintf("Warning: # atoms in current parm (%s, %i) > # atoms in initial frame (%i)\n",
               setup.Top().c_str(), n_atom_, initial_natom );
      mprintf("Warning: This may lead to segmentation faults.\n");
    }
  }

  // Based on the calculation type:
  //   1- Reserve space for the initial and previous frame arrays
  //   2- Allocate the distance arrays
  //   3- Allocate the delta arrays
  if ( calcType_ == DEFAULT ) { // All
      initialxyz_.reserve( n_atom_*3 );
      int selectedcrd = mask_.Nselected() * 3;
      previousxyz_.reserve( selectedcrd );
      distancexyz_.resize( selectedcrd );
      distance_.resize( mask_.Nselected() );
      deltaxyz_.assign( selectedcrd, 0 ); 
  } else if ( calcType_ == COM ) { // Center of mass
    initialxyz_.reserve( 3 );
    previousxyz_.reserve( 3 );
    distancexyz_.resize( 3 );
    distance_.resize( 1 );
    deltaxyz_.resize( 3 );
  } else if ( calcType_ == DIST ) { // Region based
    int ncrd = n_atom_*3;
    initialxyz_.reserve( ncrd );
    previousxyz_.reserve( ncrd );
    distancexyz_.resize( ncrd );
    distance_.resize( n_atom_ );
    deltaxyz_.assign( ncrd, 0 );
    nInside_.resize( n_atom_ );
  }

  // Allocate the dSum arrays
  dSum1_.resize( n_atom_, 0);
  dSum2_.resize( n_atom_, 0);

  return Action::OK;
}

// Action_STFC_Diffusion::calculateMSD()
/** Calculate Mean Square displacement.
  * \author Originally by Hannes H. Loeffler.
  * \author Adapted by Daniel R. Roe.
  * \param XYZ Current coordinates.
  * \param idx1 Original atom index.
  * \param idx2 Current atom index.
  * \param box Current box coordinates.
  */
void Action_STFC_Diffusion::calculateMSD(const double* XYZ, int idx1, int idx2, Vec3 const& box)
{
  double sum = 0;
  double sum2 = 0;

  // Calculate distance
  int idx23 = idx2 * 3;
  double delx = XYZ[0] - previousxyz_[idx23  ];
  double dely = XYZ[1] - previousxyz_[idx23+1];
  double delz = XYZ[2] - previousxyz_[idx23+2];

  // If the particle moved more than half the box, assume it was imaged and adjust
  // the distance of the total movement with respect to the original frame. 
  if (box[0] > 0.0) {
    if (delx > box[0]/2.0)
      deltaxyz_[idx23  ] -= box[0];
    else if ( delx < -box[0]/2.0)
      deltaxyz_[idx23  ] += box[0];
    if (dely > box[1]/2.0)
      deltaxyz_[idx23+1] -= box[1];
    else if ( dely < -box[1]/2.0)
      deltaxyz_[idx23+1] += box[1];
    if (delz > box[2]/2.0)
      deltaxyz_[idx23+2] -= box[2];
    else if (delz < -box[2]/2.0)
      deltaxyz_[idx23+2] += box[2];
  }

  // Set the current x with reference to the un-imaged trajectory
  double xx = XYZ[0] + deltaxyz_[idx23  ];
  double yy = XYZ[1] + deltaxyz_[idx23+1];
  double zz = XYZ[2] + deltaxyz_[idx23+2];
  //mprintf("DEBUG: xx=%f yy=%f zz=%f\n", xx, yy, zz);

  // Calculate the distance between this "fixed" coordinate and the
  // reference (initial) frame
  int idx13 = idx1 * 3; 
  delx = xx - initialxyz_[idx13  ];
  dely = yy - initialxyz_[idx13+1];
  delz = zz - initialxyz_[idx13+2];
  //mprintf("DEBUG: delx=%f dely=%f delz=%f   ix=%f iy=%f iz=%f\n", delx, dely, delz, initialxyz_[idx13  ], initialxyz_[idx13+1], initialxyz_[idx13+2]);

  // store the distance for this atom
  distancexyz_[idx23  ] = delx*delx;
  distancexyz_[idx23+1] = dely*dely;
  distancexyz_[idx23+2] = delz*delz;

  switch (direction_) {
    case DX: sum = distancexyz_[idx23  ]; sum2 = xx * xx; break;
    case DY: sum = distancexyz_[idx23+1]; sum2 = yy * yy; break;
    case DZ: sum = distancexyz_[idx23+2]; sum2 = zz * zz; break;
    case DXY: 
      sum = distancexyz_[idx23  ] + distancexyz_[idx23+1];
      sum2 = (xx * xx) + (yy * yy);
      break;
    case DXZ: 
      sum = distancexyz_[idx23  ] + distancexyz_[idx23+2];
      sum2 = (xx * xx) + (zz * zz);
      break;
    case DYZ: 
      sum = distancexyz_[idx23+1] + distancexyz_[idx23+2];
      sum2 = (yy * yy) + (zz * zz);
      break;
    case DXYZ:
      sum = distancexyz_[idx23  ] + distancexyz_[idx23+1] + distancexyz_[idx23+2];
      sum2 = (xx * xx) + (yy * yy) + (zz * zz);
      break;
  }
  //mprintf("DEBUG: sum=%f sum2=%f\n", sum, sum2);

  distance_[idx2] = sum;
  dSum1_[idx2] += sum2;       // sum r^2
  dSum2_[idx2] += sqrt(sum2); // sum r

  // Update previous coords
  previousxyz_[idx23  ] = XYZ[0];
  previousxyz_[idx23+1] = XYZ[1];
  previousxyz_[idx23+2] = XYZ[2];
}

// Action_STFC_Diffusion::DoAction()
Action::RetType Action_STFC_Diffusion::DoAction(int frameNum, ActionFrame& frm) {
  double Time, average, avgx, avgy, avgz;

  // ----- Load initial frame if necessary -------
  if ( initialxyz_.empty() ) {
    //mprintf("DEBUG: Initial frame is empty, mode %i\n", (int)calcType_);
    if ( calcType_ == DEFAULT ) { // All
      // Save all initial coords TODO use copy
      for (int aidx = 0; aidx != frm.Frm().Natom(); ++aidx) {
        const double* XYZ = frm.Frm().XYZ( aidx );
        initialxyz_.push_back( XYZ[0] );
        initialxyz_.push_back( XYZ[1] );
        initialxyz_.push_back( XYZ[2] );
      }
      // Save selected previous coords.
      for (AtomMask::const_iterator selected_atom = mask_.begin();
                                    selected_atom != mask_.end(); ++selected_atom)
      {
        const double* XYZ = frm.Frm().XYZ(*selected_atom);
        previousxyz_.push_back( XYZ[0] );
        previousxyz_.push_back( XYZ[1] );
        previousxyz_.push_back( XYZ[2] );
      }
    } else if ( calcType_ == COM ) { // Center of Mass
      // Save initial COM
      Vec3 XYZ = frm.Frm().VCenterOfMass( mask_ );
      initialxyz_.push_back( XYZ[0] );
      previousxyz_.push_back( XYZ[0] );
      initialxyz_.push_back( XYZ[1] );
      previousxyz_.push_back( XYZ[1] );
      initialxyz_.push_back( XYZ[2] );
      previousxyz_.push_back( XYZ[2] );
    } else if (calcType_ == DIST ) { // Region Based
      // Save all coords 
      for (int atnum = 0; atnum < n_atom_; ++atnum) {
        const double* XYZ = frm.Frm().XYZ(atnum);
        initialxyz_.push_back( XYZ[0] );
        previousxyz_.push_back( XYZ[0] );
        initialxyz_.push_back( XYZ[1] );
        previousxyz_.push_back( XYZ[1] );
        initialxyz_.push_back( XYZ[2] );
        previousxyz_.push_back( XYZ[2] );
      }
    }
    return Action::OK;
  } 

  // ----- Initial frame is loaded ---------------
  ++elapsedFrames_;
  Time = elapsedFrames_ * time_;
  // For accumulating averages
  average = 0;
  avgx = 0;
  avgy = 0;
  avgz = 0;

  Vec3 Box = frm.Frm().BoxCrd().Lengths();
  if ( calcType_ == DEFAULT ) 
  {
    int idx2 = 0;
    int idx23 = 0;
    for ( AtomMask::const_iterator atom = mask_.begin(); atom != mask_.end(); ++atom)
    {
      calculateMSD( frm.Frm().XYZ(*atom), *atom, idx2, Box );
      average += distance_[idx2];
      avgx += distancexyz_[idx23++];
      avgy += distancexyz_[idx23++];
      avgz += distancexyz_[idx23++];
      ++idx2;
    }
    average /= (double)mask_.Nselected();
    avgx /= (double)mask_.Nselected();
    avgy /= (double)mask_.Nselected();
    avgz /= (double)mask_.Nselected();
  } 
  else if (calcType_ == COM) 
  {
    Vec3 XYZ = frm.Frm().VCenterOfMass( mask_ );
    //mprintf("CDBG:\tXYZ[%i] = %f %f %f\n", elapsedFrames_,XYZ[0], XYZ[1], XYZ[2]);
    calculateMSD( XYZ.Dptr(), 0, 0, Box );
    average = distance_[0];
    avgx = distancexyz_[0];
    avgy = distancexyz_[1];
    avgz = distancexyz_[2];
  } 
  else if (calcType_ == DIST) 
  {
    // Check which atoms from mask1 are inside or outside the region defined
    // by lower and upper.
    //mprintf("DEBUG: distances mode frame %i\n", elapsedFrames_);
    nInside_.assign( n_atom_, 0 );
    for ( AtomMask::const_iterator atom1 = mask_.begin(); atom1 != mask_.end(); ++atom1)
    {
      double minDist = upperCutoff_;
      double dist2 = upperCutoff_;
      const double* XYZ1 = frm.Frm().XYZ(*atom1);
      for ( AtomMask::const_iterator atom2 = mask2_.begin(); atom2 != mask2_.end(); ++atom2)
      {
        const double* XYZ2 = frm.Frm().XYZ(*atom2);
        Matrix_3x3 ucell, recip;
        switch ( image_.ImageType() ) {
          case NONORTHO:
            frm.Frm().BoxCrd().ToRecip(ucell, recip);
            dist2 = DIST2_ImageNonOrtho(XYZ1, XYZ2, ucell, recip);
            break;
          case ORTHO:   dist2 = DIST2_ImageOrtho(XYZ1, XYZ2, frm.Frm().BoxCrd()); break;
          case NOIMAGE: dist2 = DIST2_NoImage(XYZ1, XYZ2); break;
        }
        // Find minimum distance.
        if (dist2 < minDist) {
          minDist = dist2;
          //mprintf("DEBUG:\tMinDist^2 %i {%g %g %g} to %i {%g %g %g} is %f\n",
          //        *atom1,XYZ1[0],XYZ1[1],XYZ1[2], *atom2,XYZ2[0],XYZ2[1],XYZ2[2], dist2);
        }
      }
      //mprintf("DEBUG: MinDist^2=%f\n", minDist);
      if (minDist > lowerCutoff_ && minDist < upperCutoff_) {
        nInside_[ *atom1 ] = 1;
        calculateMSD( frm.Frm().XYZ(*atom1), *atom1, *atom1, Box );
      }
    }
    // Calc average
    int Nin = 0;
    int i3 = 0;
    for (int i=0; i < n_atom_; ++i) {
      if ( nInside_[i] == 1 ) {
        average += distance_[i];
        avgx += distancexyz_[i3  ];
        avgy += distancexyz_[i3+1];
        avgz += distancexyz_[i3+2];
        ++Nin; // NOTE: nInside only ever gets set to 1
        //mprintf("DEBUG: %i nInside= %i d= %f dx= %f dy= %f dz= %f\n",
        //        i, Nin, distance_[i], distancexyz_[i3],distancexyz_[i3+1],distancexyz_[i3+2]);
      }
      i3 += 3;
    }
    if (Nin == 0) {
      mprinterr("Error: No atoms of mask 1 left for processing.\n");
      return Action::ERR;
    }
    average /= (double)Nin;
    avgx /= (double)Nin;
    avgy /= (double)Nin;
    avgz /= (double)Nin;
    outputnw_->Printf("%9.3f %7i\n", Time, Nin);
  }

  // Output
  output_->Printf("%10.3f %10.3f %10.3f %10.3f %10.3f",
                 Time, avgx, avgy, avgz, average);

  // Write un-imaged distances if requested.
  if (printDistances_) {
    int i3 = 0;
    for (int i = 0; i < mask_.Nselected(); ++i) {
      output_->Printf(" %9.3f %9.3f %9.3f %9.3f",
                     distancexyz_[i3  ], distancexyz_[i3+1], 
                     distancexyz_[i3+2], distance_[i]);
      i3 += 3;
    }
  }

  output_->Printf("\n");

  return Action::OK;
}

// Action_STFC_Diffusion::Print()
void Action_STFC_Diffusion::Print() {
  if (outputad_ != 0) {
    Darray::iterator d2 = dSum2_.begin();
    int i = 0;
    for (Darray::iterator d1 = dSum1_.begin(); d1 != dSum1_.end(); ++d1)
    {
      if (*d1 > 0.0) {
        double avMSD = (*d1 - (*d2 * *d2) / elapsedFrames_) / elapsedFrames_;
        outputad_->Printf("%7i %9.3f\n", i++, avMSD);
      }
      ++d2;
    }
  }
}
