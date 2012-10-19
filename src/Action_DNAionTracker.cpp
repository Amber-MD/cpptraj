#include <cfloat> // DBL_MAX
#include <cmath> // sqrt
#include "Action_DNAionTracker.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_DNAionTracker::Action_DNAionTracker() :
  distance_(0),
  bintype_(COUNT),
  poffset_(0),
  useMass_(true)
{ }

void Action_DNAionTracker::Help() {

}

/** dnaiontracker name mask_p1 mask_p2 mask_base mask_ions 
  *               [poffset <value>] [out <filename>] [time <interval>] [noimage] 
  *               [shortest | counttopcone| countbottomcone | count]
  */
int Action_DNAionTracker::init() {
  // Get keywords
  std::string filename_ = actionArgs.GetStringKey("out");
  poffset_ = actionArgs.getKeyDouble("poffset", 5.0);
  InitImaging( !actionArgs.hasKey("noimage") );
  if (actionArgs.hasKey("shortest"))
    bintype_ = SHORTEST;
  else if (actionArgs.hasKey("counttopcone"))
    bintype_ = TOPCONE;
  else if (actionArgs.hasKey("countbottomcone"))
    bintype_ = BOTTOMCONE;
  else if (actionArgs.hasKey("count"))
    bintype_ = COUNT;

  // Get masks - 4 must be specified
  ArgList::ConstArg m1 = actionArgs.getNextMask();
  ArgList::ConstArg m2 = actionArgs.getNextMask();
  ArgList::ConstArg m3 = actionArgs.getNextMask();
  ArgList::ConstArg m4 = actionArgs.getNextMask();
  if (m1==NULL || m2==NULL || m3==NULL || m4==NULL) {
    mprinterr("Error: dnaiontracker requires 4 masks.\n");
    mprinterr("Error: mask1=%s  mask2=%s  mask3=%s  mask4=%s\n",m1,m2,m3,m4);
    return 1;
  }
  p1_.SetMaskString(m1);
  p2_.SetMaskString(m2);
  base_.SetMaskString(m3);
  ions_.SetMaskString(m4);

  // Dataset name
  ArgList::ConstArg dsetname = actionArgs.getNextString();
  if (dsetname==NULL) {
    mprinterr("Error: dnaiontracker: It is necessary to specify a unique name\n");
    mprinterr("Error:                for each specified tracking.\n");
    return 1;
  }

  // Add dataset to dataset list (and datafile list if filename specified)
  distance_ = DSL->Add(DataSet::DOUBLE, dsetname, "DNAion");
  // NOTE: Set to mode distance in PTRAJ
  distance_->SetScalar( DataSet::M_DISTANCE );
  if (distance_==NULL) return 1;
  if (!filename_.empty())
    DFL->Add( filename_.c_str(), distance_ );

  // INFO
  mprintf("    DNAIONTRACKER: Data representing the ");
  switch (bintype_) {
    case COUNT : 
      mprintf("count within the cone will be\n"); break;
    case SHORTEST: 
      mprintf("shortest distance to a phosphate or base centroid will be\n"); break;
    case TOPCONE: 
      mprintf("count in the top half of the cone (and sort-of bound) will be\n"); break;
    case BOTTOMCONE: 
      mprintf("count in the bottom half of the cone will be\n"); break;
  }
  mprintf("      saved to array named %s\n", distance_->Legend().c_str());
  mprintf("      Perpendicular offset for cone is %5.2f angstroms\n", poffset_);
  if (!UseImage())
    mprintf("      Imaging has been disabled\n");
  mprintf("\tPhosphate1 Mask [%s]\n", p1_.MaskString());
  mprintf("\tPhosphate2 Mask [%s]\n", p2_.MaskString());
  mprintf("\tBase Mask       [%s]\n", base_.MaskString());
  mprintf("\tIons Mask       [%s]\n", ions_.MaskString());
  if (!filename_.empty())
    mprintf("      Data will be printed to a file named %s\n", filename_.c_str());

  return 0;
}

int Action_DNAionTracker::setup() {
  // Setup masks
  if (currentParm->SetupIntegerMask( p1_ )) return 1; 
  if ( p1_.None() ) {
    mprinterr("Error: dnaiontracker: No atoms in mask1\n");
    return 1;
  }
  if (currentParm->SetupIntegerMask( p2_ )) return 1;  
  if ( p2_.None() ) {
    mprinterr("Error: dnaiontracker: No atoms in mask2\n");
    return 1;
  }
  if (currentParm->SetupIntegerMask( base_ )) return 1;  
  if ( base_.None() ) {
    mprinterr("Error: dnaiontracker: No atoms in mask3\n");
    return 1;
  }
  if (currentParm->SetupIntegerMask( ions_ )) return 1;  
  if ( ions_.None() ) {
    mprinterr("Error: dnaiontracker: No atoms in mask4\n");
    return 1;
  }
  SetupImaging( currentParm->BoxType() );
  mprintf("\tPhosphate1 Mask [%s] %i atoms.\n", p1_.MaskString(), p1_.Nselected());
  mprintf("\tPhosphate2 Mask [%s] %i atoms.\n", p2_.MaskString(), p2_.Nselected());
  mprintf("\t      Base Mask [%s] %i atoms.\n", base_.MaskString(), base_.Nselected());
  mprintf("\t      Ions Mask [%s] %i atoms.\n", ions_.MaskString(), ions_.Nselected());

  return 0;
}

int Action_DNAionTracker::action() {
  double ucell[9], recip[9], pp_centroid[3], P1[3], P2[3], BASE[3];
  double d_tmp, dval;
  Vec3 boxXYZ(currentFrame->BoxX(), currentFrame->BoxY(), currentFrame->BoxZ() );
  // Setup imaging info if necessary
  if (ImageType()==NONORTHO) 
    currentFrame->BoxToRecip(ucell,recip);

  // Get center for P1, P2, and Base
  if (useMass_) {
    currentFrame->CenterOfMass( P1, p1_ );
    currentFrame->CenterOfMass( P2, p2_ );
    currentFrame->CenterOfMass( BASE, base_ );
  } else {
    currentFrame->GeometricCenter( P1, p1_ );
    currentFrame->GeometricCenter( P2, p2_ );
    currentFrame->GeometricCenter( BASE, base_ );
  }
 
  // Calculate P -- P distance and centroid
  double d_pp = DIST2(P1, P2, ImageType(), boxXYZ, ucell, recip);
  pp_centroid[0] = (P1[0] + P2[0]) / 2;
  pp_centroid[1] = (P1[1] + P2[1]) / 2;
  pp_centroid[2] = (P1[2] + P2[2]) / 2;

  // Cutoff^2
  double d_cut = d_pp*0.25 + (poffset_*poffset_); // TODO: precalc offset^2

  // Calculate P -- base centroid to median point
  double d_pbase = DIST2(pp_centroid, BASE, ImageType(), boxXYZ, ucell, recip);

  //double d_min = DBL_MAX;
  if (bintype_ == SHORTEST)
    dval = DBL_MAX; //d_min;
  else
    dval = 0;
  // Loop over ion positions
  for (AtomMask::const_iterator ion = ions_.begin(); ion != ions_.end(); ++ion)
  {
    const double* ionxyz = currentFrame->XYZ(*ion);
    double d_p1ion =   DIST2(P1,   ionxyz, ImageType(), boxXYZ, ucell, recip);
    double d_p2ion =   DIST2(P2,   ionxyz, ImageType(), boxXYZ, ucell, recip);
    double d_baseion = DIST2(BASE, ionxyz, ImageType(), boxXYZ, ucell, recip);
    //mprintf("DEBUG: ion atom %i to P1 is %f\n", *ion+1, sqrt(d_p1ion));
    //mprintf("DEBUG: ion atom %i to P2 is %f\n", *ion+1, sqrt(d_p2ion));
    //mprintf("DEBUG: ion atom %i to base is %f\n", *ion+1, sqrt(d_baseion));
    //mprintf("DEBUG: d_pp is %f, poffset is %f, d_cut is %f\n", sqrt(d_pp), poffset_, sqrt(d_cut));

    int bound = 0;
    int boundLower = 0;
    int boundUpper = 0;

    if (d_p1ion < d_cut && d_p2ion < d_cut)
      bound = 1;
    if (d_baseion < d_pbase)
      boundLower = 1;
    if (bound && boundLower == 0)
      boundUpper = 1;
    //if (d_tmp > d_min)
    //  d_min = d_tmp;

    switch (bintype_) {
      case COUNT: 
        dval += (double)bound; break;
      case SHORTEST:
        if (d_p1ion < d_p2ion)
          d_tmp = d_p1ion;
        else
          d_tmp = d_p2ion;
        if (d_tmp > d_baseion)
          d_tmp = d_baseion;
        if (d_tmp < dval) 
          dval = d_tmp;
        break;
      case TOPCONE: 
        dval += (double)boundUpper; break;
      case BOTTOMCONE: 
        dval += (double)boundLower; break;
    }
  }
  if (bintype_ == SHORTEST)
    dval = sqrt(dval);
  distance_->Add(frameNum, &dval);
  
  return 0;
}
