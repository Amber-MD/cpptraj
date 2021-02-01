#include <cfloat> // DBL_MAX
#include <cmath> // sqrt
#include "Action_DNAionTracker.h"
#include "CpptrajStdio.h"
#include "DistRoutines.h"

// CONSTRUCTOR
Action_DNAionTracker::Action_DNAionTracker() :
  distance_(0),
  bintype_(COUNT),
  poffset_(0),
  useMass_(true)
{
  SetHidden(true);
}

void Action_DNAionTracker::Help() const {
  mprintf("\tname mask_p1 mask_p2 mask_base mask_ions\n"
          "\t[poffset <value>] [out <filename>] [time <interval>] [noimage]\n"
          "\t[shortest | counttopcone| countbottomcone | count]\n");
}

Action::RetType Action_DNAionTracker::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  // Get keywords
  DataFile* outfile = init.DFL().AddDataFile(actionArgs.GetStringKey("out"), actionArgs);
  poffset_ = actionArgs.getKeyDouble("poffset", 5.0);
  imageOpt_.InitImaging( !actionArgs.hasKey("noimage") );
  if (actionArgs.hasKey("shortest"))
    bintype_ = SHORTEST;
  else if (actionArgs.hasKey("counttopcone"))
    bintype_ = TOPCONE;
  else if (actionArgs.hasKey("countbottomcone"))
    bintype_ = BOTTOMCONE;
  else if (actionArgs.hasKey("count"))
    bintype_ = COUNT;

  // Get masks - 4 must be specified
  std::string m1 = actionArgs.GetMaskNext();
  std::string m2 = actionArgs.GetMaskNext();
  std::string m3 = actionArgs.GetMaskNext();
  std::string m4 = actionArgs.GetMaskNext();
  if (m1.empty() || m2.empty() || m3.empty() || m4.empty()) {
    mprinterr("Error: dnaiontracker requires 4 masks.\n");
    return Action::ERR;
  }
  if (p1_.SetMaskString(m1)) return Action::ERR;
  if (p2_.SetMaskString(m2)) return Action::ERR;
  if (base_.SetMaskString(m3)) return Action::ERR;
  if (ions_.SetMaskString(m4)) return Action::ERR;

  // Add dataset to dataset list (and datafile list if filename specified)
  distance_ = init.DSL().AddSet(DataSet::DOUBLE, MetaData(actionArgs.GetStringNext(),
                                                    MetaData::M_DISTANCE), "DNAion");
  if (distance_==0) return Action::ERR;
  if (outfile != 0)
    outfile->AddDataSet( distance_ );

  // INFO
  mprintf("Warning: DNAIONTRACKER is experimental code!\n");
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
  mprintf("      saved to array named %s\n", distance_->legend());
  mprintf("      Perpendicular offset for cone is %5.2f angstroms\n", poffset_);
  if (!imageOpt_.UseImage())
    mprintf("      Imaging has been disabled\n");
  mprintf("\tPhosphate1 Mask [%s]\n", p1_.MaskString());
  mprintf("\tPhosphate2 Mask [%s]\n", p2_.MaskString());
  mprintf("\tBase Mask       [%s]\n", base_.MaskString());
  mprintf("\tIons Mask       [%s]\n", ions_.MaskString());
  if (outfile != 0)
    mprintf("\tData will be printed to a file named %s\n", outfile->DataFilename().full());

  return Action::OK;
}

Action::RetType Action_DNAionTracker::Setup(ActionSetup& setup) {
  // Setup masks
  if (setup.Top().SetupIntegerMask( p1_ )) return Action::ERR; 
  if ( p1_.None() ) {
    mprinterr("Error: dnaiontracker: No atoms in mask1\n");
    return Action::ERR;
  }
  if (setup.Top().SetupIntegerMask( p2_ )) return Action::ERR;  
  if ( p2_.None() ) {
    mprinterr("Error: dnaiontracker: No atoms in mask2\n");
    return Action::ERR;
  }
  if (setup.Top().SetupIntegerMask( base_ )) return Action::ERR;  
  if ( base_.None() ) {
    mprinterr("Error: dnaiontracker: No atoms in mask3\n");
    return Action::ERR;
  }
  if (setup.Top().SetupIntegerMask( ions_ )) return Action::ERR;  
  if ( ions_.None() ) {
    mprinterr("Error: dnaiontracker: No atoms in mask4\n");
    return Action::ERR;
  }
  imageOpt_.SetupImaging( setup.CoordInfo().TrajBox().HasBox() );
  mprintf("\tPhosphate1 Mask [%s] %i atoms.\n", p1_.MaskString(), p1_.Nselected());
  mprintf("\tPhosphate2 Mask [%s] %i atoms.\n", p2_.MaskString(), p2_.Nselected());
  mprintf("\t      Base Mask [%s] %i atoms.\n", base_.MaskString(), base_.Nselected());
  mprintf("\t      Ions Mask [%s] %i atoms.\n", ions_.MaskString(), ions_.Nselected());

  return Action::OK;
}

Action::RetType Action_DNAionTracker::DoAction(int frameNum, ActionFrame& frm) {
  double d_tmp, dval;
  Vec3 P1, P2, BASE;

  if (imageOpt_.ImagingEnabled())
    imageOpt_.SetImageType( frm.Frm().BoxCrd().Is_X_Aligned_Ortho());
  // Get center for P1, P2, and Base
  if (useMass_) {
    P1 = frm.Frm().VCenterOfMass( p1_ );
    P2 = frm.Frm().VCenterOfMass( p2_ );
    BASE = frm.Frm().VCenterOfMass( base_ );
  } else {
    P1 = frm.Frm().VGeometricCenter( p1_ );
    P2 = frm.Frm().VGeometricCenter( p2_ );
    BASE = frm.Frm().VGeometricCenter( base_ );
  }
 
  // Calculate P -- P distance and centroid
  double d_pp = DIST2(imageOpt_.ImagingType(), P1, P2, frm.Frm().BoxCrd());
  Vec3 pp_centroid = (P1 + P2) / 2.0;

  // Cutoff^2
  double d_cut = d_pp*0.25 + (poffset_*poffset_); // TODO: precalc offset^2

  // Calculate P -- base centroid to median point
  double d_pbase = DIST2(imageOpt_.ImagingType(), pp_centroid, BASE, frm.Frm().BoxCrd());

  //double d_min = DBL_MAX;
  if (bintype_ == SHORTEST)
    dval = DBL_MAX; //d_min;
  else
    dval = 0;
  // Loop over ion positions
  for (AtomMask::const_iterator ion = ions_.begin(); ion != ions_.end(); ++ion)
  {
    const double* ionxyz = frm.Frm().XYZ(*ion);
    double d_p1ion =   DIST2(imageOpt_.ImagingType(), P1,   ionxyz, frm.Frm().BoxCrd());
    double d_p2ion =   DIST2(imageOpt_.ImagingType(), P2,   ionxyz, frm.Frm().BoxCrd());
    double d_baseion = DIST2(imageOpt_.ImagingType(), BASE, ionxyz, frm.Frm().BoxCrd());
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
  
  return Action::OK;
}
