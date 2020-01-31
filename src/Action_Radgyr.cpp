// Radius of Gyration 
#include <cmath> // sqrt
#include "Action_Radgyr.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_Radgyr::Action_Radgyr() :
  rog_(0),
  rogmax_(0),
  rogtensor_(0),
  calcRogmax_(true),
  calcTensor_(false),
  useMass_(false)
{ } 

void Action_Radgyr::Help() const {
  mprintf("\t[<name>] [<mask1>] [out <filename>] [mass] [nomax] [tensor]\n"
          "  Calculate radius of gyration of atoms in <mask>\n");
}

// Action_Radgyr::init()
Action::RetType Action_Radgyr::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  // Get keywords
  DataFile* outfile = init.DFL().AddDataFile(actionArgs.GetStringKey("out"), actionArgs);
  useMass_ = actionArgs.hasKey("mass");
  calcRogmax_ = !actionArgs.hasKey("nomax");
  calcTensor_ = actionArgs.hasKey("tensor");

  // Get Masks
  if (Mask1_.SetMaskString( actionArgs.GetMaskNext() )) return Action::ERR;

  // Datasets to store radius of gyration and max
  // Also add datasets to data file list
  rog_ = init.DSL().AddSet(DataSet::DOUBLE, actionArgs.GetStringNext(), "RoG");
  if (rog_==0) return Action::ERR;
  if (outfile != 0) outfile->AddDataSet( rog_ );
  if (calcRogmax_) {
    rogmax_ = init.DSL().AddSet(DataSet::DOUBLE, MetaData(rog_->Meta().Name(), "Max"));
    if (rogmax_ == 0) return Action::ERR; 
    if (outfile != 0) outfile->AddDataSet( rogmax_ );
  }
  if (calcTensor_) {
    rogtensor_ = init.DSL().AddSet(DataSet::VECTOR, MetaData(rog_->Meta().Name(), "Tensor"));
    if (rogtensor_ == 0) return Action::ERR;
    if (outfile != 0) outfile->AddDataSet( rogtensor_ );
  }

  mprintf("    RADGYR: Calculating for atoms in mask %s",Mask1_.MaskString());
  if (useMass_)
    mprintf(" using mass weighting");
  mprintf(".\n");
  if (!calcRogmax_)
    mprintf("\tRoG max will not be stored.\n");
  if (calcTensor_)
    mprintf("\tRoG tensor will also be calcd.\n");

  return Action::OK;
}

// Action_Radgyr::setup()
/** Set radius of gyration up for this parmtop. Get masks etc. */
// currentParm is set in Action::Setup
Action::RetType Action_Radgyr::Setup(ActionSetup& setup) {
  if ( setup.Top().SetupIntegerMask(Mask1_)) return Action::ERR;
  mprintf("\t%s (%i atoms).\n",Mask1_.MaskString(),Mask1_.Nselected());
  if (Mask1_.None()) {
    mprintf("Warning: Radgyr::setup: Mask contains 0 atoms.\n");
    return Action::SKIP;
  }

  return Action::OK;  
}

// CalcTensor()
static inline void CalcTensor(double* tsum, double dx, double dy, double dz, double mass) {
  double xx = dx * dx * mass;
  double yy = dy * dy * mass;
  double zz = dz * dz * mass;
  double xy = dx * dy * mass;
  double xz = dx * dz * mass;
  double yz = dy * dz * mass;

  tsum[0] += xx;
  tsum[1] += yy;
  tsum[2] += zz;
  tsum[3] += xy;
  tsum[4] += xz;
  tsum[5] += yz;
}

// Action_Radgyr::action()
/** Calc the radius of gyration of atoms in mask. Also record the maximum 
  * distance from center. Use center of mass if useMass is true.
  */
Action::RetType Action_Radgyr::DoAction(int frameNum, ActionFrame& frm) {
  double max = 0.0;
  double total_mass = 0.0;
  double maxMass = 1.0;
  double sumDist2 = 0.0;
  double tsum[6];
  tsum[0] = 0.0; tsum[1] = 0.0; tsum[2] = 0.0;
  tsum[3] = 0.0; tsum[4] = 0.0; tsum[5] = 0.0;

  // TODO: Make sumMass part of Frame?
  if (useMass_) {
    Vec3 mid = frm.Frm().VCenterOfMass( Mask1_ );
    for (AtomMask::const_iterator atom = Mask1_.begin(); atom != Mask1_.end(); ++atom)
    {
      const double* XYZ = frm.Frm().XYZ( *atom );
      double dx = XYZ[0] - mid[0];
      double dy = XYZ[1] - mid[1];
      double dz = XYZ[2] - mid[2];
      double mass = frm.Frm().Mass( *atom );
      if (calcTensor_) CalcTensor(tsum, dx, dy, dz, mass);
      total_mass += mass;
      double dist2 = ((dx*dx) + (dy*dy) + (dz*dz)) * mass;
      if (dist2 > max) {
        max = dist2;
        maxMass = mass;
      }
      sumDist2 += dist2;
    }
  } else {
    Vec3 mid = frm.Frm().VGeometricCenter( Mask1_ );
    total_mass = (double)Mask1_.Nselected();
    for (AtomMask::const_iterator atom = Mask1_.begin(); atom != Mask1_.end(); ++atom)
    {
      const double* XYZ = frm.Frm().XYZ( *atom );
      double dx = XYZ[0] - mid[0];
      double dy = XYZ[1] - mid[1];
      double dz = XYZ[2] - mid[2];
      if (calcTensor_) CalcTensor(tsum, dx, dy, dz, 1.0);
      double dist2 = ((dx*dx) + (dy*dy) + (dz*dz));
      if (dist2 > max)
        max = dist2;
      sumDist2 += dist2;
    }
  }

  if (total_mass == 0.0) {
    mprinterr("Error: radgyr: divide by zero.\n");
    return Action::ERR;
  }
  double Rog = sqrt( sumDist2 / total_mass );
  rog_->Add(frameNum, &Rog);
  if (calcRogmax_) {
    max = sqrt( max / maxMass );
    rogmax_->Add(frameNum, &max);
  }
  if (calcTensor_) {
    tsum[0] /= total_mass;
    tsum[1] /= total_mass;
    tsum[2] /= total_mass;
    tsum[3] /= total_mass;
    tsum[4] /= total_mass;
    tsum[5] /= total_mass;
    rogtensor_->Add(frameNum, tsum);
  }

  return Action::OK;
} 
