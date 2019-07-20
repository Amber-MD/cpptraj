#include "Action_Pucker.h"
#include "CpptrajStdio.h"
#include "Constants.h" // RADDEG
#include "TorsionRoutines.h"

// CONSTRUCTOR
Action_Pucker::Action_Pucker() :
  pucker_(0),
  amplitude_(0),
  theta_(0),
  puckerMin_(0.0),
  puckerMax_(0.0),
  offset_(0.0),
  puckerMethod_(ALTONA),
  useMass_(true)
{ } 

void Action_Pucker::Help() const {
  mprintf("\t[<name>] <mask1> <mask2> <mask3> <mask4> <mask5> [<mask6>] [geom]\n"
          "\t[out <filename>] [altona | cremer] [amplitude] [theta]\n"
          "\t[range360] [offset <offset>]\n"
          "\tCalculate pucker of atoms in masks 1-5 (or 6, 'cremer' only).\n");
}

// Action_Pucker::Init()
Action::RetType Action_Pucker::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  // Get keywords
  DataFile* outfile = init.DFL().AddDataFile( actionArgs.GetStringKey("out"), actionArgs);
  if      (actionArgs.hasKey("altona")) puckerMethod_ = ALTONA;
  else if (actionArgs.hasKey("cremer")) puckerMethod_ = CREMER;
  else                                  puckerMethod_ = UNSPECIFIED;
  bool calc_amp = actionArgs.hasKey("amplitude");
  bool calc_theta = actionArgs.hasKey("theta");
  offset_ = actionArgs.getKeyDouble("offset",0.0);
  if (actionArgs.hasKey("range360"))
    puckerMin_ = 0.0;
  else
    puckerMin_ = -180.0;
  puckerMax_ = puckerMin_ + 360.0;
  useMass_ = !actionArgs.hasKey("geom");
  // DEPRECATED - 'type pucker', now always set.
  std::string stypename = actionArgs.GetStringKey("type");
  if (!stypename.empty())
    mprintf("Warning: 'type' keyword is no longer necessary for 'pucker' and will be\n"
            "Warning:   deprecated in the future.\n");
  // Get Masks
  Masks_.clear();
  std::string mask_expression = actionArgs.GetMaskNext();
  while (!mask_expression.empty()) {
    Masks_.push_back( AtomMask( mask_expression ) );
    mask_expression = actionArgs.GetMaskNext();
  }
  if (Masks_.size() < 5 || Masks_.size() > 6) {
    mprinterr("Error: Pucker can only be calculated for 5 or 6 masks, %zu specified.\n",
              Masks_.size());
    return Action::ERR;
  }
  // Choose default method
  if (puckerMethod_ == UNSPECIFIED) {
    mprintf("Warning: Pucker method not specified");
    if (Masks_.size() == 6) {
      mprintf(" and 6 masks present, defaulting to 'cremer'.\n");
      puckerMethod_ = CREMER;
    } else {
      mprintf(", defaulting to 'altona'.\n");
      puckerMethod_ = ALTONA;
    }
  }
  // 6 masks only supported by cremer method right now
  if (Masks_.size() == 6 && puckerMethod_ != CREMER) {
    mprinterr("Error: Pucker with 6 masks only supported with 'cremer'\n");
    return Action::ERR;
  }
  // Set up array to hold coordinate vectors.
  AX_.resize( Masks_.size() );

  // Setup dataset
  MetaData md(actionArgs.GetStringNext());
  md.SetScalarMode( MetaData::M_PUCKER );
  md.SetScalarType( MetaData::PUCKER );
  pucker_ = init.DSL().AddSet(DataSet::DOUBLE, md, "Pucker");
  if (pucker_ == 0) return Action::ERR;
  amplitude_ = 0;
  theta_ = 0;
  if (calc_amp)
    amplitude_ = init.DSL().AddSet(DataSet::DOUBLE, MetaData(pucker_->Meta().Name(), "Amp"));
  if (calc_theta) {
    if ( Masks_.size() < 6 )
      mprintf("Warning: 'theta' calc. not supported for < 6 masks.\n");
    else
      theta_ = init.DSL().AddSet(DataSet::DOUBLE, MetaData(pucker_->Meta().Name(), "Theta"));
  }
  // Add dataset to datafile list
  if (outfile != 0) {
    outfile->AddDataSet( pucker_ );
    if (amplitude_ != 0) outfile->AddDataSet( amplitude_ );
    if (theta_ != 0) outfile->AddDataSet( theta_ );
  }

  mprintf("    PUCKER: ");
  for (std::vector<AtomMask>::const_iterator MX = Masks_.begin();
                                             MX != Masks_.end(); ++MX)
  {
    if (MX != Masks_.begin()) mprintf("-");
    mprintf("[%s]", MX->MaskString());
  }
  mprintf("\n");
  if (puckerMethod_==ALTONA) 
    mprintf("\tUsing Altona & Sundaralingam method.\n");
  else if (puckerMethod_==CREMER)
    mprintf("\tUsing Cremer & Pople method.\n");
  if (outfile != 0) 
    mprintf("\tData will be written to %s\n", outfile->DataFilename().base());
  if (amplitude_!=0)
    mprintf("\tAmplitudes (in degrees) will be stored.\n");
  if (theta_!=0)
    mprintf("\tThetas (in degrees) will be stored.\n");
  if (offset_!=0)
    mprintf("\tOffset: %f degrees will be added to values.\n", offset_);
  if (puckerMin_ > -180.0)
    mprintf("\tOutput range is 0 to 360 degrees.\n");
  else
    mprintf("\tOutput range is -180 to 180 degrees.\n");

  return Action::OK;
}

// Action_Pucker::Setup()
Action::RetType Action_Pucker::Setup(ActionSetup& setup) {
  mprintf("\t");
  for (std::vector<AtomMask>::iterator MX = Masks_.begin();
                                       MX != Masks_.end(); ++MX)
  {
    if ( setup.Top().SetupIntegerMask( *MX ) ) return Action::ERR;
    MX->BriefMaskInfo();
    if (MX->None()) {
      mprintf("\nWarning: Mask '%s' selects no atoms for topology '%s'\n",
              MX->MaskString(), setup.Top().c_str());
      return Action::SKIP;
    }
  }
  mprintf("\n");

  return Action::OK;  
}

// Action_Pucker::DoAction()
Action::RetType Action_Pucker::DoAction(int frameNum, ActionFrame& frm) {
  double pval, aval, tval;
  std::vector<Vec3>::iterator ax = AX_.begin(); 

  if (useMass_) {
    for (std::vector<AtomMask>::const_iterator MX = Masks_.begin();
                                               MX != Masks_.end(); ++MX, ++ax)
      *ax = frm.Frm().VCenterOfMass( *MX );
  } else {
     for (std::vector<AtomMask>::const_iterator MX = Masks_.begin();
                                               MX != Masks_.end(); ++MX, ++ax)
      *ax = frm.Frm().VGeometricCenter( *MX );
  }

  switch (puckerMethod_) {
    case ALTONA: 
      pval = Pucker_AS( AX_[0].Dptr(), AX_[1].Dptr(), AX_[2].Dptr(), 
                        AX_[3].Dptr(), AX_[4].Dptr(), aval );
      break;
    case CREMER:
      pval = Pucker_CP( AX_[0].Dptr(), AX_[1].Dptr(), AX_[2].Dptr(), 
                        AX_[3].Dptr(), AX_[4].Dptr(), AX_[5].Dptr(), 
                        AX_.size(), aval, tval );
      break;
    case UNSPECIFIED : // Sanity check
      return Action::ERR;
  }
  if ( amplitude_ != 0 ) {
    aval *= Constants::RADDEG;
    amplitude_->Add(frameNum, &aval);
  }
  if ( theta_ != 0 ) {
    tval *= Constants::RADDEG;
    theta_->Add(frameNum, &tval);
  }
  pval = (pval * Constants::RADDEG) + offset_;
  if (pval > puckerMax_)
    pval -= 360.0;
  else if (pval < puckerMin_)
    pval += 360.0;
  pucker_->Add(frameNum, &pval);

  return Action::OK;
} 
