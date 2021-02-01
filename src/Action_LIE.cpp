// Action_LIE
#include <cmath> // sqrt
#include "Action_LIE.h"
#include "CpptrajStdio.h"
#include "Constants.h" // ELECTOAMBER
#include "DistRoutines.h"

// CONSTRUCTOR
Action_LIE::Action_LIE() :
  dovdw_(false),
  doelec_(false),
  cut2vdw_(0.0),
  dielc_(1.0),
  cut2elec_(0.0),
  onecut2_(0.0)
{ }

void Action_LIE::Help() const {
  mprintf("\t[<name>] <Ligand Mask> [<Surroundings Mask>] [out filename] [nopbc]\n"
          "\t[noelec] [novdw] [cutvdw <cutoff>] [cutelec <cutoff>] [diel <dielc>]\n"
          "  Calculate linear interaction energy between <Ligand Mask> and <Surroundings Mask>\n");
}

// Action_LIE::init()
Action::RetType Action_LIE::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  // Get Keywords
  doelec_ = !(actionArgs.hasKey("noelec"));
  dovdw_ = !(actionArgs.hasKey("novdw"));
  usepbc_ = !(actionArgs.hasKey("nopbc"));
  DataFile* datafile = init.DFL().AddDataFile(actionArgs.GetStringKey("out"), actionArgs);
  dielc_ = actionArgs.getKeyDouble("diel", 1.0);
  double cut = actionArgs.getKeyDouble("cutvdw", 12.0);
  cut2vdw_ = cut * cut; // store square of cut for computational efficiency
  cut = actionArgs.getKeyDouble("cutelec", 12.0);
  cut2elec_ = cut * cut; // store square of cut for computational efficiency
  onecut2_ = 1 / cut2elec_;
  bool has_mask2 = false;

  // Use imaged distances unless requested otherwise
  imageOpt_.InitImaging(usepbc_);

  if (!doelec_ && !dovdw_) {
    mprinterr("Error: LIE: Cannot skip both ELEC and VDW calcs\n");
    return Action::ERR;
  }

  // Get Masks
  if (Mask1_.SetMaskString( actionArgs.GetMaskNext() )) return Action::ERR;
  std::string refmask = actionArgs.GetMaskNext();
  if (!refmask.empty()) {
    if (Mask2_.SetMaskString(refmask)) return Action::ERR;
    has_mask2 = true;
  }
  else {
    Mask2_ = Mask1_;
    Mask2_.InvertMaskExpression();
  }

  // Get data set name
  std::string ds_name = actionArgs.GetStringNext();
  if (ds_name.empty())
    ds_name = init.DSL().GenerateDefaultName("LIE");

  // Datasets
  if (doelec_) {
    elec_ = init.DSL().AddSet(DataSet::DOUBLE, MetaData(ds_name, "EELEC"));
    if (elec_ == 0) return Action::ERR;
    if (datafile != 0) datafile->AddDataSet(elec_);
  }
  if (dovdw_) {
    vdw_ = init.DSL().AddSet(DataSet::DOUBLE, MetaData(ds_name, "EVDW"));
    if (vdw_ == 0) return Action::ERR;
    if (datafile != 0) datafile->AddDataSet(vdw_);
  }

  mprintf("    LIE: Ligand mask is %s. Surroundings are ", Mask1_.MaskString());
  if (!has_mask2)
    mprintf("everything else. ");
  else
    mprintf("atoms in mask %s. ", Mask2_.MaskString());
  mprintf("Cutoff is %.3lf Ang. ", cut);
  if (!doelec_)
    mprintf("Skipping Electrostatic Calc. ");
  if (!dovdw_)
    mprintf("Skipping VDW Calc. ");
  if (usepbc_)
    mprintf("Using PBC.");
  else
    mprintf("NOT using PBC. Make sure the trajectory is properly imaged.");
  mprintf("\n");

  return Action::OK;
}

// Action_LIE::setup()
/** Determine what atoms each mask pertains to for the current parm file.
  */
Action::RetType Action_LIE::Setup(ActionSetup& setup) {
  if (setup.Top().SetupIntegerMask( Mask1_ )) return Action::ERR;
  if (setup.Top().SetupIntegerMask( Mask2_ )) return Action::ERR;

  mprintf("\tLIE: %i Ligand Atoms, %i Surrounding Atoms\n",
          Mask1_.Nselected(), Mask2_.Nselected());

  if (!setup.CoordInfo().TrajBox().HasBox()) {
    mprinterr("Error: LIE: Must have explicit solvent system with box info\n");
    return Action::ERR;
  }
  imageOpt_.SetupImaging( setup.CoordInfo().TrajBox().HasBox() );

  if (Mask1_.None() || Mask2_.None()) {
    mprintf("Warning: LIE: One or both masks have no atoms.\n");
    return Action::SKIP;
  }

  if (SetupParms(setup.Top()))
    return Action::ERR;

  // Back up the parm
  CurrentParm_ = setup.TopAddress();

  return Action::OK;
}

// Action_LIE::SetupParms
/** Sets the temporary charge array and makes sure that we have the necessary
  * parameters in our topology to calculate nonbonded energy terms
  */
int Action_LIE::SetupParms(Topology const& ParmIn) {
  if (!ParmIn.Nonbond().HasNonbond()) {
    mprinterr("Error: Topology does not have LJ information.\n");
    return 1;
  }
  // Store the charges
  atom_charge_.clear();
  atom_charge_.reserve( ParmIn.Natom() );
  for (Topology::atom_iterator atom = ParmIn.begin();
                               atom != ParmIn.end(); ++atom)
    atom_charge_.push_back( atom->Charge() * Constants::ELECTOAMBER / sqrt(dielc_) );
  return 0;
}

double Action_LIE::Calculate_LJ(Frame const& frameIn, Topology const& parmIn) const {
  double result = 0;
  // Loop over ligand atoms
  AtomMask::const_iterator mask1_end = Mask1_.end();
  AtomMask::const_iterator mask2_end = Mask2_.end();
  for (AtomMask::const_iterator maskatom1 = Mask1_.begin();
       maskatom1 != mask1_end; maskatom1++) {

    int crdidx1 = (*maskatom1) * 3; // index into coordinate array
    const double* atm1 = frameIn.CRD(crdidx1);

    for (AtomMask::const_iterator maskatom2 = Mask2_.begin();
         maskatom2 != mask2_end; maskatom2++) {

      int crdidx2 = (*maskatom2) * 3; // index into coordinate array
      const double* atm2 = frameIn.CRD(crdidx2);

      // Get imaged distance
      double dist2 = DIST2(imageOpt_.ImagingType(), atm1, atm2, frameIn.BoxCrd());

      if (dist2 > cut2vdw_) continue;
      // Here we add to our nonbonded (VDW) energy
      NonbondType const& LJ = parmIn.GetLJparam(*maskatom1, *maskatom2);
      double r2 = 1 / dist2;
      double r6 = r2 * r2 * r2;
      result += LJ.A() * r6 * r6 - LJ.B() * r6;
    }
  }

  return result;
}

double Action_LIE::Calculate_Elec(Frame const& frameIn) const {
  double result = 0;
  // Loop over ligand atoms
  AtomMask::const_iterator mask1_end = Mask1_.end();
  AtomMask::const_iterator mask2_end = Mask2_.end();
  for (AtomMask::const_iterator maskatom1 = Mask1_.begin();
       maskatom1 != mask1_end; maskatom1++) {

    int crdidx1 = (*maskatom1) * 3; // index into coordinate array
    const double* atm1 = frameIn.CRD(crdidx1);

    for (AtomMask::const_iterator maskatom2 = Mask2_.begin();
         maskatom2 != mask2_end; maskatom2++) {

      int crdidx2 = (*maskatom2) * 3; // index into coordinate array
      const double* atm2 = frameIn.CRD(crdidx2);

      // Get imaged distance
      double dist2 = DIST2(imageOpt_.ImagingType(), atm1, atm2, frameIn.BoxCrd());

      if (dist2 > cut2elec_) continue;
      // Here we add to our electrostatic energy
      double qiqj = atom_charge_[*maskatom1] * atom_charge_[*maskatom2];
      double shift = (1 - dist2 * onecut2_);
      result += qiqj / sqrt(dist2) * shift * shift;
    }
  }

  return result;
}

// Action_LIE::action()
Action::RetType Action_LIE::DoAction(int frameNum, ActionFrame& frm) {
  if (imageOpt_.ImagingEnabled()) {
    imageOpt_.SetImageType( frm.Frm().BoxCrd().Is_X_Aligned_Ortho() );
    //mprintf("DEBUG: Image type is %i\n", (int)imageOpt_.ImagingType());
  }
  if (doelec_) {
    double e = Calculate_Elec(frm.Frm());
    elec_->Add(frameNum, &e);
  }

  if (dovdw_) {
    double e = Calculate_LJ(frm.Frm(), *CurrentParm_);
    vdw_->Add(frameNum, &e);
  }

  return Action::OK;
}
