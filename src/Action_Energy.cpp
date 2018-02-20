#include "Action_Energy.h"
#include "CpptrajStdio.h"
#include "Ewald_Regular.h"
#include "Ewald_ParticleMesh.h"

/// CONSTRUCTOR
Action_Energy::Action_Energy() : currentParm_(0), debug_(0), EW_(0)
{
  std::fill(mlimits_, mlimits_+3, 0);
}

/// DESTRUCTOR
Action_Energy::~Action_Energy() {
  if (EW_ != 0) delete EW_;
}

void Action_Energy::Help() const {
  mprintf("\t[<name>] [<mask1>] [out <filename>]\n"
          "\t[bond] [angle] [dihedral] [kinetic [dt <dt>]] {[nb14] | [e14] | [v14]}\n"
          "\t{[nonbond] | [elec] [vdw]}\n"
          "\t[ etype { simple |\n"
          "\t          directsum [npoints <N>] |\n"
          "\t          ewald [cut <cutoff>] [dsumtol <dtol>] [rsumtol <rtol>]\n"
          "\t                [ewcoeff <coeff>] [maxexp <max>] [skinnb <skinnb>]\n"
          "\t                [mlimits <X>,<Y>,<Z>] [erfcdx <dx>]\n"
          "\t          pme [cut <cutoff>] [dsumtol <dtol>] [order <order>]\n"
          "\t              [ewcoeff <coeff>] [skinnb <skinnb>]\n"
          "\t              [nfft <nfft1>,<nfft2>,<nfft3>] [erfcdx <dx>]\n"
          "\t        } ]\n"
          "  Calculate energy for atoms in mask.\n");
}

/// Corresponds to Etype 
static const char* AspectStr[] = {"bond", "angle", "dih", "vdw14", "elec14",
                                  "vdw", "elec", "kinetic", "total"};

/// Corresponds to Etype
static const char* EtypeStr[] = {"Bonds", "Angles", "Dihedrals", "1-4 VDW", "1-4 Elec.",
                                 "VDW", "Elec.", "Kinetic", "Total"};

/// Corresponds to ElecType 
static const char* ElecStr[] = { "None", "Simple", "Direct Sum", "Regular Ewald",
                                 "Particle Mesh Ewald" };

// Action_Energy::AddSet()
int Action_Energy::AddSet(Etype typeIn, DataSetList& DslIn, DataFile* outfile,
                          std::string const& setnameIn)
{
  Energy_[typeIn] = DslIn.AddSet(DataSet::DOUBLE, MetaData(setnameIn, AspectStr[typeIn]));
  if (Energy_[typeIn] == 0) return 1;
  if (outfile != 0) outfile->AddDataSet( Energy_[typeIn] );
  return 0;
}

// Action_Energy::Init()
Action::RetType Action_Energy::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  debug_ = debugIn;
  ENE_.SetDebug( debug_ );
  DataFile* outfile = init.DFL().AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  // Determine which energy terms are active 
  std::vector<bool> termEnabled((int)TOTAL+1, false);
  termEnabled[BOND] = actionArgs.hasKey("bond");
  termEnabled[ANGLE] = actionArgs.hasKey("angle");
  termEnabled[DIHEDRAL] = actionArgs.hasKey("dihedral");
  termEnabled[V14] = actionArgs.hasKey("v14");
  termEnabled[Q14] = actionArgs.hasKey("e14");
  termEnabled[VDW] = actionArgs.hasKey("vdw");
  termEnabled[ELEC] = actionArgs.hasKey("elec");
  if (actionArgs.hasKey("nb14")) {
    termEnabled[V14] = true;
    termEnabled[Q14] = true;
  }
  if (actionArgs.hasKey("nonbond")) {
    termEnabled[VDW] = true;
    termEnabled[ELEC] = true;
  }
  termEnabled[KE] = actionArgs.hasKey("kinetic");
  int nactive = 0;
  for (std::vector<bool>::const_iterator it = termEnabled.begin(); it != termEnabled.end(); ++it)
    if (*it) ++nactive;
  // If no terms specified, enabled everything. TODO disable KE?
  if (nactive == 0) termEnabled.assign((int)TOTAL+1, true);
  // If more than one term enabled ensure total will be calculated.
  if (nactive > 1) termEnabled[TOTAL] = true;
  // If KE enabled get time step
  if (termEnabled[KE])
    dt_ = actionArgs.getKeyDouble("dt", 0.001);
  // Electrostatics type.
  std::string etypearg = actionArgs.GetStringKey("etype");
  elecType_ = NO_ELE;
  EW_ = 0;
  if (!etypearg.empty()) {
    termEnabled[ELEC] = true;
    if (etypearg == "directsum") {
      // Direct sum method
      elecType_ = DIRECTSUM;
      npoints_ = actionArgs.getKeyInt("npoints", 0);
    } else if (etypearg == "ewald") {
      // Ewald method
      elecType_ = EWALD;
      cutoff_ = actionArgs.getKeyDouble("cut", 8.0);
      dsumtol_ = actionArgs.getKeyDouble("dsumtol", 1E-5);
      rsumtol_ = actionArgs.getKeyDouble("rsumtol", 5E-5);
      ewcoeff_ = actionArgs.getKeyDouble("ewcoeff", 0.0);
      maxexp_ = actionArgs.getKeyDouble("maxexp", 0.0);
      skinnb_ = actionArgs.getKeyDouble("skinnb", 2.0);
      erfcDx_ = actionArgs.getKeyDouble("erfcdx", 0.0);
      std::string marg = actionArgs.GetStringKey("mlimits");
      if (!marg.empty()) {
        ArgList mlim(marg, ",");
        if (mlim.Nargs() != 3) {
          mprinterr("Error: Need 3 integers in comma-separated list for 'mlimits'\n");
          return Action::ERR;
        }
        mlimits_[0] = mlim.getNextInteger(0);
        mlimits_[1] = mlim.getNextInteger(0);
        mlimits_[2] = mlim.getNextInteger(0);
      } else
        std::fill(mlimits_, mlimits_+3, 0);
      EW_ = (Ewald*)new Ewald_Regular();
    } else if (etypearg == "pme") {
      // particle mesh Ewald method
#     ifdef LIBPME
      elecType_ = PME;
      cutoff_ = actionArgs.getKeyDouble("cut", 8.0);
      dsumtol_ = actionArgs.getKeyDouble("dsumtol", 1E-5);
      ewcoeff_ = actionArgs.getKeyDouble("ewcoeff", 0.0);
      skinnb_ = actionArgs.getKeyDouble("skinnb", 2.0);
      erfcDx_ = actionArgs.getKeyDouble("erfcdx", 0.0);
      npoints_ = actionArgs.getKeyInt("order", 6);
      std::string marg = actionArgs.GetStringKey("nfft");
      if (!marg.empty()) {
        ArgList mlim(marg, ",");
        if (mlim.Nargs() != 3) {
          mprinterr("Error: Need 3 integers in comma-separated list for 'nfft'\n");
          return Action::ERR;
        }
        mlimits_[0] = mlim.getNextInteger(0);
        mlimits_[1] = mlim.getNextInteger(0);
        mlimits_[2] = mlim.getNextInteger(0);
      } else
        std::fill(mlimits_, mlimits_+3, -1);
      EW_ = (Ewald*)new Ewald_ParticleMesh();
#     else
      mprinterr("Error: 'pme' requires compiling with FFTW3 and C++11 support.\n");
      return Action::ERR;
#     endif
    } else if (etypearg == "simple") {
      // Simple method
      elecType_ = SIMPLE;
    } else {
      mprinterr("Error: Unrecognized option for 'etype': %s\n", etypearg.c_str());
      return Action::ERR;
    }
  }
  // If electrostatics enabled but type not specified, default to SIMPLE
  if (termEnabled[ELEC] && elecType_ == NO_ELE) elecType_ = SIMPLE;
  // Set up calculations
  Ecalcs_.clear();
  if (termEnabled[BOND])
    Ecalcs_.push_back(C_BND);
  if (termEnabled[ANGLE])
    Ecalcs_.push_back(C_ANG);
  if (termEnabled[DIHEDRAL])
    Ecalcs_.push_back(C_DIH);
  if (termEnabled[KE]) {
    calc_ke_ = true;
    Ecalcs_.push_back(C_KE);
  } else
    calc_ke_ = false;
  if (termEnabled[V14] || termEnabled[Q14])
    Ecalcs_.push_back(C_N14);
  // Determine which nonbonded calc to use if any.
  bool lj_longrange_correction = false;
  need_lj_params_ = false;
  if (termEnabled[ELEC] || termEnabled[VDW]) {
    // NOTE: if elecType_ is not NO_ELE then ELEC term is enabled by default
    if (elecType_ == SIMPLE) {
      if (termEnabled[ELEC] && termEnabled[VDW]) {
        Ecalcs_.push_back(C_NBD);
        need_lj_params_ = true;
      } else if (termEnabled[ELEC] && !termEnabled[VDW]) {
        Ecalcs_.push_back(C_COULOMB);
      }
    } else if (elecType_ == DIRECTSUM) {
      if (termEnabled[ELEC] && termEnabled[VDW]) {
        Ecalcs_.push_back(C_LJ);
        Ecalcs_.push_back(C_DIRECT);
        need_lj_params_ = true;
      } else if (termEnabled[ELEC] && !termEnabled[VDW]) {
        Ecalcs_.push_back(C_DIRECT);
      }
    } else if (elecType_ == EWALD) {
      Ecalcs_.push_back(C_EWALD);
      lj_longrange_correction = true;
      need_lj_params_ = true;
    } else if (elecType_ == PME) {
      Ecalcs_.push_back(C_PME);
      lj_longrange_correction = true;
      need_lj_params_ = true;
    } else if (elecType_ == NO_ELE) {
      Ecalcs_.push_back(C_LJ);
      need_lj_params_ = true;
    }
  }

  // Get Masks
  Mask1_.SetMaskString( actionArgs.GetMaskNext() );

  // DataSet
  setname_ = actionArgs.GetStringNext();
  if (setname_.empty())
    setname_ = init.DSL().GenerateDefaultName("ENE");
  Energy_.clear();
  Energy_.resize( (int)TOTAL + 1, 0 );
  for (int i = 0; i != (int)TOTAL+1; i++)
  {
    if (termEnabled[i]) {
      if (AddSet((Etype)i, init.DSL(), outfile, setname_)) return Action::ERR;
    }
  }
      
  mprintf("    ENERGY: Calculating energy for atoms in mask '%s'\n", Mask1_.MaskString());
  mprintf("\tCalculating terms:");
  for (int i = 0; i != (int)TOTAL+1; i++)
    if (termEnabled[i]) mprintf(" '%s'", EtypeStr[i]);
  mprintf("\n");
  if (elecType_ != NO_ELE)
    mprintf("\tElectrostatics method: %s\n", ElecStr[elecType_]);
  if (elecType_ == DIRECTSUM) {
    if (npoints_ < 0)
      mprintf("\tDirect sum energy for up to %i unit cells in each direction will be calculated.\n",
              -npoints_);
    else
      mprintf("\tDirect sum energy for %i unit cells in each direction will be calculated.\n",
              npoints_);
  } else if (elecType_ == EWALD) {
    mprintf("\tDirect space cutoff= %.4f\n", cutoff_);
    if (dsumtol_ != 0.0)
      mprintf("\tDirect sum tolerance= %g\n", dsumtol_);
    if (rsumtol_ != 0.0)
      mprintf("\tReciprocal sum tolerance= %g\n", rsumtol_);
    if (ewcoeff_ == 0.0)
      mprintf("\tWill determine Ewald coefficient from cutoff and direct sum tolerance.\n");
    else
      mprintf("\tEwald coefficient= %.4f\n", ewcoeff_);
    if (maxexp_ == 0.0)
      mprintf("\tWill determine MaxExp from Ewald coefficient and direct sum tolerance.\n");
    else
      mprintf("\tMaxExp= %g\n", maxexp_);
    if (mlimits_[0] < 1 && mlimits_[1] < 1 && mlimits_[2] < 1)
      mprintf("\tWill determine number of reciprocal vectors from MaxExp.\n");
    else
      mprintf("\tNumber of reciprocal vectors in each direction= {%i,%i,%i}\n",
              mlimits_[0], mlimits_[1], mlimits_[2]);
    if (erfcDx_ > 0.0)
      mprintf("\tERFC table dx= %g\n", erfcDx_);
  } else if (elecType_ == PME) {
    mprintf("\tDirect space cutoff= %.4f\n", cutoff_);
    if (dsumtol_ != 0.0)
      mprintf("\tDirect sum tolerance= %g\n", dsumtol_);
    mprintf("\tSpline order= %i\n", npoints_);
    if (ewcoeff_ == 0.0)
      mprintf("\tWill determine Ewald coefficient from cutoff and direct sum tolerance.\n");
    else
      mprintf("\tEwald coefficient= %.4f\n", ewcoeff_);
    if (mlimits_[0] < 1 && mlimits_[1] < 1 && mlimits_[2] < 1)
      mprintf("\tWill determine number of FFT grid points from box size.\n");
    else
      mprintf("\tNumber of FFT grid points in each direction= {%i,%i,%i}\n",
              mlimits_[0], mlimits_[1], mlimits_[2]);
    if (erfcDx_ > 0.0)
      mprintf("\tERFC table dx= %g\n", erfcDx_);
  }
  if (termEnabled[VDW] && lj_longrange_correction)
    mprintf("\tUsing long range correction for nonbond VDW calc.\n");
  if (calc_ke_)
    mprintf("\tTime step for KE calculation if forces present: %g ps\n", dt_);

  return Action::OK;
}

// Action_Energy::Setup()
/** Set angle up for this parmtop. Get masks etc.
  */
Action::RetType Action_Energy::Setup(ActionSetup& setup) {
  if (setup.Top().SetupCharMask(Mask1_)) return Action::ERR;
  if (Mask1_.None()) {
    mprintf("Warning: Mask '%s' selects no atoms.\n", Mask1_.MaskString());
    return Action::SKIP;
  }
  Mask1_.MaskInfo();
  Imask_ = AtomMask(Mask1_.ConvertToIntMask(), Mask1_.Natom());
  // Check for LJ terms
  if (need_lj_params_ && !setup.Top().Nonbond().HasNonbond())
  {
    mprinterr("Error: LJ energy calc requested but topology '%s'\n"
              "Error:   does not have LJ parameters.\n", setup.Top().c_str());
    return Action::ERR;
  }
  // Set up Ewald if necessary.
  if (elecType_ == EWALD) {
    if (((Ewald_Regular*)EW_)->Init(setup.CoordInfo().TrajBox(), cutoff_, dsumtol_, rsumtol_,
                                    ewcoeff_, maxexp_, skinnb_, erfcDx_, debug_, mlimits_))
      return Action::ERR;
    EW_->Setup( setup.Top(), Imask_ );
  }
  // For KE, check for velocities/forces
  if (calc_ke_) {
    if (!setup.CoordInfo().HasVel()) {
      mprintf("Warning: Coordinates have no velocities - kinetic energy will be zero.\n");
    } else {
      if (setup.CoordInfo().HasForce())
        mprintf("\tForce info present. Assuming plus-half time step velocities.\n"
                "\tVelocities at time 't' will be estimated using force info.\n");
      else
        mprintf("\tForce info not present. Assuming velocities are at same time\n"
                "\tstep as coordinates.\n");
    }
  }
# ifdef LIBPME
  else if (elecType_ == PME) {
    if (((Ewald_ParticleMesh*)EW_)->Init(setup.CoordInfo().TrajBox(), cutoff_, dsumtol_,
                                         ewcoeff_, skinnb_, erfcDx_, npoints_, debug_, mlimits_))
      return Action::ERR;
    EW_->Setup( setup.Top(), Imask_ );
  }
# endif
  currentParm_ = setup.TopAddress();
  return Action::OK;
}

/// For debugging the direct sum convergence
double Action_Energy::Dbg_Direct(Frame const& frameIn, int maxpoints) {
  // DEBUG
  double lastEQ = 0.0;
  for (int npoints = 0; npoints < maxpoints; npoints++) {
    double EQ = ENE_.E_DirectSum(frameIn, *currentParm_, Imask_, npoints);
    mprintf("DEBUG: %i points DirectSum= %12.4f", npoints, EQ);
    if (npoints > 0) {
      mprintf(" delta= %g", EQ - lastEQ);
    }
    mprintf("\n");
    lastEQ = EQ;
  }
  return lastEQ;
}

// Action_Energy::DoAction()
Action::RetType Action_Energy::DoAction(int frameNum, ActionFrame& frm) {
  etime_.Start();
  double Etot = 0.0, ene, ene2;
  typedef std::vector<CalcType>::const_iterator calc_it;
  for (calc_it calc = Ecalcs_.begin(); calc != Ecalcs_.end(); ++calc)
  {
    switch (*calc) {
      case C_BND:
        ene = ENE_.E_bond(frm.Frm(), *currentParm_, Mask1_);
        Energy_[BOND]->Add(frameNum, &ene);
        Etot += ene;
        break;
      case C_ANG:
        ene = ENE_.E_angle(frm.Frm(), *currentParm_, Mask1_);
        Energy_[ANGLE]->Add(frameNum, &ene);
        Etot += ene;
        break;
      case C_DIH:
        ene = ENE_.E_torsion(frm.Frm(), *currentParm_, Mask1_);
        Energy_[DIHEDRAL]->Add(frameNum, &ene);
        Etot += ene;
        break;
      case C_N14:
        ene = ENE_.E_14_Nonbond(frm.Frm(), *currentParm_, Mask1_, ene2);
        if (Energy_[V14] != 0) Energy_[V14]->Add(frameNum, &ene);
        if (Energy_[Q14] != 0) Energy_[Q14]->Add(frameNum, &ene2);
        Etot += (ene + ene2);
        break;
      case C_NBD: // Both nonbond terms must be enabled 
        ene = ENE_.E_Nonbond(frm.Frm(), *currentParm_, Imask_, ene2);
        Energy_[VDW]->Add(frameNum, &ene);
        Energy_[ELEC]->Add(frameNum, &ene2);
        Etot += (ene + ene2);
        break;
      case C_LJ:
        ene = ENE_.E_VDW(frm.Frm(), *currentParm_, Imask_);
        Energy_[VDW]->Add(frameNum, &ene);
        Etot += ene;
        break;
      case C_COULOMB:
        ene = ENE_.E_Elec(frm.Frm(), *currentParm_, Imask_);
        Energy_[ELEC]->Add(frameNum, &ene);
        Etot += ene;
        break;
      case C_DIRECT:
        if (npoints_ < 0)
          ene = Dbg_Direct(frm.Frm(), (-npoints_)+1);
        else
          ene = ENE_.E_DirectSum(frm.Frm(), *currentParm_, Imask_, npoints_);
        Energy_[ELEC]->Add(frameNum, &ene);
        Etot += ene;
        break;
      case C_EWALD:
      case C_PME: // Elec must be enabled, vdw may not be
        ene = EW_->CalcEnergy(frm.Frm(), Imask_, ene2);
        Energy_[ELEC]->Add(frameNum, &ene);
        if (Energy_[VDW] != 0) Energy_[VDW]->Add(frameNum, &ene2);
        Etot += (ene + ene2);
        break;
      case C_KE:
        if (frm.Frm().HasVelocity()) {
          if (frm.Frm().HasForce())
            ene = ENE_.E_Kinetic_VV(frm.Frm(), *currentParm_, Imask_, dt_);
          else
            ene = ENE_.E_Kinetic(frm.Frm(), *currentParm_, Imask_);
          Energy_[KE]->Add(frameNum, &ene);
        }
        break;
    }
  }
  if (Energy_[TOTAL] != 0)
    Energy_[TOTAL]->Add(frameNum, &Etot);
  etime_.Stop();
  return Action::OK;
}

void Action_Energy::Print() {
  mprintf("Timing for energy: '%s' ('%s')\n", setname_.c_str(), Mask1_.MaskString());
  etime_.WriteTiming(0, " Total:");
  ENE_.PrintTiming(etime_.Total());
  if (elecType_ == EWALD || elecType_ == PME)
    EW_->Timing(etime_.Total());
}
