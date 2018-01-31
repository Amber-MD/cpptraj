#include "Action_Energy.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_Energy::Action_Energy() : currentParm_(0), debug_(0)
{
  std::fill(mlimits_, mlimits_+3, 0);
}


void Action_Energy::Help() const {
  mprintf("\t[<name>] [<mask1>] [out <filename>]\n"
          "\t[bond] [angle] [dihedral] [nb14] {[nonbond] | [elec] [vdw]}\n"
          "\t[ etype {simple | directsum [npoints <N>] |\n"
          "\t         ewald [cut <cutoff>] [dsumtol <dtol>] [rsumtol <rtol>]\n"
          "\t               [ewcoeff <coeff>] [maxexp <max>] [skinnb <skinnb>]\n"
          "\t               [mlimits <X>,<Y>,<Z>]} ]\n"
          "  Calculate energy for atoms in mask.\n");
}

/// DataSet aspects
static const char* Estring[] = {"bond", "angle", "dih", "vdw14", "elec14", "vdw", "elec", "total"};

/// Calculation types (CalcType)
static const char* Cstring[] = {"Bond", "Angle", "Torsion", "1-4 Nonbond", "Nonbond",
                                "Electrostatics", "van der Waals", "Electrostatics (Direct Sum)",
                                "Electrostatics (Ewald)", "Electrostatics (PME)" };

// Action_Energy::AddSet()
int Action_Energy::AddSet(Etype typeIn, DataSetList& DslIn, DataFile* outfile,
                          std::string const& setname)
{
  Energy_[typeIn] = DslIn.AddSet(DataSet::DOUBLE, MetaData(setname, Estring[typeIn]));
  if (Energy_[typeIn] == 0) return 1;
  if (outfile != 0) outfile->AddDataSet( Energy_[typeIn] );
  return 0;
}

// Action_Energy::Init()
Action::RetType Action_Energy::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  debug_ = debugIn;
  ENE_.SetDebug( debug_ );
  // Get keywords
  DataFile* outfile = init.DFL().AddDataFile( actionArgs.GetStringKey("out"), actionArgs );

  // Which terms will be calculated?
  bool calc_vdw  = actionArgs.hasKey("vdw" );
  bool calc_elec = actionArgs.hasKey("elec");
  bool calc_nb   = actionArgs.hasKey("nonbond");
  if (calc_vdw && calc_elec)
    calc_nb = true;
  if (calc_nb) {
    calc_vdw = false;
    calc_elec = false;
  }
  // Electrostatics type. If specified always split the Elec/VDW calc.
  etype_ = SIMPLE;
  std::string etypearg = actionArgs.GetStringKey("etype");
  if (!etypearg.empty()) {
    if (calc_nb) {
      calc_nb = false;
      calc_vdw = true;
    }
    if (etypearg == "directsum") {
      // Direct sum method
      etype_ = DIRECTSUM;
      calc_elec = true;
      npoints_ = actionArgs.getKeyInt("npoints", 0);
    } else if (etypearg == "ewald") {
      // Ewald method
      etype_ = EW;
      calc_elec = true;
      cutoff_ = actionArgs.getKeyDouble("cut", 8.0);
      dsumtol_ = actionArgs.getKeyDouble("dsumtol", 1E-5);
      rsumtol_ = actionArgs.getKeyDouble("rsumtol", 5E-5);
      ewcoeff_ = actionArgs.getKeyDouble("ewcoeff", 0.0);
      maxexp_ = actionArgs.getKeyDouble("maxexp", 0.0);
      skinnb_ = actionArgs.getKeyDouble("skinnb", 2.0);
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
    } else if (etypearg == "pme") {
      // Ewald method
      etype_ = PME;
      calc_elec = true;
      cutoff_ = actionArgs.getKeyDouble("cut", 8.0);
      dsumtol_ = actionArgs.getKeyDouble("dsumtol", 1E-5);
      rsumtol_ = actionArgs.getKeyDouble("rsumtol", 5E-5);
      ewcoeff_ = actionArgs.getKeyDouble("ewcoeff", 0.0);
      skinnb_ = actionArgs.getKeyDouble("skinnb", 2.0);
    } else if (etypearg == "simple") {
      // Simple method
      etype_ = SIMPLE;
      if (!calc_nb && !calc_elec) calc_elec = true;
    } else {
      mprinterr("Error: Unrecognized option for 'etype': %s\n", etypearg.c_str());
      return Action::ERR;
    }
  }
  // Set up calculations
  Ecalcs_.clear();
  if (actionArgs.hasKey("bond"))     Ecalcs_.push_back(BND);
  if (actionArgs.hasKey("angle"))    Ecalcs_.push_back(ANG);
  if (actionArgs.hasKey("dihedral")) Ecalcs_.push_back(DIH);
  if (actionArgs.hasKey("nb14"))     Ecalcs_.push_back(N14);
  if (calc_nb)                       Ecalcs_.push_back(NBD);
  if (calc_vdw)                      Ecalcs_.push_back(LJ);
  if (calc_elec) {
    switch (etype_) {
      case SIMPLE:    Ecalcs_.push_back(COULOMB); break;
      case DIRECTSUM: Ecalcs_.push_back(DIRECT); break;
      case EW:        Ecalcs_.push_back(EWALD); break;
      case PME:       Ecalcs_.push_back(PMEWALD); break;
    }
  }
  // If nothing is selected, select all.
  if (Ecalcs_.empty()) {
    for (int c = 0; c <= (int)NBD; c++)
      Ecalcs_.push_back( (CalcType)c );
  }

  // Get Masks
  Mask1_.SetMaskString( actionArgs.GetMaskNext() );

  // DataSet
  std::string setname = actionArgs.GetStringNext();
  if (setname.empty())
    setname = init.DSL().GenerateDefaultName("ENE");
  Energy_.clear();
  Energy_.resize( (int)TOTAL + 1, 0 );
  for (calc_it calc = Ecalcs_.begin(); calc != Ecalcs_.end(); ++calc)
  {
    switch (*calc) {
      case BND: if (AddSet(BOND, init.DSL(), outfile, setname)) return Action::ERR; break;
      case ANG: if (AddSet(ANGLE, init.DSL(), outfile, setname)) return Action::ERR; break;
      case DIH: if (AddSet(DIHEDRAL, init.DSL(), outfile, setname)) return Action::ERR; break;
      case N14:
        if (AddSet(V14, init.DSL(), outfile, setname)) return Action::ERR;
        if (AddSet(Q14, init.DSL(), outfile, setname)) return Action::ERR;
        break;
      case NBD:
        if (AddSet(VDW, init.DSL(), outfile, setname)) return Action::ERR;
        if (AddSet(ELEC, init.DSL(), outfile, setname)) return Action::ERR;
        break;
      case LJ:
        if (AddSet(VDW, init.DSL(), outfile, setname)) return Action::ERR; break;
      case COULOMB:
      case DIRECT:
      case EWALD:
      case PMEWALD:
        if (AddSet(ELEC, init.DSL(), outfile, setname)) return Action::ERR; break;
    }
  }
//  if (Ecalcs_.size() > 1) {
    if (AddSet(TOTAL, init.DSL(), outfile, setname)) return Action::ERR;
//  }
      
  mprintf("    ENERGY: Calculating energy for atoms in mask '%s'\n", Mask1_.MaskString());
  mprintf("\tCalculating terms:");
  for (calc_it calc = Ecalcs_.begin(); calc != Ecalcs_.end(); ++calc) {
    if (calc != Ecalcs_.begin()) mprintf(",");
    mprintf(" %s", Cstring[*calc]);
  }
  mprintf("\n");
  if (etype_ == DIRECTSUM) {
    if (npoints_ < 0)
      mprintf("\tDirect sum energy for up to %i unit cells in each direction will be calculated.\n",
              -npoints_);
    else
      mprintf("\tDirect sum energy for %i unit cells in each direction will be calculated.\n",
              npoints_);
  } else if (etype_ == EW) {
    mprintf("\tCalculating electrostatics with Ewald method.\n");
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
  } else if (etype_ == PME) {
    mprintf("\tCalculating electrostatics with particle mesh Ewald method.\n");
    mprintf("\tDirect space cutoff= %.4f\n", cutoff_);
    if (dsumtol_ != 0.0)
      mprintf("\tDirect sum tolerance= %g\n", dsumtol_);
    if (rsumtol_ != 0.0)
      mprintf("\tReciprocal sum tolerance= %g\n", rsumtol_);
    if (ewcoeff_ == 0.0)
      mprintf("\tWill determine Ewald coefficient from cutoff and direct sum tolerance.\n");
    else
      mprintf("\tEwald coefficient= %.4f\n", ewcoeff_);
  }
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
  for (calc_it calc = Ecalcs_.begin(); calc != Ecalcs_.end(); ++calc)
    if ((*calc == N14 || *calc == NBD) && !setup.Top().Nonbond().HasNonbond())
    {
      mprinterr("Error: Nonbonded energy calc requested but topology '%s'\n"
                "Error:   does not have non-bonded parameters.\n", setup.Top().c_str());
      return Action::ERR;
    }
  // Set up Ewald if necessary.
  if (etype_ == EW) { // TODO erfc table dx option
    if (EW_.EwaldInit(setup.CoordInfo().TrajBox(), cutoff_, dsumtol_, rsumtol_,
                      ewcoeff_, maxexp_, skinnb_, 0.0, debug_, mlimits_))
      return Action::ERR;
    EW_.EwaldSetup( setup.Top(), Imask_ );
  } else if (etype_ == PME) {
    if (EW_.PME_Init(setup.CoordInfo().TrajBox(), cutoff_, dsumtol_, rsumtol_,
                     ewcoeff_, skinnb_, 0.0, debug_))
      return Action::ERR;
    EW_.PME_Setup( setup.Top(), Imask_ );
  }
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
  for (calc_it calc = Ecalcs_.begin(); calc != Ecalcs_.end(); ++calc)
  {
    switch (*calc) {
      case BND:
        ene = ENE_.E_bond(frm.Frm(), *currentParm_, Mask1_);
        Energy_[BOND]->Add(frameNum, &ene);
        Etot += ene;
        break;
      case ANG:
        ene = ENE_.E_angle(frm.Frm(), *currentParm_, Mask1_);
        Energy_[ANGLE]->Add(frameNum, &ene);
        Etot += ene;
        break;
      case DIH:
        ene = ENE_.E_torsion(frm.Frm(), *currentParm_, Mask1_);
        Energy_[DIHEDRAL]->Add(frameNum, &ene);
        Etot += ene;
        break;
      case N14:
        ene = ENE_.E_14_Nonbond(frm.Frm(), *currentParm_, Mask1_, ene2);
        Energy_[V14]->Add(frameNum, &ene);
        Energy_[Q14]->Add(frameNum, &ene2);
        Etot += (ene + ene2);
        break;
      case NBD:
        ene = ENE_.E_Nonbond(frm.Frm(), *currentParm_, Imask_, ene2);
        Energy_[VDW]->Add(frameNum, &ene);
        Energy_[ELEC]->Add(frameNum, &ene2);
        Etot += (ene + ene2);
        break;
      case LJ:
        ene = ENE_.E_VDW(frm.Frm(), *currentParm_, Imask_);
        Energy_[VDW]->Add(frameNum, &ene);
        Etot += ene;
        break;
      case COULOMB:
        ene = ENE_.E_Elec(frm.Frm(), *currentParm_, Imask_);
        Energy_[ELEC]->Add(frameNum, &ene);
        Etot += ene;
        break;
      case DIRECT:
        if (npoints_ < 0)
          ene = Dbg_Direct(frm.Frm(), (-npoints_)+1);
        else
          ene = ENE_.E_DirectSum(frm.Frm(), *currentParm_, Imask_, npoints_);
        Energy_[ELEC]->Add(frameNum, &ene);
        Etot += ene;
        break;
      case EWALD:
        ene = EW_.CalcEnergy(frm.Frm(), Imask_);
        Energy_[ELEC]->Add(frameNum, &ene);
        Etot += ene;
        break;
      case PMEWALD:
        ene = EW_.CalcPmeEnergy(frm.Frm(), *currentParm_, Imask_);
        Energy_[ELEC]->Add(frameNum, &ene);
        Etot += ene;
        break;
    }
  }

  Energy_[TOTAL]->Add(frameNum, &Etot);
  etime_.Stop();
  return Action::OK;
}

void Action_Energy::Print() {
  mprintf("Timing for energy: '%s' ('%s')\n", Energy_[TOTAL]->legend(),
           Mask1_.MaskString());
  etime_.WriteTiming(0, " Total:");
  ENE_.PrintTiming(etime_.Total());
  if (etype_ == EW || etype_ == PME)
    EW_.Timing(etime_.Total());
}
