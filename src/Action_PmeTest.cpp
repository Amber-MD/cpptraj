#include "Action_PmeTest.h"
#include "CpptrajStdio.h"

Action_PmeTest::Action_PmeTest() {
  SetHidden(true);
}

// Action_PmeTest::Help()
void Action_PmeTest::Help() const {
  mprintf("\t[ %s\n", EwaldOptions::KeywordsPME());
  mprintf("\t  %s\n", EwaldOptions::KeywordsCommon1());
  mprintf("\t  %s ]\n", EwaldOptions::KeywordsCommon2());
  mprintf("\t{{orginal|new|newljpme} | ewald {original|new}} [out <file>] [<mask>] [<setname>]\n");
}

// Action_PmeTest::Init()
Action::RetType Action_PmeTest::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  debug_ = debugIn;
  EwaldOptions::OptType opt = EwaldOptions::PME;
  if (actionArgs.hasKey("ewald"))
    opt = EwaldOptions::REG_EWALD;
  if (ewaldOpts_.GetOptions(opt, actionArgs, "pmetest")) {
    mprinterr("Error: Could not get Ewald options.\n");
    return Action::ERR;
  }
  if (opt == EwaldOptions::PME) {
    if (actionArgs.hasKey("original"))
      method_ = 0;
    else if (actionArgs.hasKey("new"))
      method_ = 1;
    else if (actionArgs.hasKey("newljpme"))
      method_ = 2;
    else {
      mprinterr("Error: Specify original, new, or newljpme.\n");
      return Action::ERR;
    }
  } else {
    if (actionArgs.hasKey("original"))
      method_ = 3;
    else if (actionArgs.hasKey("new"))
      method_ = 4;
    else {
      mprinterr("Error: Specify original or new.\n");
      return Action::ERR;
    }
  }

  DataFile* outfile = init.DFL().AddDataFile( actionArgs.GetStringKey("out"), actionArgs );

  if (Mask1_.SetMaskString( actionArgs.GetMaskNext() )) return Action::ERR;
  std::string setname_ = actionArgs.GetStringNext();
  if (setname_.empty())
    setname_ = init.DSL().GenerateDefaultName("PME");
  ele_ = init.DSL().AddSet( DataSet::DOUBLE, MetaData(setname_, "ELE") );
  vdw_ = init.DSL().AddSet( DataSet::DOUBLE, MetaData(setname_, "VDW") );
  if (ele_ == 0 || vdw_ == 0) {
    mprinterr("Error: Could not allocate data sets.\n");
    return Action::ERR;
  }
  if (outfile != 0) {
    outfile->AddDataSet(ele_);
    outfile->AddDataSet(vdw_);
  }

  mprintf("    ENERGY: Calculating energy for atoms in mask '%s'\n", Mask1_.MaskString());
  if (method_ == 0)
    mprintf("\tOriginal PME method.\n");
  else if (method_ == 1)
    mprintf("\tNew PME method.\n");
  else if (method_ == 2)
    mprintf("\tNew LJPME method.\n");
  else if (method_ == 3)
    mprintf("\tOriginal Ewald method.\n");
  else if (method_ == 4)
    mprintf("\tNew Ewald method.\n");

  ewaldOpts_.PrintOptions();

  return Action::OK;
}

// Action_PmeTest::Setup()
Action::RetType Action_PmeTest::Setup(ActionSetup& setup)
{
  using namespace Cpptraj::Energy;
  // Set up mask
  if (setup.Top().SetupIntegerMask(Mask1_)) return Action::ERR;
  if (Mask1_.None()) {
    mprintf("Warning: Mask '%s' selects no atoms.\n", Mask1_.MaskString());
    return Action::SKIP;
  }
  Mask1_.MaskInfo();

  if (method_ == 0) {
    if (PME0_.Init(setup.CoordInfo().TrajBox(), ewaldOpts_, debug_))
      return Action::ERR;
    PME0_.Setup( setup.Top(), Mask1_ );
  } else if (method_ == 3) {
    if (EWALD0_.Init(setup.CoordInfo().TrajBox(), ewaldOpts_, debug_))
      return Action::ERR;
    EWALD0_.Setup( setup.Top(), Mask1_ );
  } else {
    Ecalc_Nonbond::CalcType opt = Ecalc_Nonbond::UNSPECIFIED;
    if (method_ == 1)
      opt = Ecalc_Nonbond::PME;
    else if (method_ == 2)
      opt = Ecalc_Nonbond::LJPME;
    else if (method_ == 4)
      opt = Ecalc_Nonbond::REGULAR_EWALD;

    if (NB_.InitNonbondCalc(opt, false, setup.CoordInfo().TrajBox(), ewaldOpts_, debug_))
      return Action::ERR;
    if (NB_.SetupNonbondCalc( setup.Top(), Mask1_ ))
      return Action::ERR;
  } 
  return Action::OK;
}

// Action_PmeTest::DoAction()
Action::RetType Action_PmeTest::DoAction(int frameNum, ActionFrame& frm)
{
  t_nb_.Start();
  double ene, ene2;
  int err = 1;
  if (method_ == 0) {
    err = PME0_.CalcNonbondEnergy(frm.Frm(), Mask1_, ene, ene2);
  } else if (method_ == 3) {
    err = EWALD0_.CalcNonbondEnergy(frm.Frm(), Mask1_, ene, ene2);
  } else {
    err = NB_.NonbondEnergy(frm.Frm(), Mask1_, ene, ene2);
  }
  if (err != 0) return Action::ERR;
  ele_->Add(frameNum, &ene);
  vdw_->Add(frameNum, &ene2);
  t_nb_.Stop();
  return Action::OK;
}

void Action_PmeTest::Print() {
  if (method_ == 0)
    PME0_.Timing(t_nb_.Total());
  else if (method_ == 3)
    EWALD0_.Timing(t_nb_.Total());
  else 
    NB_.PrintTiming(t_nb_.Total());
}
