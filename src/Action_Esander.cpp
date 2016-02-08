#include "Action_Esander.h"
#include "CpptrajStdio.h"


#ifdef USE_SANDERLIB
// CONSTRUCTOR
Action_Esander::Action_Esander() : currentParm_(0) {}

int Action_Esander::AddSet(Energy_Sander::Etype typeIn, DataSetList& DslIn, DataFile* outfile,
                           std::string const& setname)
{
  Esets_[typeIn] = DslIn.AddSet(DataSet::DOUBLE, MetaData(setname, Energy_Sander::Easpect(typeIn)));
  if (Esets_[typeIn] == 0) return 1;
  if (outfile != 0) outfile->AddDataSet( Esets_[typeIn] );
  return 0;
}
#else
Action_Esander::Action_Esander() {}
#endif

void Action_Esander::Help() const {
# ifdef USE_SANDERLIB
  mprintf("\t[<name>] [out <filename>]\n"
          "  Calculate energy for atoms in mask using Sander energy routines.\n");
# else
  mprintf("Warning: CPPTRAJ was compiled without libsander. This Action is disabled.\n");
# endif
}
// Action_Esander::Init()
Action::RetType Action_Esander::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
# ifdef USE_SANDERLIB
  //ENE_.SetDebug( debugIn );
  // Get keywords
  DataFile* outfile = init.DFL().AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  ReferenceFrame REF = init.DSL().GetReferenceFrame( actionArgs );
  if (REF.error()) return Action::ERR;
  if (!REF.empty()) {
    refFrame_ = REF.Coord();
    currentParm_ = REF.ParmPtr();
  }
  if (SANDER_.SetInput( actionArgs )) return Action::ERR;
  // DataSet
  std::string setname = actionArgs.GetStringNext();
  if (setname.empty())
    setname = init.DSL().GenerateDefaultName("ENE");
  Esets_.clear();
  int Nsets = (int)Energy_Sander::N_ENERGYTYPES; // TODO only add sets that will be calcd
  Esets_.resize( Nsets, 0 );
  for (int ie = 0; ie != Nsets; ie++)
    if (AddSet((Energy_Sander::Etype)ie, init.DSL(), outfile, setname)) return Action::ERR;
      
  mprintf("    ESANDER: Calculating energy using Sander.\n");
  mprintf("\tReference for initialization");
  if (!REF.empty())
    mprintf(" is '%s'\n", REF.refName());
  else
    mprintf(" will be first frame.\n");
  mprintf("\tCalculating terms:");
  for (Earray::const_iterator it = Esets_.begin(); it != Esets_.end(); ++it)
    mprintf(" %s", (*it)->legend());
  mprintf("\n");
  return Action::OK;
# else
  mprinterr("Error: CPPTRAJ was compiled without libsander. This Action is disabled.\n");
  return Action::ERR;
# endif
}

// Action_Esander::Setup()
/** Set angle up for this parmtop. Get masks etc.
  */
Action::RetType Action_Esander::Setup(ActionSetup& setup) {
# ifdef USE_SANDERLIB
  if (currentParm_ != 0 && currentParm_->Pindex() != setup.Top().Pindex())
  {
    mprintf("Warning: Current topology is %i:%s but reference is %i:%s. Skipping.\n",
            setup.Top().Pindex(), setup.Top().c_str(),
            currentParm_->Pindex(), currentParm_->c_str());
    return Action::SKIP;
  }
  // Check for LJ terms
  if (!setup.Top().Nonbond().HasNonbond())
  {
    mprinterr("Error: Topology '%s' does not have non-bonded parameters.\n", setup.Top().c_str());
    return Action::ERR;
  }
  // If reference specified, init now. Otherwise using first frame.
  if (currentParm_ != 0 ) {
    if ( SANDER_.Initialize( *currentParm_, refFrame_ ) ) return Action::ERR;
  } else
    currentParm_ = setup.TopAddress();
  return Action::OK;
# else
  return Action::ERR;
# endif
}

void Action_Esander::AddEne(Energy_Sander::Etype t, int fn) {
  Esets_[t]->Add(fn, SANDER_.Eptr(t));
}

// Action_Esander::DoAction()
Action::RetType Action_Esander::DoAction(int frameNum, ActionFrame& frm) {
# ifdef USE_SANDERLIB
  if (refFrame_.empty()) {
    refFrame_ = frm.Frm();
    if ( SANDER_.Initialize( *currentParm_, refFrame_ ) ) return Action::ERR;
  }
  // FIXME: Passing in ModifyFrm() to give CalcEnergy access to non-const pointers
  SANDER_.CalcEnergy( frm.ModifyFrm() );

  AddEne(Energy_Sander::BOND, frameNum);
  AddEne(Energy_Sander::ANGLE, frameNum);
  AddEne(Energy_Sander::DIHEDRAL, frameNum);
  AddEne(Energy_Sander::VDW14, frameNum);
  AddEne(Energy_Sander::ELEC14, frameNum);
  AddEne(Energy_Sander::VDW, frameNum);
  AddEne(Energy_Sander::ELEC, frameNum);
  AddEne(Energy_Sander::TOTAL, frameNum);

  return Action::OK;
# else
  return Action::ERR;
# endif
}

void Action_Esander::Print() {}
