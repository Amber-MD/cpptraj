#include "Action_Esander.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_Esander::Action_Esander() : currentParm_(0) {}


void Action_Esander::Help() const {
# ifdef USE_SANDERLIB
  mprintf("\t[<name>] [out <filename>]\n"
          "  Calculate energy for atoms in mask using Sander energy routines.\n");
# else
  mprintf("Warning: CPPTRAJ was compiled without libsander. This Action is disabled.\n");
# endif
}

/// DataSet aspects
static const char* Estring[] = {"bond", "angle", "dih", "vdw14", "elec14", "vdw", "elec", "total"};

int Action_Esander::AddSet(Etype typeIn, DataSetList& DslIn, DataFile* outfile,
                           std::string const& setname)
{
  Esets_[typeIn] = DslIn.AddSet(DataSet::DOUBLE, MetaData(setname, Estring[typeIn]));
  if (Esets_[typeIn] == 0) return 1;
  if (outfile != 0) outfile->AddDataSet( Esets_[typeIn] );
  return 0;
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
  // DataSet
  std::string setname = actionArgs.GetStringNext();
  if (setname.empty())
    setname = init.DSL().GenerateDefaultName("ENE");
  Esets_.clear();
  int Nsets = (int)TOTAL + 1;
  Esets_.resize( Nsets, 0 );
  for (int ie = 0; ie != Nsets; ie++)
    if (AddSet((Etype)ie, init.DSL(), outfile, setname)) return Action::ERR;
      
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
    mprinterr("Error: Nonbonded energy calc requested but topology '%s'\n"
              "Error:   does not have non-bonded parameters.\n", setup.Top().c_str());
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

// Action_Esander::DoAction()
Action::RetType Action_Esander::DoAction(int frameNum, ActionFrame& frm) {
# ifdef USE_SANDERLIB
  if (refFrame_.empty()) {
    refFrame_ = frm.Frm();
    if ( SANDER_.Initialize( *currentParm_, refFrame_ ) ) return Action::ERR;
  }
  // FIXME: Passing in ModifyFrm() to give CalcEnergy access to non-const pointers
  SANDER_.CalcEnergy( frm.ModifyFrm() );

  Esets_[BOND]->Add(frameNum, SANDER_.EbondPtr());
  Esets_[ANGLE]->Add(frameNum, SANDER_.EanglePtr());
  Esets_[DIHEDRAL]->Add(frameNum, SANDER_.EdihedralPtr());
  Esets_[V14]->Add(frameNum, SANDER_.Evdw14Ptr());
  Esets_[Q14]->Add(frameNum, SANDER_.Eelec14Ptr());
  Esets_[VDW]->Add(frameNum, SANDER_.EvdwPtr());
  Esets_[ELEC]->Add(frameNum, SANDER_.EelecPtr());
  double totalE = SANDER_.Etotal();
  Esets_[TOTAL]->Add(frameNum, &totalE);

  //Energy_[TOTAL]->Add(frameNum, &Etot);
  //sanderout_.Printf("%8i %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f\n", frameNum+1,
  //        SANDER_.Ebond(), SANDER_.Eangle(), SANDER_.Edihedral(),
  //        SANDER_.Evdw14(), SANDER_.Eelec14(), SANDER_.Evdw(), SANDER_.Eelec(),
  //        SANDER_.Ebond() + SANDER_.Eangle() + SANDER_.Edihedral() +
  //        SANDER_.Evdw14() + SANDER_.Eelec14() + SANDER_.Evdw() + SANDER_.Eelec());
  return Action::OK;
# else
  return Action::ERR;
# endif
}

void Action_Esander::Print() {}
