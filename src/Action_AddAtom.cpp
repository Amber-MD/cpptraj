#include "Action_AddAtom.h"
#include "CpptrajStdio.h"

/** CONSTRUCTOR */
Action_AddAtom::Action_AddAtom() :
  newParm_(0)
{}

/** DESTRUCTOR */
Action_AddAtom::~Action_AddAtom() {
  if (newParm_ != 0) delete newParm_;
}

// Action_AddAtom::Help()
void Action_AddAtom::Help() const {
  mprintf("\taname <name> [elt <element>] [rname <res name>]\n"
          "\t[xyz <X> <Y> <Z>] [mass <mass>] [charge <charge>]\n");
  mprintf("%s", ActionTopWriter::Keywords());
  mprintf("  Add an atom to current topology/coordinates.\n");
}

// Action_AddAtom::Init()
Action::RetType Action_AddAtom::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  // Get output stripped parm filename
  topWriter_.InitTopWriter(actionArgs, "modified", debugIn);

  std::string aname = actionArgs.GetStringKey("aname");
  if (aname.empty()) {
    mprinterr("Error: Must specify atom name with 'aname'.\n");
    return Action::ERR;
  }
  NameType atomName( aname );

  std::string elt = actionArgs.GetStringKey("elt");
  if (elt.empty())
    elt.assign("H");
  if (elt.size() > 2) {
    mprinterr("Error: Element name '%s' is too big; should be 2 characters max.\n", elt.c_str());
    return Action::ERR;
  }

  std::string rname = actionArgs.GetStringKey("rname");
  if (rname.empty())
    rname.assign("TMP");
  residueName_ = NameType( rname );

  bool has_mass = false;
  double mass = 0;
  if (actionArgs.Contains("mass")) {
    has_mass = true;
    mass = actionArgs.getKeyDouble("mass", 0);
  }

  bool has_charge = false;
  double charge = 0;
  if (actionArgs.Contains("charge")) {
    has_charge = true;
    charge = actionArgs.getKeyDouble("charge", 0);
  }

  // ----- Last arg to process -----
  if (actionArgs.Contains("xyz")) {
    ArgList xyzargs = actionArgs.GetNstringKey("xyz", 3);
    if (xyzargs.Nargs() != 3) {
      mprinterr("Error: Expected 3 arguments after 'xyz', got '%i'\n", xyzargs.Nargs());
      return Action::ERR;
    }
    xyz_[0] = xyzargs.getNextDouble(0);
    xyz_[1] = xyzargs.getNextDouble(0);
    xyz_[2] = xyzargs.getNextDouble(0);
  } else
    xyz_ = Vec3(0.0);

  // ----- No more args after here -----
  newAtom_ = Atom(atomName, elt.c_str());
  if (has_mass) newAtom_.SetMass( mass );
  if (has_charge) newAtom_.SetCharge( charge );

  mprintf("    ADDATOM: Adding atom named '%s', element %s, residue name '%s'\n",
          *atomName, elt.c_str(), *residueName_);
  if (has_mass) mprintf("\tAtom mass = %g\n", mass);
  if (has_charge) mprintf("\tAtom charge = %g\n", charge);
  mprintf("\tAtom will be placed at XYZ= %g %g %g\n", xyz_[0], xyz_[1], xyz_[2]);
  topWriter_.PrintOptions();

  return Action::OK;
}

// Action_AddAtom::Setup()
Action::RetType Action_AddAtom::Setup(ActionSetup& setup)
{
  // Copy existing topology
  if (newParm_ != 0) delete newParm_;
  newParm_ = new Topology( setup.Top() );
  // Set up the new residue
  Residue const& lastRes = setup.Top().Res( setup.Top().Nres()-1 );
  int newResNum = lastRes.OriginalResNum() + 1;
  Residue newRes( residueName_, newResNum, ' ', "" );
  // Add to new topology
  newParm_->AddTopAtom( newAtom_, newRes );
  // If topology has PDB info, set that as well TODO handle in Topology?
  if (!newParm_->AtomAltLoc().empty()) newParm_->AddAtomAltLoc(' ');
  if (!newParm_->Occupancy().empty()) newParm_->AddOccupancy(1.0);
  if (!newParm_->Bfactor().empty()) newParm_->AddBfactor(0.0);
  if (!newParm_->PdbSerialNum().empty())
    newParm_->AddPdbSerialNum( newParm_->PdbSerialNum().back() + 1 );
  // If topology has extra Amber info, add placeholder values
  if (!newParm_->TreeChainClassification().empty()) newParm_->AddTreeChainClassification("BLA");
  if (!newParm_->JoinArray().empty()) newParm_->AddJoinArray(0);
  if (!newParm_->RotateArray().empty()) newParm_->AddRotateArray(0);
  // Final setup
  newParm_->CommonSetup();
  mprintf("\tAdded '%s'\n", newParm_->AtomMaskName( setup.Top().Natom() ).c_str());

  setup.SetTopology( newParm_ );
  // Remove box information if asked
  if (topWriter_.ModifyActionState(setup, newParm_))
    return Action::ERR;
  // Allocate space for new frame
  newFrame_.SetupFrameV(setup.Top().Atoms(), setup.CoordInfo());

  // If prefix given then output stripped topology
  topWriter_.WriteTops( setup.Top() );

  return Action::MODIFY_TOPOLOGY;
}

// Action_AddAtom::DoAction()
Action::RetType Action_AddAtom::DoAction(int frameNum, ActionFrame& frm)
{

  newFrame_.CopyFrom( frm.Frm(), 0, newParm_->Natom()-1 );
  int idx = (newParm_->Natom()-1) * 3;
  newFrame_[idx  ] = xyz_[0];
  newFrame_[idx+1] = xyz_[1];
  newFrame_[idx+2] = xyz_[2];

  // Set frame
  frm.SetFrame( &newFrame_ );

  return Action::MODIFY_COORDS;
}
