// Action_FixAtomOrder 
#include "Action_FixAtomOrder.h"
#include "CpptrajStdio.h"
#include "ParmFile.h"

// CONSTRUCTOR
Action_FixAtomOrder::Action_FixAtomOrder() : debug_(0), newParm_(0) {} 

Action_FixAtomOrder::~Action_FixAtomOrder() {
  if (newParm_!=0) delete newParm_;
}

void Action_FixAtomOrder::Help() const {
  mprintf("\t[outprefix <name>]\n"
          "  Fix atom ordering so that all atoms in molecules are sequential.\n");
}

// Action_FixAtomOrder::init()
Action::RetType Action_FixAtomOrder::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  debug_ = debugIn;
  prefix_ = actionArgs.GetStringKey("outprefix");
  mprintf("    FIXATOMORDER: Will attempt to fix atom ordering when atom numbering\n"
          "                  in molecules is non-sequential.\n");
  if (!prefix_.empty())
    mprintf("\tRe-ordered topology will be output with prefix %s\n",prefix_.c_str());
  return Action::OK;
}

void Action_FixAtomOrder::VisitAtom(int atomnum, int mol, Topology const& Parm) {
  // Return if this atom already has a molecule number
  if (molNums_[atomnum]!=-1) return;
  // Mark this atom as visited
  molNums_[atomnum] = mol;
  // Visit each atom bonded to this atom
  for (Atom::bond_iterator bondedatom = Parm[atomnum].bondbegin();
                           bondedatom != Parm[atomnum].bondend(); bondedatom++)
    VisitAtom(*bondedatom, mol, Parm);
}


// Action_FixAtomOrder::setup()
Action::RetType Action_FixAtomOrder::Setup(ActionSetup& setup) {
  // If topology already has molecule info assume no need to reorder.
  if (setup.Top().Nmol() > 0) {
    mprintf("Warning: %s already has molecule information. No reordering will occur.\n"
            "Warning: This indicates that there is no need to fix atom ordering in this topology.\n",
            setup.Top().c_str());
    return Action::SKIP;
  }
  molNums_.resize( setup.Top().Natom(), -1 );
  // Perform recursive search along bonds of each atom.
  int Nmol = 0;
  for (int atomnum = 0; atomnum < setup.Top().Natom(); ++atomnum)
  {
    if (molNums_[atomnum] == -1) {
      VisitAtom( atomnum, Nmol, setup.Top() );
      ++Nmol;
    }
  }
  mprintf("\tDetected %i molecules.\n", Nmol);
  if (Nmol < 1) {
    mprinterr("Error: No molecules detected in %s\n", setup.Top().c_str());
    return Action::ERR;
  }
  if (debug_ > 0) {
    for (MapType::const_iterator mnum = molNums_.begin(); mnum != molNums_.end(); ++mnum)
      mprintf("\t\tAtom %li assigned to molecule %i\n", mnum - molNums_.begin() + 1, *mnum + 1);
  }
  // Figure out which atoms should go in which molecules 
  std::vector<MapType> molecules(Nmol);
  for (int atomnum = 0; atomnum < setup.Top().Natom(); ++atomnum)
    molecules[molNums_[atomnum]].push_back( atomnum );
  atomMap_.clear();
  atomMap_.reserve( setup.Top().Natom() );
  // Place all atoms in molecule 0 first, molecule 1 next and so on
  for (std::vector<MapType>::const_iterator mol = molecules.begin();
                                            mol != molecules.end(); ++mol)
    for (MapType::const_iterator atom = mol->begin(); atom != mol->end(); ++atom)
      atomMap_.push_back( *atom );
  if (debug_ > 0) {
    mprintf("\tNew atom mapping:\n");
    for (MapType::const_iterator atom = atomMap_.begin();
                                 atom != atomMap_.end(); ++atom)
      mprintf("\t\tNew atom %8li => old atom %8i\n", atom - atomMap_.begin() + 1, *atom + 1);
  }
  // Create new topology based on map
  if (newParm_ != 0) delete newParm_;
  newParm_ = setup.Top().ModifyByMap( atomMap_ );
  if (newParm_ == 0) {
    mprinterr("Error: Could not create re-ordered topology.\n");
    return Action::ERR;
  }
  newParm_->Brief("Re-ordered parm:");
  setup.SetTopology( newParm_ );
  // Allocate space for new frame
  newFrame_.SetupFrameV( setup.Top().Atoms(), setup.CoordInfo() );

  // If prefix given then output stripped parm
  if (!prefix_.empty()) {
    ParmFile pfile;
    if ( pfile.WritePrefixTopology( setup.Top(), prefix_, ParmFile::AMBERPARM, debug_ ) )
      mprinterr("Error: Could not write out reordered parm file.\n");
  }

  return Action::MODIFY_TOPOLOGY;
}

// Action_FixAtomOrder::action()
Action::RetType Action_FixAtomOrder::DoAction(int frameNum, ActionFrame& frm) {
  // Reorder atoms in the frame
  newFrame_.SetCoordinatesByMap( frm.Frm(), atomMap_ );
  frm.SetFrame( &newFrame_ );
  return Action::MODIFY_COORDS;
}
