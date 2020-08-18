// Action_FixAtomOrder 
#include "Action_FixAtomOrder.h"
#include "CpptrajStdio.h"
#include "AtomTopType.h"
#include <algorithm> // std::sort

/** CONSTRUCTOR */
Action_FixAtomOrder::Action_FixAtomOrder() :
  debug_(0),
  newParm_(0)
{} 

/** DESTRUCTOR */
Action_FixAtomOrder::~Action_FixAtomOrder() {
  if (newParm_ != 0) delete newParm_;
}

// void Action_FixAtomOrder::Help()
void Action_FixAtomOrder::Help() const {
  mprintf("%s", ActionTopWriter::Keywords());
  mprintf("  Fix atom ordering so that all atoms in molecules are sequential.\n");
  mprintf("%s", ActionTopWriter::Options());
}

// Action_FixAtomOrder::init()
Action::RetType Action_FixAtomOrder::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  debug_ = debugIn;
  topWriter_.InitTopWriter(actionArgs, "re-ordered", debug_);
  mode_ = FIX_MOLECULES;
  if (actionArgs.hasKey("pdborder"))
    mode_ = PDB_ORDER;

  if (mode_ == FIX_MOLECULES)
    mprintf("    FIXATOMORDER: Will attempt to fix atom ordering when atom numbering\n"
            "                  in molecules is non-sequential.\n");
  else
    mprintf("    FIXATOMORDER: Will attempt to re-order according to PDB info.\n");
  topWriter_.PrintOptions();

  return Action::OK;
}

/** Mark atom as visited, visit all bonded atoms and mark them as well. */
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
  Action::RetType ret;
  if (mode_ == FIX_MOLECULES)
    ret = FixMolecules(setup);
  else if (mode_ == PDB_ORDER)
    ret = PdbOrder(setup);
  else
    return Action::ERR;

  // If not OK, means we should bail now.
  if (ret != Action::OK)
    return ret;

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
  topWriter_.WriteTops( setup.Top() );
 
  return Action::MODIFY_TOPOLOGY;
}

/** Try to make the order match original PDB info. */
Action::RetType Action_FixAtomOrder::PdbOrder(ActionSetup& setup) {
  // Create array with PDB info.
  std::vector<AtomTopType> atoms;
  for (int idx = 0; idx != setup.Top().Natom(); idx++)
  {
    Residue const& res = setup.Top().Res( setup.Top()[idx].ResNum() );
    atoms.push_back( AtomTopType(idx, res.OriginalResNum(),
                                 res.Icode(), res.ChainId()) );
  }
  // Sort by PDB info
  std::sort(atoms.begin(), atoms.end());
  // Put indices in atomMap_
  atomMap_.clear();
  atomMap_.reserve( setup.Top().Natom() );
  for (std::vector<AtomTopType>::const_iterator it = atoms.begin();
                                                it != atoms.end(); ++it)
    atomMap_.push_back( it->AtomIdx() );

  return Action::OK;
}

/** Fix molecules that are not contiguous. */
Action::RetType Action_FixAtomOrder::FixMolecules(ActionSetup& setup) {
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

  return Action::OK;
}

// Action_FixAtomOrder::action()
Action::RetType Action_FixAtomOrder::DoAction(int frameNum, ActionFrame& frm) {
  // Reorder atoms in the frame
  newFrame_.SetCoordinatesByMap( frm.Frm(), atomMap_ );
  frm.SetFrame( &newFrame_ );
  return Action::MODIFY_COORDS;
}
