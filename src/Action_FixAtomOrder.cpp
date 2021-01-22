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
  mprintf("\t[pdborder [hetatm <mask>] (EXPERIMENTAL)]\n");
  mprintf("  Fix atom ordering so that all atoms in molecules are sequential.\n"
          "  If 'pdborder' is specified, attempt to organize atoms by PDB\n"
          "  information (i.e. Chain ID, original residue numbering, and\n"
          "  insertion codes). Atoms optionally specified by 'hetatm <mask>'\n"
          "  will be placed after all other atoms. **Note that the 'pdborder'\n"
          "  keyword is still experimental, and requires that the Topology have\n"
          "  PDB-type information present.**\n");
  mprintf("%s", ActionTopWriter::Options());
}

// Action_FixAtomOrder::init()
Action::RetType Action_FixAtomOrder::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  debug_ = debugIn;
  topWriter_.InitTopWriter(actionArgs, "re-ordered", debug_);
  mode_ = FIX_MOLECULES;
  if (actionArgs.hasKey("pdborder")) {
    mode_ = PDB_ORDER;
    std::string mexp = actionArgs.GetStringKey("hetatm");
    if (!mexp.empty()) {
      if (hetatm_.SetMaskString(mexp)) return Action::ERR;
    }
  }

  if (mode_ == FIX_MOLECULES)
    mprintf("    FIXATOMORDER: Re-ordering atoms so molecules are contiguous.\n");
  else {
    mprintf("    FIXATOMORDER: Re-ordering atoms according to PDB info.\n");
    if (hetatm_.MaskStringSet())
      mprintf("\tAtoms selected by %s will be considered HETATM.\n", hetatm_.MaskString());
  }
  topWriter_.PrintOptions();

  return Action::OK;
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
  topWriter_.WriteTops( setup.Top() );
 
  return Action::MODIFY_TOPOLOGY;
}

/** Try to make the order match original PDB info. */
Action::RetType Action_FixAtomOrder::PdbOrder(ActionSetup& setup) {
  if (hetatm_.MaskStringSet()) {
    if (setup.Top().SetupCharMask( hetatm_ )) return Action::ERR;
    hetatm_.MaskInfo();
  } else {
    hetatm_ = CharMask( setup.Top().Natom() );
  }
  // Create array with PDB info.
  std::vector<AtomTopType> atoms;
  for (int idx = 0; idx != setup.Top().Natom(); idx++)
  {
    Residue const& res = setup.Top().Res( setup.Top()[idx].ResNum() );
    AtomTopType::PdbType pt;
    if (hetatm_.AtomInCharMask(idx))
      pt = AtomTopType::HETATM;
    else
      pt = AtomTopType::ATOM;
    atoms.push_back( AtomTopType(pt, idx, res.OriginalResNum(),
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
  atomMap_.clear();
  atomMap_.reserve( setup.Top().Natom() );
  for (Topology::mol_iterator mol = setup.Top().MolStart();
                              mol != setup.Top().MolEnd(); ++mol)
  {
    for (Unit::const_iterator seg = mol->MolUnit().segBegin();
                              seg != mol->MolUnit().segEnd(); ++seg)
    {
      for (int idx = seg->Begin(); idx != seg->End(); ++idx)
        atomMap_.push_back( idx );
    }
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
