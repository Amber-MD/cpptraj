#include "Exec_Change.h"
#include "CpptrajStdio.h"

// Exec_Change::Help()
void Exec_Change::Help() const
{
  mprintf("\t[ {%s |\n"
          "\t   crdset <COORDS set> ]\n"
          "\t{ resname from <mask> to <value> |\n"
          "\t  chainid of <mask> to <value> |\n"
          "\t  atomname from <mask> to <value> |\n"
          "\t  addbond <mask1> <mask2> [req <length> <rk> <force constant>] }\n"
          "  Change specified parts of topology or topology of a COORDS data set.\n",
          DataSetList::TopArgs);
}

// Exec_Change::Execute()
Exec::RetType Exec_Change::Execute(CpptrajState& State, ArgList& argIn)
{
  // Change type
  enum ChangeType { UNKNOWN = 0, RESNAME, CHAINID, ATOMNAME, ADDBOND };
  ChangeType type = UNKNOWN;
  if (argIn.hasKey("resname"))
    type = RESNAME;
  else if (argIn.hasKey("chainid"))
    type = CHAINID;
  else if (argIn.hasKey("atomname"))
    type = ATOMNAME;
  else if (argIn.hasKey("addbond"))
    type = ADDBOND;
  if (type == UNKNOWN) {
    mprinterr("Error: No change type specified.\n");
    return CpptrajState::ERR;
  }
  // Are we doing a COORDS set or a topology?
  Topology* parm = 0;
  std::string crdset = argIn.GetStringKey("crdset");
  if (!crdset.empty()) {
    DataSet_Coords* cset = (DataSet_Coords*)State.DSL().FindSetOfGroup( crdset, DataSet::COORDINATES );
    if (cset == 0) {
      mprinterr("Error: No COORDS set with name '%s' found.\n", crdset.c_str());
      return CpptrajState::ERR;
    }
    parm = cset->TopPtr();
  } else
    parm = State.DSL().GetTopology( argIn );
  if (parm == 0) return CpptrajState::ERR;
  int err = 0;
  switch (type) {
    case RESNAME  : err = ChangeResidueName(*parm, argIn); break;
    case CHAINID  : err = ChangeChainID(*parm, argIn); break;
    case ATOMNAME : err = ChangeAtomName(*parm, argIn); break;
    case ADDBOND  : err = AddBond(*parm, argIn); break;
    case UNKNOWN  : err = 1; // sanity check
  }
  if (err != 0) return CpptrajState::ERR;
  return CpptrajState::OK;
}

// Exec_Change::ChangeResidueName()
int Exec_Change::ChangeResidueName(Topology& topIn, ArgList& argIn)
const
{
  // Name to change to.
  std::string name = argIn.GetStringKey("to");
  if (name.empty()) {
    mprinterr("Error: Specify residue name to change to ('to <name>').\n");
    return 1;
  }
  NameType rname( name );
  // Residues to change
  std::string mexpr = argIn.GetStringKey("from");
  if (mexpr.empty()) {
    mprinterr("Error: Specify residue(s) to change names of ('from <mask>').\n");
    return 1;
  }
  AtomMask mask(mexpr);
  if (topIn.SetupIntegerMask( mask )) return 1;
  if (mask.None()) {
    mprintf("Warning: No atoms selected by mask.\n");
    return 0;
  }
  std::vector<int> resNums = topIn.ResnumsSelectedBy( mask );
  for (std::vector<int>::const_iterator rnum = resNums.begin(); rnum != resNums.end(); ++rnum)
  {
    mprintf("\tChanging residue %s to %s\n", topIn.Res(*rnum).c_str(), *rname);
    topIn.SetRes(*rnum).SetName( rname );
  }
  return 0;
}

// Exec_Change::ChangeChainID()
int Exec_Change::ChangeChainID(Topology& topIn, ArgList& argIn)
const
{
  // ID to change to.
  std::string name = argIn.GetStringKey("to");
  if (name.empty()) {
    mprinterr("Error: Specify chain ID to change to ('to <ID>').\n");
    return 1;
  }
  if (name.size() != 1) {
    mprinterr("Error: Chain ID can only be a single character.\n");
    return 1;
  }
  char cid = name[0];
  // Residues to change
  std::string mexpr = argIn.GetStringKey("of");
  if (mexpr.empty()) {
    mprinterr("Error: Specify residue(s) to change chain IDs of ('of <mask>').\n");
    return 1;
  }
  AtomMask mask(mexpr);
  if (topIn.SetupIntegerMask( mask )) return 1;
  if (mask.None()) {
    mprintf("Warning: No atoms selected by mask.\n");
    return 0;
  }
  std::vector<int> resNums = topIn.ResnumsSelectedBy( mask );
  for (std::vector<int>::const_iterator rnum = resNums.begin(); rnum != resNums.end(); ++rnum)
  {
    mprintf("\tChanging chain ID of residue %s from %c to %c\n", topIn.Res(*rnum).c_str(),
            topIn.Res(*rnum).ChainID(), cid);
    topIn.SetRes(*rnum).SetChainID( cid );
  }
  return 0;
}

// Exec_Change::ChangeAtomName()
int Exec_Change::ChangeAtomName(Topology& topIn, ArgList& argIn)
const
{
  // Name to change to.
  std::string name = argIn.GetStringKey("to");
  if (name.empty()) {
    mprinterr("Error: Specify atom name to change to ('to <name>').\n");
    return 1;
  }
  NameType aname( name );
  // Atoms to change
  std::string mexpr = argIn.GetStringKey("from");
  if (mexpr.empty()) {
    mprinterr("Error: Specify atom(s) to change names of ('from <mask>').\n");
    return 1;
  }
  AtomMask mask(mexpr);
  if (topIn.SetupIntegerMask( mask )) return 1;
  mask.MaskInfo();
  if (mask.None()) {
    mprintf("Warning: No atoms selected by mask.\n");
    return 1;
  }
  for (AtomMask::const_iterator it = mask.begin(); it != mask.end(); ++it)
  {
    mprintf("\tChanging atom %s to %s\n", topIn[*it].c_str(), *aname);
    topIn.SetAtom(*it).SetName( aname );
  }
  return 0;
}

/** Ensure given mask once set up selects exactly one atom. */
int Exec_Change::Setup1atomMask(AtomMask& mask1, Topology const& topIn,
                                std::string const& mask1expr)
{
  if (mask1expr.empty())
    return -1;
  mask1 = AtomMask(mask1expr);
  if (topIn.SetupIntegerMask( mask1 )) return 1;
  if (mask1.Nselected() != 1) {
    mprinterr("Error: Mask must specify only 1 atom, '%s' specifies %i\n",
              mask1.MaskString(), mask1.Nselected());
    return 1;
  }
  return 0;
}

/** \return parameter index of tgtType if found in bonds, -1 otherwise. */
int Exec_Change::FindBondTypeIdx(Topology const& topIn, BondArray const& bonds,
                                 AtomTypeHolder const& tgtType)
{
  int bidx = -1;
  for (BondArray::const_iterator bnd = bonds.begin(); bnd != bonds.end(); ++bnd)
  {
    AtomTypeHolder thisType(2);
    thisType.AddName( topIn[bnd->A1()].Type() );
    thisType.AddName( topIn[bnd->A2()].Type() );
    if (thisType == tgtType) {
      bidx = bnd->Idx();
      break;
    }
  }
    
  return bidx;
}

// Exec_Change::AddBond()
int Exec_Change::AddBond(Topology& topIn, ArgList& argIn) const {
  AtomMask mask1, mask2;
  // Mask1
  int err1 = Setup1atomMask(mask1, topIn, argIn.GetMaskNext());
  int err2 = Setup1atomMask(mask2, topIn, argIn.GetMaskNext());
  if (err1 == -1 || err2 == -1) {
    mprinterr("Error: Must specify 2 masks for 'addbond'\n");
    return 1;
  }
  if (err1 == 1 || err2 == 1) return 1;
  // Get parameter if specified.
  BondParmType bp;
  bool hasBondParm = false;
  double req = argIn.getKeyDouble("req", -1.0);
  if (req > 0.0) {
    hasBondParm = true;
    double rk = argIn.getKeyDouble("rk", 0.0);
    bp = BondParmType(rk, req);
  }
  // Check if bond already exists
  for (Atom::bond_iterator ba = topIn[mask1[0]].bondbegin();
                           ba != topIn[mask1[0]].bondend(); ++ba)
    if (*ba == mask2[0]) {
      mprintf("Warning: Bond already exists between %s and %s\n",
              topIn.TruncResAtomNameNum(mask1[0]).c_str(),
              topIn.TruncResAtomNameNum(mask2[0]).c_str());
      return 0;
  }
  mprintf("\tCreating bond between %s and %s\n",
          topIn.TruncResAtomNameNum(mask1[0]).c_str(),
          topIn.TruncResAtomNameNum(mask2[0]).c_str());
  // See if parameters exist for this bond type
  if (!hasBondParm) {
    int bpidx = -1;
    if (!topIn.BondParm().empty()) {
      AtomTypeHolder tgtType(2);
      tgtType.AddName( topIn[mask1[0]].Type() );
      tgtType.AddName( topIn[mask2[0]].Type() );
      if (topIn[mask1[0]].Element() == Atom::HYDROGEN ||
          topIn[mask2[0]].Element() == Atom::HYDROGEN)
        bpidx = FindBondTypeIdx( topIn, topIn.BondsH(), tgtType );
      else
        bpidx = FindBondTypeIdx( topIn, topIn.Bonds(), tgtType );
    }
    if (bpidx == -1) {
      bp = BondParmType(0.0, Atom::GetBondLength(topIn[mask1[0]].Element(),
                                                 topIn[mask2[0]].Element()));
      mprintf("Warning: No bond current bond parameters exist for these atoms.\n"
              "Warning: Using length %g.\n", bp.Req());
    } else {
      bp = topIn.BondParm()[bpidx];
      mprintf("\tUsing existing bond parameters: length= %g, K= %g\n", bp.Req(), bp.Rk());
    }
  } else {
    mprintf("\tUsing specified bond parameters: length= %g, K= %g\n", bp.Req(), bp.Rk());
  }
  // Try to add the bond
  topIn.AddBond( mask1[0], mask2[0], bp );
  // Regenerate molecule info
  topIn.DetermineMolecules();

  return 0;
}
