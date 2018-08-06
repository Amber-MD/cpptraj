#include "Exec_Change.h"
#include "CpptrajStdio.h"

// Exec_Change::Help()
void Exec_Change::Help() const
{
  mprintf("\t[ {%s |\n"
          "\t   crdset <COORDS set>\n"
          "\t{ resname from <mask> to <value> |\n"
          "\t  atomname from <mask> to <value> }\n"
          "  Change specified parts of topology or topology of a COORDS data set.\n",
          DataSetList::TopIdxArgs);
}

// Exec_Change::Execute()
Exec::RetType Exec_Change::Execute(CpptrajState& State, ArgList& argIn)
{
  // Change type
  enum ChangeType { UNKNOWN = 0, RESNAME, ATOMNAME, ADDBOND };
  ChangeType type = UNKNOWN;
  if (argIn.hasKey("resname"))
    type = RESNAME;
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
    DataSet_Coords* cset = (DataSet_Coords*)State.DSL().FindCoordsSet( crdset );
    if (cset == 0) {
      mprinterr("Error: No COORDS set with name '%s' found.\n", crdset.c_str());
      return CpptrajState::ERR;
    }
    parm = cset->TopPtr();
  } else
    parm = State.DSL().GetTopByIndex( argIn );
  if (parm == 0) return CpptrajState::ERR;
  int err = 0;
  switch (type) {
    case RESNAME  : err = ChangeResidueName(*parm, argIn); break;
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
  CharMask mask(mexpr);
  if (topIn.SetupCharMask( mask )) return 1;
  mask.MaskInfo();
  if (mask.None()) {
    mprintf("Warning: No atoms selected by mask.\n");
    return 0;
  }
  for (int res = 0; res != topIn.Nres(); res++)
    if ( mask.AtomsInCharMask( topIn.Res(res).FirstAtom(), topIn.Res(res).LastAtom()-1 ) )
    {
      mprintf("\tChanging residue %s to %s\n", topIn.Res(res).c_str(), *rname);
      topIn.SetRes(res).SetName( rname );
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
    mprinterr("Error: No atoms selected by mask.\n");
    return 1;
  }
  for (AtomMask::const_iterator it = mask.begin(); it != mask.end(); ++it)
  {
    mprintf("\tChanging atom %s to %s\n", topIn[*it].c_str(), *aname);
    topIn.SetAtom(*it).SetName( aname );
  }
  return 0;
}

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
  BondParmType bp;
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
  // Try to add the bond
  topIn.AddBond( mask1[0], mask2[0], bp );
  // Regenerate molecule info
  topIn.DetermineMolecules();

  return 0;
}
