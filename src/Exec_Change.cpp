#include "Exec_Change.h"
#include "CharMask.h"
#include "CpptrajStdio.h"
#include "TypeNameHolder.h"
#include "ParameterTypes.h"
#include "DataSet_1D.h"

// Exec_Change::Help()
void Exec_Change::Help() const
{
  mprintf("\t[ {%s |\n"
          "\t   crdset <COORDS set> ]\n"
          "\t{ resname from <mask> to <value> |\n"
          "\t  chainid of <mask> to <value> |\n"
          "\t  oresnums of <mask> min <range min> max <range max> |\n"
          "\t  icodes of <mask> min <char min> max <char max> resnum <#> |\n"
          "\t  atomname from <mask> to <value> |\n"
          "\t  addbond <mask1> <mask2> [req <length> <rk> <force constant>] |\n"
          "\t  removebonds <mask1> [<mask2>] [out <file>] |\n"
          "\t  bondparm <mask1> [<mask2>] {setrk|scalerk|setreq|scalereq} <value> |\n"
          "\t  {mass|charge} [of <mask>] {to <value> |by <offset> |\n"
          "\t                             byfac <factor> |fromset <data set>} |\n"
          "\t  mergeres firstres <start res#> lastres <stop res#>\n"
          "\t}\n"
          "  Change specified parts of topology or topology of a COORDS data set.\n",
          DataSetList::TopArgs);
}

// Exec_Change::Execute()
Exec::RetType Exec_Change::Execute(CpptrajState& State, ArgList& argIn)
{
  // Change type
  enum ChangeType { UNKNOWN = 0, RESNAME, CHAINID, ORESNUMS, ICODES,
                    ATOMNAME, ADDBOND, REMOVEBONDS, SPLITRES, BONDPARM,
                    MASS, CHARGE, MERGERES };
  ChangeType type = UNKNOWN;
  if (argIn.hasKey("resname"))
    type = RESNAME;
  else if (argIn.hasKey("chainid"))
    type = CHAINID;
  else if (argIn.hasKey("oresnums"))
    type = ORESNUMS;
  else if (argIn.hasKey("icodes"))
    type = ICODES;
  else if (argIn.hasKey("atomname"))
    type = ATOMNAME;
  else if (argIn.hasKey("addbond"))
    type = ADDBOND;
  else if (argIn.hasKey("removebonds"))
    type = REMOVEBONDS;
  else if (argIn.hasKey("bondparm"))
    type = BONDPARM;
  else if (argIn.hasKey("splitres"))
    type = SPLITRES;
  else if (argIn.hasKey("mass"))
    type = MASS;
  else if (argIn.hasKey("charge"))
    type = CHARGE;
  else if (argIn.hasKey("mergeres"))
    type = MERGERES;
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
    mprintf("\tUsing topology from COORDS set '%s'\n", cset->legend());
  } else {
    parm = State.DSL().GetTopology( argIn );
    mprintf("\tUsing topology: %s\n", parm->c_str());
  }
  if (parm == 0) return CpptrajState::ERR;
  int err = 0;
  switch (type) {
    case RESNAME  : err = ChangeResidueName(*parm, argIn); break;
    case CHAINID  : err = ChangeChainID(*parm, argIn); break;
    case ORESNUMS : err = ChangeOresNums(*parm, argIn); break;
    case ICODES   : err = ChangeIcodes(*parm, argIn); break;
    case ATOMNAME : err = ChangeAtomName(*parm, argIn); break;
    case ADDBOND  : err = AddBond(*parm, argIn); break;
    case REMOVEBONDS : err = RemoveBonds(State, *parm, argIn); break;
    case BONDPARM    : err = ChangeBondParameters(*parm, argIn); break;
    case SPLITRES    : err = ChangeSplitRes(*parm, argIn); break;
    case MASS        : err = ChangeMassOrCharge(*parm, argIn, State.DSL(), 0); break;
    case CHARGE      : err = ChangeMassOrCharge(*parm, argIn, State.DSL(), 1); break;
    case MERGERES    : err = ChangeMergeRes(*parm, argIn); break;
    case UNKNOWN  : err = 1; // sanity check
  }
  if (err != 0) return CpptrajState::ERR;
  return CpptrajState::OK;
}

// Exec_Change::ChangeSplitRes()
int Exec_Change::ChangeSplitRes(Topology& topIn, ArgList& argIn)
const
{
  // New residue name
  std::string newname = argIn.GetStringKey("newname");
  if (newname.empty()) {
    mprinterr("Error: splitres: No new name specified.\n");
    return 1;
  }
  // Atoms in residue to split
  std::string maskStr = argIn.GetMaskNext();
  if (maskStr.empty()) {
    mprinterr("Error: splitres: No mask specified.\n");
    return 1;
  }
  AtomMask toSplit;
  if (toSplit.SetMaskString(maskStr)) return 1;
  if (topIn.SetupIntegerMask( toSplit)) return 1;
  if (toSplit.None()) {
    mprinterr("Error: splitres: No atoms selected.\n");
    return 1;
  }
  if (topIn.SplitResidue(toSplit, newname)) {
    mprinterr("Error: splitres failed.\n");
    return 1;
  }
  return 0;
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

// Exec_Change::ChangeOresNums()
int Exec_Change::ChangeOresNums(Topology& topIn, ArgList& argIn)
const
{
  // Residues to change
  std::string mexpr = argIn.GetStringKey("of");
  if (mexpr.empty()) {
    mprinterr("Error: Specify residue(s) to change chain IDs of ('of <mask>').\n");
    return 1;
  }
  AtomMask mask;
  if (mask.SetMaskString(mexpr)) return 1;
  if (topIn.SetupIntegerMask( mask )) return 1;
  std::vector<int> tResIdxs = topIn.ResnumsSelectedBy( mask );
  mprintf("\t%s selects %zu residues.\n", mask.MaskString(), tResIdxs.size());
  if (tResIdxs.empty()) {
    mprinterr("Error: No residues selected by %s\n", mask.MaskString());
    return 1;
  }
  // Number range to change to
  int omin = argIn.getKeyInt("min", 0);
  int omax = argIn.getKeyInt("max", 0);
  unsigned int num_o = (unsigned int)(omax - omin) + 1;
  mprintf("\tOriginal (output) res #s: %i to %i (%u)\n", omin, omax, num_o);
  if (omin > omax) {
    mprinterr("Error: min must be <= max.\n");
    return 1;
  }
/*
  std::string rangearg = argIn.GetStringKey("to");
  if (rangearg.empty()) {
    mprinterr("Error: Specify number range to set for residues.\n");
    return 1;
  }
  Range oResNums;
  if (oResNums.SetRange( rangearg )) {
    mprinterr("Error: Could not set range '%s'\n", rangearg.c_str());
    return 1;
  }
*/
  if (num_o != tResIdxs.size()) {
    mprinterr("Error: # selected residues (%zu) != # provided residue numbers (%u).\n",
              tResIdxs.size(), num_o);
    return 1;
  }
  int currentOnum = omin;
  for (std::vector<int>::const_iterator rnum = tResIdxs.begin();
                                        rnum != tResIdxs.end(); ++rnum, currentOnum++)
  {
    mprintf("\tChanging original res# of residue %s from %i to %i\n",
            topIn.TruncResNameNum(*rnum).c_str(),
            topIn.Res(*rnum).OriginalResNum(), currentOnum);
    topIn.SetRes(*rnum).SetOriginalNum( currentOnum );
  }


  return 0;
}

// Exec_Change::ChangeIcodes()
int Exec_Change::ChangeIcodes(Topology& topIn, ArgList& argIn)
const
{
  // Residues to change
  std::string mexpr = argIn.GetStringKey("of");
  if (mexpr.empty()) {
    mprinterr("Error: Specify residue(s) to change residue insertion codes of ('of <mask>').\n");
    return 1;
  }
  AtomMask mask(mexpr);
  if (topIn.SetupIntegerMask( mask )) return 1;
  if (mask.None()) {
    mprintf("Warning: No atoms selected by mask.\n");
    return 0;
  }
  std::vector<int> resNums = topIn.ResnumsSelectedBy( mask );
  // Original residue number to set
  if (!argIn.Contains("resnum")) {
    mprinterr("Error: Original residue number must be specified: 'resnum <#>'\n");
    return 1;
  }
  int oresnum = argIn.getKeyInt("resnum", 0);
  // Character range to change to
  std::string charstr = argIn.GetStringKey("min");
  if (charstr.empty()) {
    mprinterr("Error: Specify min character to use.\n");
    return 1;
  }
  char cmin = charstr[0];
  charstr = argIn.GetStringKey("max");
  if (charstr.empty()) {
    mprinterr("Error: Specify max character to use.\n");
    return 1;
  }
  char cmax = charstr[0];
  if (cmin == cmax) {
    mprinterr("Error: Min char must be different than max char.\n");
    return 1;
  }
  int dir;
  if (cmin < cmax)
    dir = 1;
  else
    dir = -1;
  int num_c = (int)cmax - (int)cmin + 1;
  if (num_c < 0) num_c = -num_c;
  mprintf("\tOutput residue insertion codes from %c to %c (%u, dir=%i)\n", cmin, cmax, num_c, dir);

  char currentChar = cmin;
  for (std::vector<int>::const_iterator rnum = resNums.begin(); rnum != resNums.end(); ++rnum)
  {
    mprintf("\tChanging insertion code of residue %s from '%i%c' to '%i%c'\n",
            topIn.TruncResNameNum(*rnum).c_str(),
            topIn.Res(*rnum).OriginalResNum(), topIn.Res(*rnum).Icode(), 
            oresnum, currentChar);
    topIn.SetRes(*rnum).SetOriginalNum( oresnum );
    topIn.SetRes(*rnum).SetIcode( currentChar );
    currentChar = (char)((int)currentChar + dir);
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
    mprintf("Warning: Chain IDs with multiple characters are not supported by all output formats.\n");
  }
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
    mprintf("\tChanging chain ID of residue %s from '%s' to '%s'\n",
            topIn.TruncResNameNum(*rnum).c_str(),
            topIn.Res(*rnum).chainID(), name.c_str());
    topIn.SetRes(*rnum).SetChainID( name );
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
                                 TypeNameHolder const& tgtType)
{
  int bidx = -1;
  for (BondArray::const_iterator bnd = bonds.begin(); bnd != bonds.end(); ++bnd)
  {
    TypeNameHolder thisType(2);
    thisType.AddName( topIn[bnd->A1()].Type() );
    thisType.AddName( topIn[bnd->A2()].Type() );
    if (thisType == tgtType) {
      bidx = bnd->Idx();
      break;
    }
  }
    
  return bidx;
}

// Exec_Change::RemoveBonds()
int Exec_Change::RemoveBonds(CpptrajState& State, Topology& topIn, ArgList& argIn) const {
  AtomMask mask1, mask2;
  std::string str1 = argIn.GetMaskNext();
  if (str1.empty()) {
    mprinterr("Error: Must specify at least 1 atom mask.\n");
    return 1;
  }
  if (mask1.SetMaskString( str1 )) return 1;
  if (topIn.SetupIntegerMask( mask1 )) return 1;
  if (mask1.None()) {
    mprinterr("Error: %s selects no atoms.\n", str1.c_str());
    return 1;
  }
  std::string str2 = argIn.GetMaskNext();
  if (!str2.empty()) {
    if (mask2.SetMaskString( str2 )) return 1;
    if (topIn.SetupIntegerMask( mask2 )) return 1;
    if (mask2.None()) {
      mprinterr("Error: %s selects no atoms.\n", str2.c_str());
      return 1;
    }
  }
  CpptrajFile* outfile = State.DFL().AddCpptrajFile(argIn.GetStringKey("out"), "RemovedBonds",
                                                    DataFileList::TEXT, true);
  if (outfile == 0) {
    mprinterr("Internal Error: RemoveBonds could not get an output file.\n");
    return 1;
  }
  const char* prefix = "";
  if (outfile->IsStream())
    prefix = "\t\t";

  if (str2.empty()) {
    mprintf("\tRemoving bonds to atoms selected by %s (%i atoms).\n", 
            str1.c_str(), mask1.Nselected());
    for (AtomMask::const_iterator atm = mask1.begin(); atm != mask1.end(); ++atm) {
      std::string atmStr = topIn.ResNameNumAtomNameNum(*atm);
      // Make a copy of the atoms bonds array because it will be modified.
      std::vector<int> atoms = topIn[*atm].BondIdxArray();
      for (std::vector<int>::const_iterator bnd = atoms.begin(); bnd != atoms.end(); ++bnd)
      {
        int ret = topIn.RemoveBond(*atm, *bnd);
        if (ret == 0)
          outfile->Printf("%s%s to %s\n", prefix, atmStr.c_str(),
                          topIn.ResNameNumAtomNameNum(*bnd).c_str());
      }
    }
  } else {
    mprintf("\tRemoving any bonds between atoms selected by %s (%i atoms)\n"
            "\tand %s (%i atoms).\n", str1.c_str(), mask1.Nselected(),
            str2.c_str(), mask2.Nselected());
    for (AtomMask::const_iterator atm1 = mask1.begin(); atm1 != mask1.end(); ++atm1) {
      std::string atmStr = topIn.ResNameNumAtomNameNum(*atm1);
      for (AtomMask::const_iterator atm2 = mask2.begin(); atm2 != mask2.end(); ++atm2) {
        int ret = topIn.RemoveBond(*atm1, *atm2);
        if (ret == 0)
          outfile->Printf("%s%s to %s\n", prefix,atmStr.c_str(),
                          topIn.ResNameNumAtomNameNum(*atm2).c_str());
      }
    }
  }
  // Since molecule info has likely changed, re-determine
  topIn.DetermineMolecules();

  return 0;
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
      TypeNameHolder tgtType(2);
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

/** Change/scale bond force constants. */
int Exec_Change::ChangeBondParameters(Topology& topIn, ArgList& argIn) const {
  CharMask mask1, mask2;
  std::string str1 = argIn.GetMaskNext();
  if (str1.empty()) {
    mprinterr("Error: Must specify at least 1 atom mask.\n");
    return 1;
  }
  if (mask1.SetMaskString( str1 )) return 1;
  if (topIn.SetupCharMask( mask1 )) return 1;
  if (mask1.None()) {
    mprinterr("Error: %s selects no atoms.\n", str1.c_str());
    return 1;
  }
  std::string str2 = argIn.GetMaskNext();
  if (!str2.empty()) {
    if (mask2.SetMaskString( str2 )) return 1;
    if (topIn.SetupCharMask( mask2 )) return 1;
    if (mask2.None()) {
      mprinterr("Error: %s selects no atoms.\n", str2.c_str());
      return 1;
    }
  }
  // Determine if we are setting or scaling
  enum ModeType { SET_RK = 0, SET_REQ, SCALE_RK, SCALE_REQ };
  ModeType mode;
  double dval = 0;
  if (argIn.Contains("setrk")) {
    mode = SET_RK;
    dval = argIn.getKeyDouble("setrk", 0);
  } else if (argIn.Contains("scalerk")) {
    mode = SCALE_RK;
    dval = argIn.getKeyDouble("scalerk", 0);
  } else if (argIn.Contains("setreq")) {
    mode = SET_REQ;
    dval = argIn.getKeyDouble("setreq", 0);
  } else if (argIn.Contains("scalereq")) {
    mode = SCALE_REQ;
    dval = argIn.getKeyDouble("scalereq", 0);
  } else {
    mprinterr("Error: Must specify one of 'setrk', 'scalerk', 'setreq', 'scalereq'.\n");
    return 1;
  }

  static const char* ModeTypeStr[] = {"set RK to", "set REQ to", "scale RK by", "scale REQ by"};
  mprintf("\tWill %s %f\n", ModeTypeStr[mode], dval);

  // Get the list of bonds to modify
  BondArray bondsToModify;
  if (mask2.MaskStringSet()) {
    // Atom 1 must be in mask1, atom 2 in mask2
    for (BondArray::const_iterator it = topIn.Bonds().begin(); it != topIn.Bonds().end(); ++it)
      if (mask1.AtomInCharMask(it->A1()) && mask2.AtomInCharMask(it->A2()))
        bondsToModify.push_back( *it );
     for (BondArray::const_iterator it = topIn.BondsH().begin(); it != topIn.BondsH().end(); ++it)
      if (mask1.AtomInCharMask(it->A1()) && mask2.AtomInCharMask(it->A2()))
        bondsToModify.push_back( *it );
  } else {
    for (BondArray::const_iterator it = topIn.Bonds().begin(); it != topIn.Bonds().end(); ++it)
      if (mask1.AtomInCharMask(it->A1()) && mask1.AtomInCharMask(it->A2()))
        bondsToModify.push_back( *it );
     for (BondArray::const_iterator it = topIn.BondsH().begin(); it != topIn.BondsH().end(); ++it)
      if (mask1.AtomInCharMask(it->A1()) && mask1.AtomInCharMask(it->A2()))
        bondsToModify.push_back( *it );
  }
  mprintf("\tBonds to modify (%zu):\n", bondsToModify.size());
  for (BondArray::const_iterator it = bondsToModify.begin(); it != bondsToModify.end(); ++it)
    mprintf("\t\t%s to %s\n",
            topIn.TruncResAtomName( it->A1() ).c_str(),
            topIn.TruncResAtomName( it->A2() ).c_str());

  // Remove, modify, add back.
  for (BondArray::const_iterator it = bondsToModify.begin(); it != bondsToModify.end(); ++it)
  {
    BondParmType bp = topIn.BondParm()[it->Idx()];
    topIn.RemoveBond(it->A1(), it->A2());
    switch (mode) {
      case SET_RK    : bp.SetRk( dval ); break;
      case SCALE_RK  : bp.SetRk( bp.Rk() * dval ); break;
      case SET_REQ   : bp.SetReq( dval ); break;
      case SCALE_REQ : bp.SetReq( bp.Req() * dval ); break;
    }
    topIn.AddBond(it->A1(), it->A2(), bp);
  }

  return 0;
}

/** Function to change specific value in topology. */
void Exec_Change::changeTopVal(Topology& topIn, int atnum, ChangeType typeIn, double newVal)
{
  Atom& currentAtom = topIn.SetAtom(atnum);
  double oldVal;
  switch (typeIn) {
    case MASS_TO :
      oldVal = currentAtom.Mass();
      mprintf("\tChanging mass of atom '%s' from %g to %g\n", topIn.AtomMaskName(atnum).c_str(), oldVal, newVal);
      currentAtom.SetMass( newVal );
      break;
    case CHARGE_TO :
      oldVal = currentAtom.Charge();
      mprintf("\tChanging charge of atom '%s' from %g to %g\n", topIn.AtomMaskName(atnum).c_str(), oldVal, newVal);
      currentAtom.SetCharge( newVal );
      break;
    case MASS_BY :
      oldVal = currentAtom.Mass();
      mprintf("\tChanging mass of atom '%s' by %g to %g\n", topIn.AtomMaskName(atnum).c_str(), newVal, oldVal+newVal);
      currentAtom.SetMass( oldVal + newVal );
      break;
    case CHARGE_BY :
      oldVal = currentAtom.Charge();
      mprintf("\tChanging charge of atom '%s' by %g to %g\n", topIn.AtomMaskName(atnum).c_str(), newVal, oldVal+newVal);
      currentAtom.SetCharge( oldVal + newVal );
      break;
    case MASS_BYFAC :
      oldVal = currentAtom.Mass();
      mprintf("\tChanging mass of atom '%s' by factor of %g to %g\n", topIn.AtomMaskName(atnum).c_str(), newVal, oldVal*newVal);
      currentAtom.SetMass( oldVal * newVal );
      break;
    case CHARGE_BYFAC :
      oldVal = currentAtom.Charge();
      mprintf("\tChanging charge of atom '%s' by factor of %g to %g\n", topIn.AtomMaskName(atnum).c_str(), newVal, oldVal*newVal);
      currentAtom.SetCharge( oldVal * newVal );
      break;
  }
}

/** Change mass/charge in topology.
  * \param typeIn: 0=mass, 1=charge
  */
int Exec_Change::ChangeMassOrCharge(Topology& topIn, ArgList& argIn,
                                    DataSetList const& DSL, int typeIn)
const
{
  // sanity check
  if (typeIn > 1 || typeIn < 0) {
    mprinterr("Internal Error: typeIn is not 0 or 1.\n");
    return 1;
  }
  static const char* desc[] = { "mass", "charge" };

  std::string maskExpression = argIn.GetStringKey("of");
  AtomMask atomsToChange;
  if (atomsToChange.SetMaskString( maskExpression )) {
    mprinterr("Error: Could not set mask expression.\n");
    return 1;
  }
  if (topIn.SetupIntegerMask( atomsToChange )) {
    mprinterr("Error: Could not set up mask.\n");
    return 1;
  }
  if (atomsToChange.None()) {
    mprintf("Warning: Mask '%s' selects no atoms.\n", atomsToChange.MaskString());
    return 0;
  }
  atomsToChange.MaskInfo();

  std::string fromSet = argIn.GetStringKey("fromset");
  if (!fromSet.empty()) {
    // Get charges from a data set
    DataSet* ds = DSL.GetDataSet( fromSet );
    if (ds == 0) {
      mprinterr("Error: No set selected by '%s'\n", fromSet.c_str());
      return 1;
    }
    if (ds->Group() != DataSet::SCALAR_1D) {
      mprinterr("Error: Data set '%s' is not scalar 1D.\n", ds->legend());
      return 1;
    }
    mprintf("\tUsing data from '%s' for %s.\n", ds->legend(), desc[typeIn]);
    if (ds->Size() != (unsigned int)atomsToChange.Nselected()) {
      mprinterr("Error: %i atoms to change %s of, but set '%s' has %zu elements.\n",
                atomsToChange.Nselected(), desc[typeIn], ds->legend(), ds->Size());
      return 1;
    }
    ChangeType ctype;
    if (typeIn == 0)
      ctype = MASS_TO;
    else
      ctype = CHARGE_TO;
    DataSet_1D const& dset = static_cast<DataSet_1D const&>( *ds );
    for (int idx = 0; idx != atomsToChange.Nselected(); idx++) {
      changeTopVal(topIn, atomsToChange[idx], ctype, dset.Dval(idx));
    }
  } else {
    // Get charge/offset from command line
    ChangeType ctype;
    bool change_to = argIn.Contains("to");
    bool change_by = argIn.Contains("by");
    bool change_byfac = argIn.Contains("byfac");
    int nchange = 0;
    if (change_to) nchange++;
    if (change_by) nchange++;
    if (change_byfac) nchange++;
    if (nchange < 1) {
      mprinterr("Error: Expected either 'fromset', 'to', 'by', or 'byfac' for 'change %s'.\n", desc[typeIn]);
      return 1;
    } else if (nchange > 1) {
      mprinterr("Error: Specify either 'to' or 'by', not both.\n");
      return 1;
    }
    if (change_to) {
      if (typeIn == 0)
        ctype = MASS_TO;
      else
        ctype = CHARGE_TO;
      double newVal = argIn.getKeyDouble("to", 0.0);
      for (AtomMask::const_iterator at = atomsToChange.begin(); at != atomsToChange.end(); ++at) {
        changeTopVal(topIn, *at, ctype, newVal);
      }
    } else if (change_by) {
      if (typeIn == 0)
        ctype = MASS_BY;
      else
        ctype = CHARGE_BY;
      double offset = argIn.getKeyDouble("by", 0.0);
      for (AtomMask::const_iterator at = atomsToChange.begin(); at != atomsToChange.end(); ++at) {
        changeTopVal(topIn, *at, ctype, offset);
      }
    } else if (change_byfac) {
      if (typeIn == 0)
        ctype = MASS_BYFAC;
      else
        ctype = CHARGE_BYFAC;
      double fac = argIn.getKeyDouble("byfac", 0.0);
      for (AtomMask::const_iterator at = atomsToChange.begin(); at != atomsToChange.end(); ++at) {
        changeTopVal(topIn, *at, ctype, fac);
      }
    }
  }

  return 0;
}

/** Merge residues in a range. */
int Exec_Change::ChangeMergeRes(Topology& topIn, ArgList& argIn) const {
  int firstres = argIn.getKeyInt("firstres", 0) - 1;
  if (firstres < 0) {
    mprinterr("Error: 'firstres' not specified or < 1.\n");
    return 1;
  }
  int lastres = argIn.getKeyInt("lastres", 0) - 1;
  if (lastres < 0) {
    mprinterr("Error: 'lastres' not specified or < 1.\n");
    return 1;
  }
  return topIn.MergeResidues(firstres, lastres);
}
