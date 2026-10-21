#include "Exec_FixAtomNames.h"
#include "CpptrajStdio.h"

// Exec_FixAtomNames::Help()
void Exec_FixAtomNames::Help() const
{

}

// Exec_FixAtomNames::Execute()
Exec::RetType Exec_FixAtomNames::Execute(CpptrajState& State, ArgList& argIn)
{
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


  Topology const& currentTop = static_cast<Topology const&>(*parm);
  // Select residues to check
  std::string maskStr = argIn.GetStringKey("mask");
  AtomMask mask;
  if (mask.SetMaskString( maskStr )) {
    mprinterr("Error: Invalid mask string '%s'\n", maskStr.c_str());
    return CpptrajState::ERR;
  }
  if (currentTop.SetupIntegerMask( mask )) {
    mprinterr("Error: Could not set up mask '%s'\n", mask.MaskString());
    return CpptrajState::ERR;
  }

  std::vector<int> resnums = currentTop.ResnumsSelectedBy( mask );

  for (std::vector<int>::const_iterator rnum = resnums.begin(); rnum != resnums.end(); ++rnum)
  {
    Residue const& res = currentTop.Res(*rnum);
    std::set<NameType> anames;
    for (int iat = res.FirstAtom(); iat != res.LastAtom(); iat++)
    {
      std::pair<std::set<NameType>::iterator, bool> ret = anames.insert( currentTop[iat].Name() );
      if (!ret.second) {
        mprintf("Warning: '%s': In residue %s duplicate atom name detected: %s\n",
                currentTop.c_str(), currentTop.TruncResNameNum(*rnum).c_str(), *(currentTop[iat].Name()) );
      }
    }
  }

  return CpptrajState::OK;
}
