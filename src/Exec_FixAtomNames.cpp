#include "Exec_FixAtomNames.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // convertToInteger
#include <cctype> // isalpha, isdigit

// -----------------------------------------------------------------------------
void Exec_FixAtomNames::AtomName::PrintName() const {
  mprintf("\t");
  for (NameTokenArray::const_iterator it = tokens_.begin(); it != tokens_.end(); ++it) {
    if (it->TokenType() == NameToken::INT)
      mprintf("%i", it->Ival());
    else
      mprintf("%s", it->Sval().c_str());
  }
  mprintf("\n");
    
  for (NameTokenArray::const_iterator it = tokens_.begin(); it != tokens_.end(); ++it) {
    switch (it->TokenType()) {
      case NameToken::INT    : mprintf("\t  INTEGER : %i\n", it->Ival()); break;
      case NameToken::STRING : mprintf("\t  STRING  : %s\n", it->Sval().c_str()); break;
      case NameToken::OTHER  : mprintf("\t  OTHER   : %s\n", it->Sval().c_str()); break;
      default: break;
    }
  }
}
// -----------------------------------------------------------------------------

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
        fixDuplicate(*parm, iat, *rnum);
      }
    }
  }

  return CpptrajState::OK;
}

Exec_FixAtomNames::NameToken::Type Exec_FixAtomNames::get_type(char cIn)
{
  if (isalpha(cIn))      return NameToken::STRING;
  else if (isdigit(cIn)) return NameToken::INT;
  else                   return NameToken::OTHER;
}

/** Dissect name into an array of tokens (AtomName) */
Exec_FixAtomNames::AtomName Exec_FixAtomNames::dissectName(NameType const& anameIn)
{
  std::string aname = anameIn.Truncated();
  std::string currentToken;
  AtomName atomName;
  // Sanity check
  if (aname.empty()) {
    mprinterr("Internal Error: Exec_FixAtomNames::dissectName(): Atom name is empty.\n");
    return atomName;
  }
  // Figure out the type of the first token
  std::string::const_iterator it = aname.begin();
  NameToken::Type currentType = get_type( *it );
  for (std::string::const_iterator it = aname.begin(); it != aname.end(); ++it)
  {
    currentToken.push_back( *it );
    // If the next char is a different token type or the end, add the token
    NameToken::Type nextType = NameToken::UNKNOWN;
    bool push_token = false;
    if ( (it+1) == aname.end())
      push_token = true;
    else {
      nextType = get_type( *(it+1) );
      if (nextType != currentType)
        push_token = true;
    }
    if (push_token) {
      if (currentType == NameToken::STRING)
        atomName.AddToken( NameToken( currentToken ) );
      else if (currentType == NameToken::INT)
        atomName.AddToken( NameToken( convertToInteger(currentToken) ) );
      else
        atomName.AddToken( NameToken( NameToken::OTHER, currentToken ) );
      currentToken.clear();
    }
    // Update the curren type with next type
    currentType = nextType;
  }
  return atomName;
}

/** Given the index of an atom that has a duplicated name in a reside, change
  * it to a unique name in that residue.
  */
int Exec_FixAtomNames::fixDuplicate(Topology& topIn, int duplicateIdx, int rnum)
const
{
  Residue const& res = topIn.Res(rnum);
  // Dissect this name.
  AtomName duplicateName = dissectName( topIn[duplicateIdx].Name() );
  // DEBUG
  duplicateName.PrintName();
  // Create a set of current atom names in the residue. They will all be unique.
  std::set<NameType> anames;
  for (int iat = res.FirstAtom(); iat != res.LastAtom(); iat++)
    anames.insert( topIn[iat].Name() );
  // DEBUG TEST dissect all names
  for (std::set<NameType>::const_iterator it = anames.begin(); it != anames.end(); ++it) {
    AtomName atomName = dissectName( *it );
    atomName.PrintName();
  }

  return 0;
}
