#include "Exec_FixAtomNames.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // convertToInteger
#include <cctype> // isalpha, isdigit

// -----------------------------------------------------------------------------
/** \return name width in characters */
int Exec_FixAtomNames::AtomName::NameWidth() const {
  int width = 0;
  for (NameTokenArray::const_iterator it = tokens_.begin(); it != tokens_.end(); ++it)
  {
    if (it->TokenType() == NameToken::INT)
      width += DigitWidth( it->Ival() );
    else
      width += it->Sval().size();
  }
  return width;
}

/** Print name to stdout */
void Exec_FixAtomNames::AtomName::PrintName() const {
  mprintf("\t Current Name: ");
  for (NameTokenArray::const_iterator it = tokens_.begin(); it != tokens_.end(); ++it) {
    if (it->TokenType() == NameToken::INT)
      mprintf("%i", it->Ival());
    else
      mprintf("%s", it->Sval().c_str());
  }
  mprintf(" (%i)\n", NameWidth());
    
  for (NameTokenArray::const_iterator it = tokens_.begin(); it != tokens_.end(); ++it) {
    switch (it->TokenType()) {
      case NameToken::INT    : mprintf("\t  INTEGER : %i\n", it->Ival()); break;
      case NameToken::STRING : mprintf("\t  STRING  : %s\n", it->Sval().c_str()); break;
      case NameToken::OTHER  : mprintf("\t  OTHER   : %s\n", it->Sval().c_str()); break;
      default: break;
    }
  }
}

/** Convert to a string */
std::string Exec_FixAtomNames::AtomName::NameString() const {
  std::string aname;
  for (NameTokenArray::const_iterator it = tokens_.begin(); it != tokens_.end(); ++it) {
    if (it->TokenType() == NameToken::INT)
      aname.append( integerToString(it->Ival()) );
    else
      aname.append( it->Sval() );
  }
  return aname;
}

/** \return Token type from character */
Exec_FixAtomNames::NameToken::Type Exec_FixAtomNames::get_type(char cIn)
{
  if (isalpha(cIn))      return NameToken::STRING;
  else if (isdigit(cIn)) return NameToken::INT;
  else                   return NameToken::OTHER;
}

/** Dissect name into an array of tokens (AtomName) */
int Exec_FixAtomNames::AtomName::DissectName(NameType const& anameIn)
{
  std::string aname = anameIn.Truncated();
  std::string currentToken;
  tokens_.clear(); 
  // Sanity check
  if (aname.empty()) {
    mprinterr("Internal Error: Exec_FixAtomNames::AtomName::DissectName(): Atom name is empty.\n");
    return 1;
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
        tokens_.push_back( NameToken( currentToken ) );
      else if (currentType == NameToken::INT)
        tokens_.push_back( NameToken( convertToInteger(currentToken) ) );
      else
        tokens_.push_back( NameToken( NameToken::OTHER, currentToken ) );
      currentToken.clear();
    }
    // Update the curren type with next type
    currentType = nextType;
  }
  return 0;
}

/** Naively increment the atom name. */
int Exec_FixAtomNames::AtomName::IncrementName(int maxWidth) {
  int currentWidth = NameWidth();
  int errval = 1;
  // If we have no integer token, add one before any OTHER tokens.
  bool has_int = false;
  for (NameTokenArray::const_iterator it = tokens_.begin(); it != tokens_.end(); ++it) {
    if (it->TokenType() == NameToken::INT) {
      has_int = true;
      break;
    }
  }
  if (!has_int) {
    //FIXME shouldnt put this after a prime
    if (currentWidth == maxWidth) {
      mprinterr("Error: Name %s is already at max width (%i), cannot add digit.\n", NameString().c_str(), maxWidth);
      return 1;
    }
    tokens_.push_back( NameToken(1) );
    return 0;
  }

  // Start at the last token and work up
  for (NameTokenArray::reverse_iterator rit = tokens_.rbegin(); rit != tokens_.rend(); ++rit)
  {
    if ( rit->TokenType() == NameToken::INT) {
      // What is the max digit width we can accomodate?
      int digitWidth = DigitWidth( rit->Ival() );
      int otherWidth = currentWidth - digitWidth;
      int availableWidth = maxWidth - otherWidth;
      int maxDigit = pow(10, availableWidth) - 1;
      //mprintf("DEBUG: digit=%i MaxDigit = %i\n", rit->Ival(), maxDigit);
      int newDigit = rit->Ival() + 1;
      if (newDigit <= maxDigit) {
        // Increment should be ok. Can bail out now.
        rit->SetIval( newDigit );
        errval = 0;
        break;
      } else {
        // Reset to 1. Have to increment something else.
        rit->SetIval( 1 );
      }
    }
  }

  return errval;
}

// -----------------------------------------------------------------------------

static inline void WARNING() {
  mprintf("WARNING: THIS COMMAND IS STILL EXPERIMENTAL. USE AT YOUR OWN RISK.\n");
}

// Exec_FixAtomNames::Help()
void Exec_FixAtomNames::Help() const
{
  mprintf("\t{ %s | crdset <COORDS SET>} [<mask>]\n"
          "  Attempt to fix duplicated atom names in residues of the specified topology.\n",
          DataSetList::TopArgs);
  WARNING();
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
  WARNING();

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
  mask.MaskInfo();

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
        if (fixDuplicate(*parm, iat, *rnum)) {
          mprinterr("Error: Could not fix duplicated name %s\n", *(*(ret.first)));
          return CpptrajState::ERR;
        }
        // Add the fixed (not duplicated) name
        anames.insert( currentTop[iat].Name() );
      }
    }
    // DEBUG TEST dissect all names
//    for (std::set<NameType>::const_iterator it = anames.begin(); it != anames.end(); ++it) {
//      AtomName atomName;
//      if (atomName.DissectName( *it )) return CpptrajState::ERR;
//      atomName.PrintName();
//    }
  }

  return CpptrajState::OK;
}

/** Given the index of an atom that has a duplicated name in a reside, change
  * it to a unique name in that residue.
  */
int Exec_FixAtomNames::fixDuplicate(Topology& topIn, int duplicateIdx, int rnum)
const
{
  Residue const& res = topIn.Res(rnum);
  // Dissect this name.
  AtomName duplicateName;
  if (duplicateName.DissectName( topIn[duplicateIdx].Name() )) return 1;
  // DEBUG
  duplicateName.PrintName();
  // Create a set of current atom names in the residue. They will all be unique.
  std::set<NameType> anames;
  for (int iat = res.FirstAtom(); iat != res.LastAtom(); iat++)
    anames.insert( topIn[iat].Name() );
  // Try to increment the name until we do not get a clash
  bool tryNewNames = true;
  while (tryNewNames) {
    if (duplicateName.IncrementName(4) != 0) {
      mprinterr("Error: Ran out of names to increment for %s.\n", topIn.AtomMaskName(duplicateIdx).c_str());
      return 1;
    }
    NameType newName( duplicateName.NameString() );
    std::set<NameType>::iterator it = anames.find( newName );
    if (it == anames.end()) {
      // Name is no longer duplicated!
      mprintf("\tNew name: %s\n", *newName);
      topIn.SetAtom(duplicateIdx).SetName( newName );
      tryNewNames = false;
    } else
      mprintf("\tName %s already exists.\n", *newName);
  }

  return 0;
}
