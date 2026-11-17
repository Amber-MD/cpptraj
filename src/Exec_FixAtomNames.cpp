#include "Exec_FixAtomNames.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // convertToInteger
#include <cctype> // isalpha, isdigit
#include <algorithm> // sort

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
void Exec_FixAtomNames::AtomName::PrintName(const char* desc) const {
  mprintf("\t %s: ", desc);
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

void Exec_FixAtomNames::AtomName::PrintName() const {
  PrintName("Current Name");
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

/** Remove any leading digits */
int Exec_FixAtomNames::AtomName::RemoveLeadingDigit() {
  std::vector<bool> toRemove(tokens_.size(), false);
  unsigned int nToRemove = 0;
  unsigned int idx = 0;
  for (NameTokenArray::const_iterator it = tokens_.begin(); it != tokens_.end(); ++it, idx++)
  {
    if (it->TokenType() != NameToken::INT) break;
    toRemove[idx] = true;
    nToRemove++;
  }
  if (nToRemove == tokens_.size()) {
    mprintf("Warning: Would remove all tokens. Not removing leading digit.\n");
    return 0;
  }
  NameTokenArray newTokens;
  for (idx = 0; idx < toRemove.size(); idx++)
  {
    if (!toRemove[idx])
      newTokens.push_back( tokens_[idx] );
  }
  tokens_ = newTokens;
  return 0;
}

/** Set or add terminal number */
int Exec_FixAtomNames::AtomName::SetTerminalNumber(int newDigit, int maxWidth)
{
  int currentWidth = NameWidth();
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
      if (newDigit <= maxDigit) {
        // Digit should be ok. Can bail out now.
        rit->SetIval( newDigit );
        //errval = 0;
        break;
      } else {
        mprintf("Warning: Terminal digit %i too large.\n", newDigit);
        return 1;
      }
    }
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

const int Exec_FixAtomNames::maxNameWidth_ = 4;

static inline void WARNING() {
  mprintf("WARNING: THIS COMMAND IS STILL EXPERIMENTAL. USE AT YOUR OWN RISK.\n");
}

// Exec_FixAtomNames::Help()
void Exec_FixAtomNames::Help() const
{
  mprintf("\t{ %s | crdset <COORDS SET>} [<mask>] [renameh]\n"
          "  Attempt to fix duplicated atom names in residues of the specified topology.\n",
          DataSetList::TopArgs);
  mprintf("  If a heavy atom is renamed, any bonded hydrogens will be renamed\n"
          "  according to the PDB V3 standard.\n"
          "  If 'renameh' is specified, attempt to rename all hydrogens according\n"
          "  to the PDB V3 standard whether or not the heavy atom is renamed.\n");
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

  // Do we want to rename hydrogens
  bool rename_hydrogens = argIn.hasKey("renameh");
  if (rename_hydrogens)
    mprintf("\tWill rename hydrogens after duplicates are resolved.\n");

  // Determine which residues should be checked
  Iarray resnums = currentTop.ResnumsSelectedBy( mask );
  // Loop over all residues to be checked
  for (Iarray::const_iterator rnum = resnums.begin(); rnum != resnums.end(); ++rnum)
  {
    // Hold atom name and list of duplicated indices
    if (findDuplicatesInResidue(currentTop, *rnum)) {
      // Sort by priority
      std::vector<DupName> duplicates;
      for (DuplicateMap::const_iterator it = duplicateMap_.begin(); it != duplicateMap_.end(); ++it)
      {
        if (it->second.size() > 1) {
          int wt;
          if (currentTop[it->second.front()].Element() == Atom::CARBON)
            wt = 9999;
          else
            wt = currentTop[it->second.front()].AtomicNumber();
          duplicates.push_back( DupName(it->first, it->second, wt) );
        }
      }
      std::sort( duplicates.begin(), duplicates.end() );
      // print all duplicates
      for (std::vector<DupName>::const_iterator it = duplicates.begin(); it != duplicates.end(); ++it)
      {
        mprintf("\tDuplicated name: %s @%s, indices", currentTop.TruncResNameNum(*rnum).c_str(), *(it->aname_));
        for (Iarray::const_iterator jt = it->idxs_.begin(); jt != it->idxs_.end(); ++jt)
          mprintf(" %i", *jt + 1);
          //mprintf(" %s", currentTop.AtomMaskName(*jt).c_str());
        mprintf("\n");
        if (fixDuplicateSet(*parm, it->idxs_, *rnum)) return CpptrajState::ERR;
      }
      // Sanity check - should be no more duplicates
      if (findDuplicatesInResidue(currentTop, *rnum)) {
        for (DuplicateMap::const_iterator it = duplicateMap_.begin(); it != duplicateMap_.end(); ++it)
        {
          if (it->second.size() > 1) {
            mprinterr("Error: %s is still duplicated (%zu times)\n", *(it->first), it->second.size());
          }
        }
        return CpptrajState::ERR;
      }
    } // END residue has duplicates
    // Rename hydrogens if requested
    if (rename_hydrogens) {
      if (RenameAllHydrogens(*parm, *rnum)) {
        mprinterr("Error: Could not rename hydrogens in residue %s\n", currentTop.TruncResNameNum(*rnum).c_str());
        return CpptrajState::ERR;
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

// ---------------------------------------------------------
/** Create a count of atom names in a residue, look for duplicated names.
  * \return true if any atom names are duplicated
  */
int Exec_FixAtomNames::findDuplicatesInResidue(Topology const& currentTop, int rnum)
{
  bool has_duplicates = false;
  // Hold atom name and list of duplicated indices
  duplicateMap_.clear();
  Residue const& res = currentTop.Res(rnum);
  for (int iat = res.FirstAtom(); iat != res.LastAtom(); iat++)
  {
    DuplicateMap::iterator it = duplicateMap_.find( currentTop[iat].Name() );
    if (it == duplicateMap_.end()) {
      // New name
      duplicateMap_.insert( DuplicatePair(currentTop[iat].Name(), Iarray(1, iat)) );
    } else {
      // Duplicated name
      mprintf("Warning: '%s': In residue %s duplicate atom name detected: %s\n",
              currentTop.c_str(), currentTop.TruncResNameNum(rnum).c_str(), *(currentTop[iat].Name()) );
      it->second.push_back( iat );
      has_duplicates = true;
    }
  }
  if (has_duplicates) return 1;
  return 0;
}

/** Try to rename all hydrogens in the specifed resdiue to the PDBV3 standard. */
int Exec_FixAtomNames::RenameAllHydrogens(Topology& topIn, int rnum) {
  // Create a set of current atom names in the residue. They will all be unique.
  Residue const& res = topIn.Res(rnum);
  //std::set<NameType> anames;
  //for (int iat = res.FirstAtom(); iat != res.LastAtom(); iat++)
  //  anames.insert( topIn[iat].Name() );

  for (int iat = res.FirstAtom(); iat != res.LastAtom(); iat++)
  {
    Atom const& currentAtom = topIn[iat];
    if (currentAtom.Element() != Atom::HYDROGEN) {
      // Is this atom bonded to hydrogens?
      bool has_h = false;
      for (Atom::bond_iterator bat = currentAtom.bondbegin(); bat != currentAtom.bondend(); ++bat)
      {
        if (topIn[*bat].Element() == Atom::HYDROGEN) {
          has_h = true;
          break;
        }
      }
      if (has_h) {
        mprintf("\tRenaming hydrogens bonded to %s\n", topIn.AtomMaskName(iat).c_str());
        RenameBondedHydrogens( topIn, iat );
        //if (RenameBondedHydrogens( topIn, iat )) {
          // We are not just adding hydrogen names, we may have taken some away.
          // Just regenerate the atom name list.
          //anames.clear();
          //for (int iat = res.FirstAtom(); iat != res.LastAtom(); iat++)
          //  anames.insert( topIn[iat].Name() );
        //}
      }
    } // END current element is not hydrogen
  } // END loop over residue atoms
  // Sanity check
  if (findDuplicatesInResidue(topIn, rnum)) {
    mprinterr("Error: Hydrogen renaming resulted in duplicate atoms.\n");
    return 1;
  }
  return 0;
}

/** Change atom name to another atom name.
  * If modxna info is present ensure the modxna header is changed as well.
  */
int Exec_FixAtomNames::change_aname(Topology& topIn, int atomIdx, NameType const& newName)
{
  // Sanity check. Proposed name should not already exist.
  DuplicateMap::iterator existing = duplicateMap_.find( newName );
  if (existing != duplicateMap_.end()) {
    mprinterr("Internal Error: change_aname() called with existing atom name %s\n", *newName);
    return 1;
  }
  // Get the current (existing) name
  unsigned int oldNameCount = 0;
  existing = duplicateMap_.find( topIn[atomIdx].Name() );
  if ( existing == duplicateMap_.end()) {
    mprinterr("Internal Error: Existing name %s for atom %i not found (new name %s)\n", *(topIn[atomIdx].Name()), atomIdx+1, *newName);
  } else {
    // Remove the atom index from the index array for this name
    Iarray newArray;
    for (Iarray::const_iterator it = existing->second.begin(); it != existing->second.end(); ++it)
      if (*it != atomIdx)
        newArray.push_back( *it );
    oldNameCount = newArray.size();
    // Are all the names gone?
    if (newArray.empty()) {
      mprintf("\tNo more atoms named %s\n", *(existing->first));
      duplicateMap_.erase( existing );
    } else
      existing->second = newArray;
  }
  // Check if the current name corresponds to a ModXNA mask.
  if (topIn.Modxna().HasModxna()) {
    int mtype = topIn.Modxna().FindMaskWithAtom( topIn[atomIdx].Name().Truncated() );
    if (mtype > -1) {
      if (oldNameCount < 1) {
        // Atom that was renamed was unique. Safe to rename the header.
        mprintf("\tOld atom name '%s' in ModXNA mask %s %s will be renamed to %s\n",
                *(topIn[atomIdx].Name()),
                topIn.Modxna().maskTypeStr( (ModXNA_Info::MaskType)mtype ),
                topIn.Modxna().AtomMaskString( (ModXNA_Info::MaskType)mtype ).c_str(),
                *newName);
        topIn.ChangeModxna().RenameAtomInMasks( topIn[atomIdx].Name().Truncated(), newName.Truncated() );
      } else
        // Atom that was renamed had at least 1 other duplicate name. Cannot be sure which
        // is the "right" one, so just warn.
        mprintf("Warning: Old duplicated atom name '%s' is in ModXNA mask %s %s ; may need to update ModXNA header.\n",
                *(topIn[atomIdx].Name()),
                topIn.Modxna().maskTypeStr( (ModXNA_Info::MaskType)mtype ),
                topIn.Modxna().AtomMaskString( (ModXNA_Info::MaskType)mtype ).c_str());
    }
  }
  // Set the atom with the new name
  topIn.SetAtom(atomIdx).SetName( newName );
  // Add the new name to the map
  duplicateMap_.insert( DuplicatePair(newName, Iarray(1, atomIdx)) );
  return 0;
}

/** Given the index of an atom with a duplicated name and the index of
  * an atom it is bonded to, see if we can resolve the clash by
  * making the terminal number match. This will only work if the atom
  * is the only one of its kind bonded to the bond atom.
  */
int Exec_FixAtomNames::MakeTerminalNumberMatch(Topology& topIn, int atomIdx, int bondedIdx, 
                                               AtomName& currentName, AtomName const& bondedName)
{
  Atom const& currentAtom = topIn[atomIdx];
  Atom const& bondedAtom = topIn[bondedIdx];

  mprintf("\tDuplicated Atom %s, bonded atom %s\n", topIn.AtomMaskName(atomIdx).c_str(), topIn.AtomMaskName(bondedIdx).c_str());
  int eltCount = 0;
  for (Atom::bond_iterator bat = bondedAtom.bondbegin(); bat != bondedAtom.bondend(); ++bat)
    if (topIn[*bat].Element() == currentAtom.Element())
      eltCount++;
  mprintf("\tElt count: %i\n", eltCount);
  if (eltCount > 1) return 0;

  int nameChanged = 0;
  if (currentName.MakeTerminalNumbersMatch( bondedName )) {
    currentName.PrintName("Proposed name");
    std::string currentNameStr = currentName.NameString();
    if (currentNameStr.size() > (unsigned int)maxNameWidth_) {
      mprintf("Warning: Proposed name '%s' is larger than max name width %i\n", currentNameStr.c_str(), maxNameWidth_);
      //fixChainedAtoms(topIn, it->Idx(), currentAtom.Bond(0), anames);
    } else {
      NameType proposedName( currentName.NameString() );
      DuplicateMap::const_iterator jt = duplicateMap_.find( proposedName );
      if (jt == duplicateMap_.end()) {
        // The proposed name is free.
        mprintf("\t  Changing %s to match bonded atom number.\n", topIn.AtomMaskName(atomIdx).c_str());
        change_aname( topIn, atomIdx, proposedName );
        //keepName[it - atoms.begin()] = true;
        nameChanged = 1;
        //anames.insert( proposedName );
      } else
        mprintf("Warning: Proposed name '%s' is already in use.\n", currentNameStr.c_str());
    }
  }
  return nameChanged;
}

/** Given the index of an atom with hydrogens bonded to it, 
  * rename based off the current atom name. 
  */
int Exec_FixAtomNames::RenameBondedHydrogens(Topology& topIn, int heavyIdx)
{
  Iarray hydrogenAtoms;
  Atom const& heavyAtom = topIn[heavyIdx];
  for (Atom::bond_iterator bat = heavyAtom.bondbegin(); bat != heavyAtom.bondend(); ++bat)
  {
    if (topIn[*bat].Element() == Atom::HYDROGEN)
      hydrogenAtoms.push_back( *bat );
  }
  if (hydrogenAtoms.empty()) return 0;
  // Dissect the heavy atom name, determine if it has a prime
  AtomName heavyName;
  heavyName.DissectName( heavyAtom.Name() );
  heavyName.PrintName("Heavy atom name");
  std::string elementName(heavyAtom.ElementName());
  bool has_prime = heavyName.HasPrime();
  mprintf("\tElement name: '%s'\n", elementName.c_str());
  // Determine the suffix
  std::string suffix("");
  // For non-carbons, add the element name. TODO do for 2 char names too?
  if (heavyAtom.Element() != Atom::CARBON && elementName.size() == 1)
    suffix.append( elementName );
  int terminalNumber = heavyName.TerminalNumber();
  if (terminalNumber > -1) {
    std::string tnumStr = integerToString(terminalNumber);
    // If the digit is too big just omit the element char
    unsigned int maxterm = 3;
    if (has_prime)
      maxterm = 2;
    if (tnumStr.size() < maxterm)
      suffix.append( tnumStr );
    else
      suffix = tnumStr;
  } else {
    // Does the element name match?
    std::string tmpName = heavyAtom.Name().Truncated();
    std::size_t found = tmpName.find( elementName );
    if (found != std::string::npos && found < tmpName.size()) {
      suffix.append( tmpName.substr(found+1) );
      mprintf("DEBUG: heavy name '%s' element name '%s' found=%zu  substr='%s'\n", tmpName.c_str(), elementName.c_str(), found, suffix.c_str());
    }
  }
  if (has_prime)
    suffix.append("'");
  // Loop over all bonded hydrogens and rename them
  bool allNamesValid = false;
  typedef std::vector< std::pair<int,NameType> > IdxNamePairArray;
  IdxNamePairArray idxName;
  allNamesValid = true;
  int idx = 1;
  // PDB V2 convention is to put the index first, e.g. 1H5
  // PDB V3 convention is to put the index last, e.g. H51
  for (Iarray::const_iterator it = hydrogenAtoms.begin(); it != hydrogenAtoms.end(); ++it, idx++) {
    std::string newNameStr;
    if (hydrogenAtoms.size() == 1)
      newNameStr.assign( "H" + suffix );
    else if (hydrogenAtoms.size() == 2 && has_prime) {
      newNameStr.assign( "H" + suffix );
      suffix.append("'");
    } else
      newNameStr.assign( "H" + suffix + integerToString(idx) );
    if (newNameStr.size() > (unsigned int)maxNameWidth_) {
      mprintf("Warning: New hydrogen name %s is too big.\n", newNameStr.c_str());
      allNamesValid = false;
      break;
    } else {
      NameType proposedName(newNameStr);
      mprintf("\tProposed bonded H name %s\n", *proposedName);
      if (proposedName == topIn[*it].Name())
        mprintf("\t  Name is the same as existing name.\n");
      else {
        DuplicateMap::const_iterator jt = duplicateMap_.find( proposedName );
        if (jt != duplicateMap_.end()) {
          mprintf("Warning: Proposed name already in use.\n");
          allNamesValid = false;
          break;
        }
        idxName.push_back( std::pair<int,NameType>(*it, proposedName) );
      }
    }
  }
  // Actually change the names
  if (allNamesValid) {
    for (IdxNamePairArray::const_iterator kt = idxName.begin(); kt != idxName.end(); ++kt) {
      change_aname(topIn, kt->first, kt->second);
      //anames.insert( kt->second );
    }
  }
  
  if (allNamesValid)
    return 1;
  return 0;
}

/// Recursive function to trace an atom back through single bonds to a branched atom
static inline void trace_chain(Topology const& topIn, std::vector<int>& chainAtoms, int currentIdx, unsigned int idx, unsigned int max)
{
  if (idx >= max) return;

  chainAtoms.push_back( currentIdx );
  if (topIn[currentIdx].Nbonds() < 3) {
    for (Atom::bond_iterator bat = topIn[currentIdx].bondbegin(); bat != topIn[currentIdx].bondend(); ++bat)
    {
      bool in_chain = false;
      for (std::vector<int>::const_iterator it = chainAtoms.begin(); it != chainAtoms.end(); ++it)
      {
        if (*it == *bat) {
          in_chain = true;
          break;
        }
      }
      if (!in_chain) {
        trace_chain( topIn, chainAtoms, *bat, idx+1, max );
      }
    }
  }
}

/** Given the index of an atom with a duplicated name and the index of the
  * only other atom it is bonded to, trace the bonds back to the first atom
  * that has multiple bonds. Try to rename the entire chain based on that
  * original atom.
  * \return 1 if the chain was renamed, false otherwise
  */
int Exec_FixAtomNames::fixChainedAtoms(Topology& topIn, int duplicateIdx, int bondedIdx)
{
  Iarray chainAtoms;
  chainAtoms.push_back( duplicateIdx );
  //chainAtoms.push_back( bondedIdx );
  unsigned int idx = 0;
  Residue const& res = topIn.Res( topIn[duplicateIdx].ResNum() );
  trace_chain( topIn, chainAtoms, bondedIdx, idx, res.NumAtoms() );
  mprintf("DEBUG: Bond chain:");
  for (Iarray::const_iterator it = chainAtoms.begin(); it != chainAtoms.end(); it++)
    mprintf(" %s", topIn.AtomMaskName(*it).c_str());
  mprintf("\n");
  // Use last atom in this chain as the "origin". See if we can make the chain numbering consistent.
  Iarray::const_reverse_iterator rit = chainAtoms.rbegin();
  AtomName originName;
  originName.DissectName( topIn[*rit].Name() );
  originName.PrintName("Origin name");
  int terminalNumber = originName.TerminalNumber();
  bool allNamesValid = false;
  if (terminalNumber > -1) {
    typedef std::vector< std::pair<int,NameType> > IdxNamePairArray;
    IdxNamePairArray idxName;
    ++rit;
    allNamesValid = true;
    for (; rit != chainAtoms.rend(); ++rit) {
      AtomName chainName;
      chainName.DissectName( topIn[*rit].Name() );
      chainName.RemoveLeadingDigit();
      chainName.SetTerminalNumber( terminalNumber, maxNameWidth_ );
      chainName.PrintName("Proposed name");
      NameType proposedName(chainName.NameString());
      DuplicateMap::const_iterator jt = duplicateMap_.find( proposedName );
      if (jt != duplicateMap_.end()) {
        mprintf("Warning: Proposed name already in use.\n");
        allNamesValid = false;
        break;
      }
      idxName.push_back( std::pair<int,NameType>(*rit, proposedName) );
    }
    // Actually change the names
    if (allNamesValid) {
      for (IdxNamePairArray::const_iterator kt = idxName.begin(); kt != idxName.end(); ++kt)
        change_aname(topIn, kt->first, kt->second);
    }
  }

  if (allNamesValid)
    return 1;

  return 0;
}

/** Given an array of atom indices that share an atom name in a single
  * residue, determine one which is the "true" name and rename the
  * remaining ones with unique names.
  */
int Exec_FixAtomNames::fixDuplicateSet(Topology& topIn, Iarray const& Indexes, int rnum)
{
  // TODO rnum may not be needed anymore
  // Create a set of current atom names in the residue. They will all be unique.
  //Residue const& res = topIn.Res(rnum);
  //std::set<NameType> anames;
  //for (int iat = res.FirstAtom(); iat != res.LastAtom(); iat++)
  //  anames.insert( topIn[iat].Name() );

  // It is possible that a previous iteration has resolved these particular name clashes. Check here.
//  std::map<NameType, int> nameCount;
//  for (Iarray::const_iterator it = Indexes.begin(); it != Indexes.end(); ++it)
//  {
//    std::map<NameType, int>::iterator jt = nameCount.find( topIn[*it].Name() );
//    if (jt == nameCount.end())
//      nameCount.insert( std::pair<NameType, int>(topIn[*it].Name(), 1) );
//    else
//      jt->second++;
//  }

  // Loop over all atoms with the same name
  std::vector<DupAtom> atoms;
  for (Iarray::const_iterator it = Indexes.begin(); it != Indexes.end(); ++it)
  {
    Atom const& currentAtom = topIn[*it];
    DuplicateMap::const_iterator jt = duplicateMap_.find( currentAtom.Name() );
    //std::map<NameType, int>::const_iterator jt = nameCount.find( currentAtom.Name() );
    if (jt->second.size() < 2)
      mprintf("\tAtom %s no longer has a duplicated name.\n", topIn.AtomMaskName(*it).c_str());
    else {
      mprintf("Duplicated name: %s (%zu times)\n", topIn.AtomMaskName(*it).c_str(), jt->second.size());
      // Calculate the weight of all bonded atoms
      int wt = 0;
      for (Atom::bond_iterator bat = currentAtom.bondbegin(); bat != currentAtom.bondend(); ++bat)
        wt += topIn[*bat].AtomicNumber();
      atoms.push_back( DupAtom(*it, wt) );
    }
  }
  if (atoms.size() < 2)
    mprintf("\tNo more duplicates.\n");

  // ---------------------------------------------
  // If only bonded to one other atom, see if one of the atoms number
  // prefix (if present) matches its bonded atoms prefix.
  int nChainsRenamed = 0;
  for (std::vector<DupAtom>::const_iterator it = atoms.begin(); it != atoms.end(); ++it)
  {
    Atom const& currentAtom = topIn[it->Idx()];
    if (currentAtom.Nbonds() == 1) {
      Atom const& bondedAtom = topIn[currentAtom.Bond(0)];
      if (bondedAtom.Nbonds() < 3) {
        if (fixChainedAtoms(topIn, it->Idx(), currentAtom.Bond(0)))
          nChainsRenamed++;
      }
    }
  }
  if (nChainsRenamed > 0) {
    mprintf("\t%i atom chains were renamed.\n", nChainsRenamed);
    // Do we still have name clashes?
//    std::vector<DupAtom> newAtoms = getDuplicatedAtoms(atoms, topIn, anames);
    std::vector<DupAtom> newAtoms;
    for (std::vector<DupAtom>::const_iterator it = atoms.begin(); it != atoms.end(); ++it)
    {
      DuplicateMap::const_iterator jt = duplicateMap_.find( topIn[it->Idx()].Name() );
      if (jt->second.size() > 1)
        // Still duplicated
        newAtoms.push_back( *it );
    }
    atoms = newAtoms;
  }
  if (atoms.size() < 2) {
    mprintf("\tAll name clashes resolved.\n");
    return 0;
  }

  // ---------------------------------------------
  // Special case - see if these are hydrogens bonded to a common atom.
  if ( topIn[atoms.front().Idx()].Element() == Atom::HYDROGEN ) {
    int bondedIdx = -1;
    for (std::vector<DupAtom>::const_iterator it = atoms.begin(); it != atoms.end(); ++it)
    {
      Atom const& currentAtom = topIn[it->Idx()];
      if ( currentAtom.Element() != Atom::HYDROGEN ) {
        bondedIdx = -1;
        break;
      } else if ( currentAtom.Nbonds() > 1 ) {
        bondedIdx = -1;
        break;
      } else {
        if (bondedIdx == -1)
          bondedIdx = currentAtom.Bond(0);
        else if (bondedIdx != currentAtom.Bond(0)) {
          bondedIdx = -1;
          break;
        }
      }
    } // END loop over duplicated atoms
    if (bondedIdx != -1) {
      mprintf("\tDuplicated atoms are bonded to %s\n", topIn.AtomMaskName(bondedIdx).c_str());
      if (RenameBondedHydrogens(topIn, bondedIdx)) {
        mprintf("\tBonded hydrogens have been renamed.\n");
        return 0;
      }
    }
  }

  // ---------------------------------------------
  // Try to prioritize one atom as the "true" atom, i.e. the atom
  // whose name should not be changed.
  // Sort by 1) weight of all bonded atoms 2) topology index
  std::sort( atoms.begin(), atoms.end() );

  // Do we have a "true" atom?
  std::vector<bool> keepName(atoms.size(), false);
  bool have_trueAtom = (atoms[0].Sum() > atoms[1].Sum());
  if (have_trueAtom)
    keepName[0] = true;

  // If only bonded to one other atom, see if one of the atoms number
  // prefix (if present) matches its bonded atoms prefix.
  long int numMatchIdx = -1;
  for (std::vector<DupAtom>::const_iterator it = atoms.begin(); it != atoms.end(); ++it)
  {
    Atom const& currentAtom = topIn[it->Idx()];
    if (currentAtom.Nbonds() == 1) {
      AtomName currentName, bondedName;
      if (currentName.DissectName( currentAtom.Name() )) return 1;
      if (bondedName.DissectName( topIn[currentAtom.Bond(0)].Name() )) return 1;
      currentName.PrintName("Current Name");
      bondedName.PrintName("Bonded Name");
      if (currentName.TerminalNumbersMatch( bondedName )) {
        // Terminal numbers match. If we do not already have a "true" atom,
        // make it this one.
        if (!have_trueAtom) {
          mprintf("DEBUG: Terminal numbers match.\n");
          if (numMatchIdx != -1) {
            // Both atoms have terminal numbers that match their bonded partners.
            numMatchIdx = -1;
            break;
          } else {
            numMatchIdx = (it - atoms.begin());
          }
        }
      } else {
        // Can we resolve the clash by making this number match the bonded atom number?
        if (MakeTerminalNumberMatch(topIn, it->Idx(), currentAtom.Bond(0), currentName, bondedName)) {
          keepName[it - atoms.begin()] = true;
        }
      }
    } // END atom has one bond
  } // END loop over duplicated atoms
  if (numMatchIdx != -1) {
    keepName[numMatchIdx] = true;
    have_trueAtom = true;
  }

  // If we still do not have a "true" atom, just use the first one in the array
  if (!have_trueAtom)
    keepName[0] = true;

  // Increment the name of each non-"true" atom
  unsigned int idx = 0;
  for (std::vector<DupAtom>::const_iterator it = atoms.begin(); it != atoms.end(); ++it, idx++)
  {
    mprintf("\t  Atom %s (priority=%i) (keepName=%i)\n", topIn.AtomMaskName(it->Idx()).c_str(), it->Sum(), (int)keepName[idx]);
    //for (Atom::bond_iterator bat = topIn[it->Idx()].bondbegin(); bat != topIn[it->Idx()].bondend(); ++bat)
    //{
    //  mprintf("\t\tbonded to %s\n", topIn.AtomMaskName(*bat).c_str());
    //}
    if (!keepName[idx]) {
      if (fixDuplicateByIncrement(topIn, it->Idx())) {
        mprinterr("Error: Could not fix atom %s by name increment.\n", topIn.AtomMaskName(it->Idx()).c_str());
        return 1;
      }
      // If this was not a hydrogen, it may have bonded hydrogens that need to be renamed.
      if (topIn[it->Idx()].Element() != Atom::HYDROGEN) {
        RenameBondedHydrogens(topIn, it->Idx());
      }
    }
  }
  return 0;
}

/** Given the index of an atom that has a duplicated name in a residue, change
  * it to a unique name in that residue by incrementing any terminal number.
  */
int Exec_FixAtomNames::fixDuplicateByIncrement(Topology& topIn, int duplicateIdx)
{
  // Dissect this name.
  AtomName duplicateName;
  if (duplicateName.DissectName( topIn[duplicateIdx].Name() )) return 1;
  // DEBUG
  duplicateName.PrintName("Name to increment");

  // Try to increment the name until we do not get a clash
  bool tryNewNames = true;
  while (tryNewNames) {
    if (duplicateName.IncrementName(maxNameWidth_) != 0) {
      mprinterr("Error: Ran out of names to increment for %s.\n", topIn.AtomMaskName(duplicateIdx).c_str());
      return 1;
    }
    NameType newName( duplicateName.NameString() );
    DuplicateMap::const_iterator it = duplicateMap_.find( newName );
    if (it == duplicateMap_.end()) {
      // Name is no longer duplicated!
      mprintf("\t New name: %s\n", *newName);
      change_aname(topIn, duplicateIdx, newName);
      tryNewNames = false;
      //anames.insert( newName );
    } else
      mprintf("\tName %s already exists.\n", *newName);
  }

  return 0;
}
