#include "Exec_PrepareForLeap.h"
#include "CpptrajStdio.h"
#include "DistRoutines.h"
#include "CharMask.h"
#include "TorsionRoutines.h"
#include "Constants.h"
#include "CpptrajFile.h"
#include <set>
#include <stack>

// Exec_PrepareForLeap::Help()
void Exec_PrepareForLeap::Help() const
{
  mprintf("\tcrdset <coords set> [frame <#>] out <file>\n"
          "\t[cysmask <mask>] [disulfidecut <cut>] [newcysname <name>]\n"
          "\t[sugarmask <mask>]\n"
          "\t[leapunitname <name>]\n"
         );
}

/// Used to change residue name to nameIn
static inline void ChangeResName(Residue& res, NameType const& nameIn) {
  if (res.Name() != nameIn) {
    mprintf("\tChanging residue %s to %s\n", *(res.Name()), *nameIn);
    res.SetName( nameIn );
  }
}

/// Used to change atom name to nameIn
static inline void ChangeAtomName(Atom& atm, NameType const& nameIn) {
  if (atm.Name() != nameIn) {
    mprintf("\tChanging atom %s to %s\n", *(atm.Name()), *nameIn);
    atm.SetName( nameIn );
  }
}

/// Generate leap bond command for given atoms
void Exec_PrepareForLeap::LeapBond(int at1, int at2, Topology const& topIn, CpptrajFile* outfile)
const
{
  outfile->Printf("bond %s.%i.%s %s.%i.%s\n",
                  leapunitname_.c_str(), topIn[at1].ResNum()+1, *(topIn[at1].Name()),
                  leapunitname_.c_str(), topIn[at2].ResNum()+1, *(topIn[at2].Name()));
}

/// \return Glycam linkage code for given glycam residue name and linked atoms
static std::string LinkageCode(char glycamChar, std::set<NameType> const& linkages)
{
  std::string linkcode;
  // NOTE: I dont know if this is the best way to ID linkages. Seems easiest though.
  std::string linkstr;
  for (std::set<NameType>::const_iterator it = linkages.begin(); it != linkages.end(); ++it)
    linkstr.append( it->Truncated() );
  mprintf("\t  linkstr= '%s'\n", linkstr.c_str());
  switch (glycamChar) {
    case 'Y':
      if      (linkstr == "C1") linkcode = "0";
      else if (linkstr == "C1O4") linkcode = "4";
      break;
    default:
      mprintf("Warning: Unrecognized glycam residue char: %c\n", glycamChar);
  }
  if (linkcode.empty())
    mprintf("Warning: Could not determine link code.\n");
  return linkcode;
}

/** Attempt to identify sugar residue, form, and linkages. */
int Exec_PrepareForLeap::IdentifySugar(int rnum, Topology* topIn,
                                       Frame const& frameIn, CharMask const& cmask,
                                       CpptrajFile* outfile)
const
{
  Residue& res = topIn->SetRes(rnum);
  // Try to ID the base sugar type from the input name.
  char resChar = ' ';
  if (res.Name() == "NAG") {
    resChar = 'Y';
  } else {
    mprintf("Warning: Could not identify sugar from residue name '%s'\n", *res.Name());
    return 1;
  }
  mprintf("\tSugar %s glycam name: %c\n", *res.Name(), resChar);
  // Try to identify the form
  /* The alpha form has the CH2OH substituent (C5-C6 etc in Glycam) on the 
   * opposite side of the OH on the anomeric carbon (C1 in Glycam), while
   * in the beta form it is on the same side.
   */ 
  std::string formStr; 
  int C6idx = -1;
  int C5idx = -1;
  int C1idx = -1;
  int Cxidx = -1; // Non-H1, O5, C2 substituent of C1
  // Use a set to store linkages so they are in alphabetical order for easier identification
  std::set<NameType> linkages;
  // Bonds to non sugars to be removed since these will confuse tleap
  BondArray bondsToRemove;
  // Loop over sugar atoms
  for (int at = res.FirstAtom(); at != res.LastAtom(); at++)
  {
    // Rename some common atoms TODO need to check glycam residue char?
    if ( (*topIn)[at].Name() == "C7" )
      ChangeAtomName(topIn->SetAtom(at), "C2N");
    else if ( (*topIn)[at].Name() == "O7" )
      ChangeAtomName(topIn->SetAtom(at), "O2N");
    else if ( (*topIn)[at].Name() == "C8" )
      ChangeAtomName(topIn->SetAtom(at), "CME");
    // Identify atoms for torsion
    if ( (*topIn)[at].Name() == "C6" )
      C6idx = at;
    else if ( (*topIn)[at].Name() == "C5" )
      C5idx = at;
    else if ( (*topIn)[at].Name() == "C1" ) {
      C1idx = at;
      // Check substituent of C1
      for (Atom::bond_iterator bat = (*topIn)[at].bondbegin();
                               bat != (*topIn)[at].bondend(); ++bat)
      {
        if ( (*topIn)[*bat].Element() != Atom::HYDROGEN &&
             (*topIn)[*bat].Name() != "O5" &&
             (*topIn)[*bat].Name() != "C2" )
        {
          if (Cxidx != -1) {
            mprintf("Warning: Multiple substituents for C1: %s %s\n",
                    *(*topIn)[*bat].Name(), *(*topIn)[Cxidx].Name());
            return 1;
          } else
            Cxidx = *bat;
        }
      }
    }
    // Check for bonds to other residues
    for (Atom::bond_iterator bat = (*topIn)[at].bondbegin();
                             bat != (*topIn)[at].bondend(); ++bat)
    {
      if ((*topIn)[*bat].ResNum() != rnum) {
        if (!cmask.AtomInCharMask(*bat)) {
          mprintf("\tSugar %s bonded to non-sugar %s\n",
                  topIn->ResNameNumAtomNameNum(at).c_str(),
                  topIn->ResNameNumAtomNameNum(*bat).c_str());
          linkages.insert( (*topIn)[at].Name() );
          bondsToRemove.push_back( BondType(at, *bat, -1) );
          // Check if this is a recognized linkage TODO put in another file?
          Residue& pres = topIn->SetRes( (*topIn)[*bat].ResNum() );
          if ( pres.Name() == "SER" ) {
            ChangeResName( pres, "OLS" );
          } else if ( pres.Name() == "THR" ) {
            ChangeResName( pres, "OLT" );
          } else if ( pres.Name() == "HYP" ) {
            ChangeResName( pres, "OLP" );
          } else if ( pres.Name() == "ASN" ) {
            ChangeResName( pres, "NLN" );
          } else {
            mprintf("Warning: Unrecognized link residue %s, not modifying name.\n", *pres.Name());
          }
        } else {
          mprintf("\tSugar %s bonded to sugar %s\n",
                  topIn->ResNameNumAtomNameNum(at).c_str(),
                  topIn->ResNameNumAtomNameNum(*bat).c_str());
          linkages.insert( (*topIn)[at].Name() );
        }
      }
    } // END loop over bonded atoms
  } // END loop over residue atoms
  mprintf("\t  C5= %i C6= %i C1= %i Cx= %i\n", C6idx, C5idx, C1idx, Cxidx);
  if (C6idx == -1) { mprintf("Warning: C6 index not found.\n"); return 1; }
  if (C5idx == -1) { mprintf("Warning: C5 index not found.\n"); return 1; }
  if (C1idx == -1) { mprintf("Warning: C1 index not found.\n"); return 1; }
  if (Cxidx == -1) { mprintf("Warning: Cx index not found.\n"); return 1; }
  double torsion = Torsion( frameIn.XYZ(C6idx), frameIn.XYZ(C5idx),
                            frameIn.XYZ(C1idx), frameIn.XYZ(Cxidx) );
  mprintf("\t  Torsion= %f deg\n", torsion * Constants::RADDEG);
  if (torsion < Constants::PI && torsion > -Constants::PI) {
    mprintf("\t  Beta form\n");
    formStr = "B";
  } else {
    mprintf("\t  Alpha form\n");
    formStr = "A";
  }
  // Determine linkage
  mprintf("\t  Link atoms:");
  for (std::set<NameType>::const_iterator it = linkages.begin();
                                          it != linkages.end(); ++it)
    mprintf(" %4s", *(*it));
  mprintf("\n");
  std::string linkcode = LinkageCode(resChar, linkages);
  mprintf("\t  Linkage code: %s\n", linkcode.c_str());
  if (linkcode.empty()) return 1;
  // Remove bonds to non-sugar
  for (BondArray::const_iterator bnd = bondsToRemove.begin();
                                 bnd != bondsToRemove.end(); ++bnd)
  {
    LeapBond(bnd->A1(), bnd->A2(), *topIn, outfile);
    topIn->RemoveBond(bnd->A1(), bnd->A2());
  }
  // Set new residue name
  NameType newResName( linkcode + std::string(1,resChar) + formStr );
  mprintf("\t  Glycam resname: %s\n", *newResName);
  ChangeResName(res, newResName);
  return 0;
}

/** Determine where molecules end based on connectivity. */
int Exec_PrepareForLeap::FindTerByBonds(Topology* topIn, CharMask const& maskIn)
const
{
  // NOTE: this code is the same algorithm from Topology::NonrecursiveMolSearch
  // TODO use a common molecule search backend
  std::stack<unsigned int> nextAtomToSearch;
  bool unassignedAtomsRemain = true;
  unsigned int currentAtom = 0;
  unsigned int currentMol = 0;
  unsigned int lowestUnassignedAtom = 0;
  std::vector<int> atomMolNum( topIn->Natom(), -1 );
  while (unassignedAtomsRemain) {
    // This atom is in molecule.
    atomMolNum[currentAtom] = currentMol;
    //mprintf("DEBUG:\tAssigned atom %u to mol %u\n", currentAtom, currentMol);
    // All atoms bonded to this one are in molecule.
    for (Atom::bond_iterator batom = (*topIn)[currentAtom].bondbegin();
                             batom != (*topIn)[currentAtom].bondend(); ++batom)
    {
      if (atomMolNum[*batom] == -1) { // -1 is no molecule
        if ((*topIn)[*batom].Nbonds() > 1)
          // Bonded atom has more than 1 bond; needs to be searched.
          nextAtomToSearch.push( *batom );
        else {
          // Bonded atom only bonded to current atom. No more search needed.
          atomMolNum[*batom] = currentMol;
          //mprintf("DEBUG:\t\tAssigned bonded atom %i to mol %u\n", *batom, currentMol);
        }
      }
    }
    if (nextAtomToSearch.empty()) {
      //mprintf("DEBUG:\tNo atoms left in stack. Searching for next unmarked atom.\n");
      // No more atoms to search. Find next unmarked atom.
      currentMol++;
      unsigned int idx = lowestUnassignedAtom;
      for (; idx != atomMolNum.size(); idx++)
        if (atomMolNum[idx] == -1) break;
      if (idx == atomMolNum.size())
        unassignedAtomsRemain = false;
      else {
        currentAtom = idx;
        lowestUnassignedAtom = idx + 1;
      }
    } else {
      currentAtom = nextAtomToSearch.top();
      nextAtomToSearch.pop();
      //mprintf("DEBUG:\tNext atom from stack: %u\n", currentAtom);
    }
  }
  //t_nostack.Stop();
  //t_nostack.WriteTiming(1, "Non-recursive mol search:");
  //return (int)currentMol;
  // For each selected atom, find last atom in corresponding molecule,
  // set corresponding residue as TER.
  int at = 0;
  while (at < topIn->Natom()) {
    // Find the next selected atom
    while (at < topIn->Natom() && !maskIn.AtomInCharMask(at)) at++;
    if (at < topIn->Natom()) {
      int currentMol = atomMolNum[at];
      // Seek to end of molecule
      while (at < topIn->Natom() && currentMol == atomMolNum[at]) at++;
      // The previous atom is the end
      int lastRes = (*topIn)[at-1].ResNum();
      mprintf("\tSetting residue %s as terminal.\n",
              topIn->TruncResNameNum(lastRes).c_str());
      topIn->SetRes(lastRes).SetTerminal( true );
    }
  }
  return 0;
}

// Exec_PrepareForLeap::Execute()
Exec::RetType Exec_PrepareForLeap::Execute(CpptrajState& State, ArgList& argIn)
{
  std::string crdset = argIn.GetStringKey("crdset");
  if (crdset.empty()) {
    mprinterr("Error: Must specify COORDS set with 'crdset'\n");
    return CpptrajState::ERR;
  }
  DataSet* ds = State.DSL().FindSetOfGroup( crdset, DataSet::COORDINATES );
  if (ds == 0) {
    mprinterr("Error: No COORDS set found matching %s\n", crdset.c_str());
    return CpptrajState::ERR;
  }
  DataSet_Coords& coords = static_cast<DataSet_Coords&>( *((DataSet_Coords*)ds) );
  int tgtframe = argIn.getKeyInt("frame", 1) - 1;
  mprintf("\tUsing frame %i from COORDS set %s\n", tgtframe+1, coords.legend());
  if (tgtframe < 0 || tgtframe >= (int)coords.Size()) {
    mprinterr("Error: Frame is out of range.\n");
    return CpptrajState::ERR;
  }
  Frame frameIn = coords.AllocateFrame();
  coords.GetFrame(tgtframe, frameIn);

  leapunitname_ = argIn.GetStringKey("leapunitname", "m");
  mprintf("\tUsing leap unit name: %s\n", leapunitname_.c_str());

  AtomMask sugarMask;
  std::string sugarmaskstr = argIn.GetStringKey("sugarmask");
  if (!sugarmaskstr.empty()) {
    if (sugarMask.SetMaskString(sugarmaskstr)) {
      mprinterr("Error: Setting sugar mask string.\n");
      return CpptrajState::ERR;
    }
  }
  // Get masks for molecules now since topology may be modified later.
  std::vector<AtomMask> molMasks;
  std::string mstr = argIn.GetStringKey("molmask");
  while (!mstr.empty()) {
    mprintf("\tAll atoms selected by '%s' will be in same molecule.\n", mstr.c_str());
    molMasks.push_back( AtomMask() );
    if (molMasks.back().SetMaskString( mstr )) {
      mprinterr("Error: Invalid mask.\n");
      return CpptrajState::ERR;
    }
    if (coords.Top().SetupIntegerMask( molMasks.back() )) return CpptrajState::ERR;
    molMasks.back().MaskInfo();
    if (molMasks.back().None()) {
      mprinterr("Error: Nothing selected by mask.\n");
      return CpptrajState::ERR;
    }
    mstr = argIn.GetStringKey("molmask");
  }
  CharMask determineMolMask;
  mstr = argIn.GetStringKey("determinemolmask");
  if (!mstr.empty()) {
    mprintf("\tAtoms in mask '%s' will determine molecules by bonds.\n", mstr.c_str());
    if (determineMolMask.SetMaskString(mstr)) {
      mprinterr("Error: Invalid mask.\n");
      return CpptrajState::ERR;
    }
    if (coords.Top().SetupCharMask( determineMolMask )) return CpptrajState::ERR;
    determineMolMask.MaskInfo();
    if (determineMolMask.None()) {
      mprinterr("Error: Nothing selected by mask.\n");
      return CpptrajState::ERR;
    }
  }

  CpptrajFile* outfile = State.DFL().AddCpptrajFile(argIn.GetStringKey("out"),
                                                    "LEaP Input", DataFileList::TEXT, true);
  if (outfile == 0) return CpptrajState::ERR;

  // Disulfide search
  double disulfidecut = argIn.getKeyDouble("disulfidecut", 2.5);
  std::string newcysnamestr = argIn.GetStringKey("newcysname", "CYX");
  NameType newcysname(newcysnamestr);
  mprintf("\tCysteine residues involved in disulfide bonds will be changed to: %s\n", *newcysname);
  std::string cysmaskstr = argIn.GetStringKey("cysmask", ":CYS@SG");
  mprintf("\tSearching for disulfide bonds with a cutoff of %g Ang.\n", disulfidecut);
  AtomMask cysmask;
  if (cysmask.SetMaskString( cysmaskstr )) {
    mprinterr("Error: Could not set up CYS mask string %s\n", cysmaskstr.c_str());
    return CpptrajState::ERR;
  }
  if (coords.Top().SetupIntegerMask( cysmask )) return CpptrajState::ERR;
  cysmask.MaskInfo();
  if (cysmask.None())
    mprintf("Warning: No cysteine sulfur atoms selected by %s\n", cysmaskstr.c_str());
  else {
    int nDisulfides = 0;
    double cut2 = disulfidecut * disulfidecut;
    // Try to find potential disulfide sites.
    for (AtomMask::const_iterator at1 = cysmask.begin(); at1 != cysmask.end(); ++at1) {
      for (AtomMask::const_iterator at2 = at1 + 1; at2 != cysmask.end(); ++at2) {
        bool isBonded = false;
        // Check if the bond already exists
        if (coords.Top()[*at1].IsBondedTo(*at2)) {
          mprintf("\tExisting disulfide: %s to %s\n",
                  coords.Top().ResNameNumAtomNameNum(*at1).c_str(),
                  coords.Top().ResNameNumAtomNameNum(*at2).c_str());
          isBonded = true;
        } else {
          // TODO imaging?
          double r2 = DIST2_NoImage(frameIn.XYZ(*at1), frameIn.XYZ(*at2));
          if (r2 < cut2) {
            mprintf("\tPotential disulfide: %s to %s (%g Ang.)\n",
                    coords.Top().ResNameNumAtomNameNum(*at1).c_str(),
                    coords.Top().ResNameNumAtomNameNum(*at2).c_str(), sqrt(r2));
            isBonded = true;
          }
        }
        if (isBonded) {
          nDisulfides++;
          LeapBond(*at1, *at2, coords.Top(), outfile);
          ChangeResName(coords.TopPtr()->SetRes(coords.Top()[*at1].ResNum()), newcysname);
          ChangeResName(coords.TopPtr()->SetRes(coords.Top()[*at2].ResNum()), newcysname);
        }
      }
    }
    mprintf("\tDetected %i disulfide bonds.\n", nDisulfides);
  }

  // Prepare sugars
  if (sugarMask.MaskStringSet()) {
    mprintf("\tPreparing sugars selected by '%s'\n", sugarMask.MaskString());
    if (coords.Top().SetupIntegerMask( sugarMask )) return CpptrajState::ERR;
    sugarMask.MaskInfo();
    if (sugarMask.None())
      mprintf("Warning: No sugar atoms selected by %s\n", sugarMask.MaskString());
    else {
      CharMask cmask( sugarMask.ConvertToCharMask(), sugarMask.Nselected() );
      // Get sugar residue numbers
      std::vector<int> sugarResNums = coords.Top().ResnumsSelectedBy( sugarMask );
      // For each sugar residue, see if it is bonded to a non-sugar residue.
      // If it is, remove that bond but record it.
      for (std::vector<int>::const_iterator rnum = sugarResNums.begin();
                                            rnum != sugarResNums.end(); ++rnum)
      {
        //Residue const& Res = coords.Top().Res(*rnum);
        // See if we recognize this sugar.
        if (IdentifySugar(*rnum, coords.TopPtr(), frameIn, cmask, outfile))
          return CpptrajState::ERR;
      } // END loop over sugar residues
    }
    // Bonds to sugars have been removed, so regenerate molecule info
    coords.TopPtr()->DetermineMolecules();
  }

  // Try to set terminal residues
  if (!molMasks.empty() || determineMolMask.MaskStringSet()) {
    // Reset terminal status
    for (int rnum = 0; rnum != coords.Top().Nres(); rnum++)
      coords.TopPtr()->SetRes(rnum).SetTerminal(false);
    // The final residue of each molMask is terminal
    for (std::vector<AtomMask>::const_iterator mask = molMasks.begin();
                                               mask != molMasks.end(); ++mask)
    {
      //std::vector<int> Rnums = coords.Top().ResnumsSelectedBy( *mask );
      int lastAtom = mask->back();
      int lastRes = coords.Top()[lastAtom].ResNum();
      mprintf("\tSetting residue %s as terminal.\n",
        coords.Top().TruncResNameNum(lastRes).c_str());
      coords.TopPtr()->SetRes(lastRes).SetTerminal( true );
    }
    // Set ter based on connectivity
    if (determineMolMask.MaskStringSet()) {
      if (FindTerByBonds(coords.TopPtr(), determineMolMask)) {
        mprinterr("Error: Could not set TER by connectivity.\n");
        return CpptrajState::ERR;
      }
    }
  }
      

  return CpptrajState::OK;
}
