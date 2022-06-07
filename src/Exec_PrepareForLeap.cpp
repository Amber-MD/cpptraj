#include "Exec_PrepareForLeap.h"
#include "CharMask.h"
#include "CpptrajStdio.h"
#include "DistRoutines.h" // SearchForDisulfides
#include "Structure/ResStatArray.h"
#include "Structure/SugarBuilder.h"
#include <algorithm> // sort
#include <stack> // FindTerByBonds

using namespace Cpptraj::Structure;

/** CONSTRUCTOR */
Exec_PrepareForLeap::Exec_PrepareForLeap() : Exec(COORDS),
  errorsAreFatal_(true),
  debug_(0)
{
  SetHidden(false);
}

/** If file not present, use a default set of residue names. */
void Exec_PrepareForLeap::SetPdbResNames() {
  //Protein
  pdb_res_names_.insert("ACE");
  pdb_res_names_.insert("ALA");
  pdb_res_names_.insert("ARG");
  pdb_res_names_.insert("ASH");
  pdb_res_names_.insert("ASN");
  pdb_res_names_.insert("ASP");
  pdb_res_names_.insert("CYM");
  pdb_res_names_.insert("CYS");
  pdb_res_names_.insert("CYX");
  pdb_res_names_.insert("GLH");
  pdb_res_names_.insert("GLN");
  pdb_res_names_.insert("GLU");
  pdb_res_names_.insert("GLY");
  pdb_res_names_.insert("HIE");
  pdb_res_names_.insert("HIP");
  pdb_res_names_.insert("HIS");
  pdb_res_names_.insert("HYP"); // Recognized by Glycam
  pdb_res_names_.insert("ILE");
  pdb_res_names_.insert("LEU");
  pdb_res_names_.insert("LYN");
  pdb_res_names_.insert("LYS");
  pdb_res_names_.insert("MET");
  pdb_res_names_.insert("NME");
  pdb_res_names_.insert("PHE");
  pdb_res_names_.insert("PRO");
  pdb_res_names_.insert("SER");
  pdb_res_names_.insert("THR");
  pdb_res_names_.insert("TRP");
  pdb_res_names_.insert("TYR");
  pdb_res_names_.insert("VAL");
  // DNA
  pdb_res_names_.insert("DA");
  pdb_res_names_.insert("DC");
  pdb_res_names_.insert("DG");
  pdb_res_names_.insert("DT");
  // RNA
  pdb_res_names_.insert("A");
  pdb_res_names_.insert("C");
  pdb_res_names_.insert("G");
  pdb_res_names_.insert("U");
}

/** Load PDB residue names recognized by Amber FFs from file. */
int Exec_PrepareForLeap::LoadPdbResNames(std::string const& fnameIn)
{
  std::string fname = fnameIn;
  if (fnameIn.empty()) {
    // Check CPPTRAJHOME
    const char* env = getenv("CPPTRAJHOME");
    if (env != 0) {
      fname.assign(env);
      fname += "/dat/PDB_ResidueNames.txt";
      mprintf("Info: Parameter file path from CPPTRAJHOME variable: '%s'\n", fname.c_str());
    } else {
      // Check AMBERHOME
      env = getenv("AMBERHOME");
      if (env != 0) {
        fname.assign(env);
        fname += "/AmberTools/src/cpptraj/dat/PDB_ResidueNames.txt";
        mprintf("Info: Parameter file path from AMBERHOME variable: '%s'\n", fname.c_str());
      }
    }
  }
  if (fname.empty()) {
    mprintf("Warning: No PDB residue name file specified and/or CPPTRAJHOME not set.\n"
            "Warning: Using standard set of PDB residue names.\n");
    SetPdbResNames();
    return 0;
  }
  mprintf("\tReading PDB residue names from '%s'\n", fname.c_str());

  CpptrajFile infile;
  if (infile.OpenRead(fname)) {
    mprinterr("Error: Could not open PDB residue name file.\n");
    return 1;
  }
  const char* ptr = 0;
  while ( (ptr = infile.NextLine()) != 0 ) {
    ArgList argline( ptr, " " );
    if (argline.Nargs() > 0) {
      if (argline[0][0] != '#') {
        pdb_res_names_.insert( argline[0] );
      }
    }
  }
  infile.CloseFile();

  return 0;
}

/** \return True if residue name is in pdb_to_glycam_ or pdb_res_names_,
  *               or is solvent.
  */
bool Exec_PrepareForLeap::IsRecognizedPdbRes(NameType const& rname,
                                             SugarBuilder const& sugarBuilder)
const
{
  if (sugarBuilder.IsRecognizedPdbSugar(rname))
    return true;
  SetType::const_iterator amberIt = pdb_res_names_.find( rname );
  if (amberIt != pdb_res_names_.end())
    return true;
  if (rname == solventResName_)
    return true;
  return false;
}

/** \return Array of residue numbers with unrecognized PDB res names. */
Exec_PrepareForLeap::Iarray
  Exec_PrepareForLeap::GetUnrecognizedPdbResidues(Topology const& topIn,
                                                  SugarBuilder const& sugarBuilder)
const
{
  Iarray rnums;
  for (int ires = 0; ires != topIn.Nres(); ires++)
  {
    if (!IsRecognizedPdbRes( topIn.Res(ires).Name(), sugarBuilder ))
    {
      mprintf("\t%s is unrecognized.\n", topIn.TruncResNameOnumId(ires).c_str());
      rnums.push_back( ires );
    }
  }
  return rnums;
}

/** Given an array of residue numbers with unrecognized PDB res names,
  * generate an array with true for unrecognized residues that are
  * either isolated or only bound to other unrecognized residues.
  */
Exec_PrepareForLeap::Iarray
  Exec_PrepareForLeap::GetIsolatedUnrecognizedResidues(Topology const& topIn,
                                                       Iarray const& rnums)
const
{
  typedef std::vector<bool> Barray;
  Barray isRecognized(topIn.Nres(), true);
  for (Iarray::const_iterator it = rnums.begin(); it != rnums.end(); ++it)
    isRecognized[ *it ] = false;

  Iarray isolated;
  for (Iarray::const_iterator it = rnums.begin(); it != rnums.end(); ++it)
  {
    bool isIsolated = true;
    Residue const& res = topIn.Res( *it );
    for (int at = res.FirstAtom(); at != res.LastAtom(); ++at)
    {
      for (Atom::bond_iterator bat = topIn[at].bondbegin(); bat != topIn[at].bondend(); ++bat)
      {
        if (topIn[*bat].ResNum() != *it)
        {
          // This bonded atom is in another residue. Is that residue recognized?
          if ( isRecognized[ topIn[*bat].ResNum() ] ) {
            // Residue *it is bonded to a recognized residue. Not isolated.
            isIsolated = false;
            break;
          }
        }
      } // END loop over residue atoms bonded atoms
      if (!isIsolated) break;
    } // END loop over residue atoms
    if (isIsolated) {
      mprintf("\t%s is isolated and unrecognized.\n", topIn.TruncResNameOnumId(*it).c_str());
      isolated.push_back( *it );
    }
  } // END loop over unrecognized residues

  return isolated;
}

// -----------------------------------------------------------------------------
/** Determine where molecules end based on connectivity. */
int Exec_PrepareForLeap::FindTerByBonds(Topology& topIn, CharMask const& maskIn)
const
{
  // NOTE: this code is the same algorithm from Topology::NonrecursiveMolSearch
  // TODO use a common molecule search backend
  std::stack<unsigned int> nextAtomToSearch;
  bool unassignedAtomsRemain = true;
  unsigned int currentAtom = 0;
  unsigned int currentMol = 0;
  unsigned int lowestUnassignedAtom = 0;
  Iarray atomMolNum( topIn.Natom(), -1 );
  while (unassignedAtomsRemain) {
    // This atom is in molecule.
    atomMolNum[currentAtom] = currentMol;
    //mprintf("DEBUG:\tAssigned atom %u to mol %u\n", currentAtom, currentMol);
    // All atoms bonded to this one are in molecule.
    for (Atom::bond_iterator batom = topIn[currentAtom].bondbegin();
                             batom != topIn[currentAtom].bondend(); ++batom)
    {
      if (atomMolNum[*batom] == -1) { // -1 is no molecule
        if (topIn[*batom].Nbonds() > 1)
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
  while (at < topIn.Natom()) {
    // Find the next selected atom
    while (at < topIn.Natom() && !maskIn.AtomInCharMask(at)) at++;
    if (at < topIn.Natom()) {
      int currentMol = atomMolNum[at];
      // Seek to end of molecule
      while (at < topIn.Natom() && currentMol == atomMolNum[at]) at++;
      // The previous atom is the end
      int lastRes = topIn[at-1].ResNum();
      mprintf("\tSetting residue %s as terminal.\n",
              topIn.TruncResNameOnumId(lastRes).c_str());
      topIn.SetRes(lastRes).SetTerminal( true );
    }
  }
  return 0;
}

// -----------------------------------------------------------------------------
/** Search for disulfide bonds. */
int Exec_PrepareForLeap::SearchForDisulfides(ResStatArray& resStat,
                                             double disulfidecut,
                                             std::string const& newcysnamestr,
                                             std::string const& cysmaskstr,
                                             bool searchForNewDisulfides,
                                             Topology& topIn, Frame const& frameIn,
                                             CpptrajFile* outfile)
{
  // Disulfide search
  NameType newcysname(newcysnamestr);
  mprintf("\tCysteine residues involved in disulfide bonds will be changed to: %s\n", *newcysname);
  if (searchForNewDisulfides)
    mprintf("\tSearching for disulfide bonds with a cutoff of %g Ang.\n", disulfidecut);
  else
    mprintf("\tOnly using existing disulfide bonds, will not search for new ones.\n");

  AtomMask cysmask;
  if (cysmask.SetMaskString( cysmaskstr )) {
    mprinterr("Error: Could not set up CYS mask string %s\n", cysmaskstr.c_str());
    return 1;
  }
  if (topIn.SetupIntegerMask( cysmask )) return 1; 
  cysmask.MaskInfo();
  if (cysmask.None())
    mprintf("Warning: No cysteine sulfur atoms selected by %s\n", cysmaskstr.c_str());
  else {
    // Sanity check - warn if non-sulfurs detected
    for (AtomMask::const_iterator at = cysmask.begin(); at != cysmask.end(); ++at)
    {
      if (topIn[*at].Element() != Atom::SULFUR)
        mprintf("Warning: Atom '%s' does not appear to be sulfur.\n",
                topIn.ResNameNumAtomNameNum(*at).c_str());
    }

    int nExistingDisulfides = 0;
    int nDisulfides = 0;
    double cut2 = disulfidecut * disulfidecut;
    // Try to find potential disulfide sites.
    // Keep track of which atoms will be part of disulfide bonds.
    Iarray disulfidePartner( cysmask.Nselected(), -1 );
    // First, check for existing disulfides.
    for (AtomMask::const_iterator at1 = cysmask.begin(); at1 != cysmask.end(); ++at1)
    {
      for (AtomMask::const_iterator at2 = at1 + 1; at2 != cysmask.end(); ++at2)
      {
        if (topIn[*at1].IsBondedTo(*at2)) {
          if (debug_ > 0)
            mprintf("\tExisting disulfide: %s to %s\n",
                    topIn.ResNameNumAtomNameNum(*at1).c_str(),
                    topIn.ResNameNumAtomNameNum(*at2).c_str());
          nExistingDisulfides++;
          int idx1 = (int)(at1 - cysmask.begin());
          int idx2 = (int)(at2 - cysmask.begin());
          disulfidePartner[idx1] = idx2;
          disulfidePartner[idx2] = idx1;
        }
      }
    }
    mprintf("\t%i existing disulfide bonds.\n", nExistingDisulfides);
    // DEBUG - Print current array
    if (debug_ > 1) {
      mprintf("DEBUG: Disulfide partner array after existing:\n");
      for (Iarray::const_iterator it = disulfidePartner.begin(); it != disulfidePartner.end(); ++it)
      {
        mprintf("  S %i [%li]", cysmask[it-disulfidePartner.begin()]+1, it-disulfidePartner.begin());
        if (*it == -1)
          mprintf(" None.\n");
        else
          mprintf(" to S %i [%i]\n", cysmask[*it]+1, *it);
      }
    }
    // Second, search for new disulfides from remaining sulfurs.
    if (searchForNewDisulfides) {
      // Only search with atoms that do not have an existing partner.
      Iarray s_idxs; // Indices into cysmask/disulfidePartner
      for (int idx = 0; idx != cysmask.Nselected(); idx++)
        if (disulfidePartner[idx] == -1)
          s_idxs.push_back( idx );
      mprintf("\t%zu sulfur atoms do not have a partner.\n", s_idxs.size());
      if (!s_idxs.empty()) {
        // In some structures, there may be 2 potential disulfide partners
        // within the cutoff. In that case, only choose the shortest.
        // To try to do this as directly as possible, calculate all possible S-S
        // distances, save the ones below the cutoff, sort them from shortest to
        // longest, then assign each that to a disulfide if not already assigned.
        // In this way, shorter S-S distances will be prioritized.
        typedef std::pair<int,int> IdxPair;
        typedef std::pair<double, IdxPair> D2Pair;
        typedef std::vector<D2Pair> D2Array;
        D2Array D2;

        for (Iarray::const_iterator it1 = s_idxs.begin(); it1 != s_idxs.end(); ++it1)
        {
          int at1 = cysmask[*it1];
          for (Iarray::const_iterator it2 = it1 + 1; it2 != s_idxs.end(); ++it2)
          {
            int at2 = cysmask[*it2];
            double r2 = DIST2_NoImage(frameIn.XYZ(at1), frameIn.XYZ(at2));
            if (r2 < cut2)
              D2.push_back( D2Pair(r2, IdxPair(*it1, *it2)) );
          }
        }
        std::sort(D2.begin(), D2.end());
        if (debug_ > 1) {
          mprintf("DEBUG: Sorted S-S array:\n");
          for (D2Array::const_iterator it = D2.begin(); it != D2.end(); ++it)
          {
            int at1 = cysmask[it->second.first];
            int at2 = cysmask[it->second.second];
            mprintf("  %8i - %8i = %g Ang.\n", at1+1, at2+2, sqrt(it->first));
          }
        }
        // All distances in D2 are below the cutoff
        for (D2Array::const_iterator it = D2.begin(); it != D2.end(); ++it)
        {
          if (disulfidePartner[it->second.first] == -1 &&
              disulfidePartner[it->second.second] == -1)
          {
            // Neither index has a partner yet
            int at1 = cysmask[it->second.first];
            int at2 = cysmask[it->second.second];
            mprintf("\t  Potential disulfide: %s to %s (%g Ang.)\n",
                    topIn.ResNameNumAtomNameNum(at1).c_str(),
                    topIn.ResNameNumAtomNameNum(at2).c_str(), sqrt(it->first));
            disulfidePartner[it->second.first ] = it->second.second;
            disulfidePartner[it->second.second] = it->second.first;
          }
        } // END loop over sorted distances
      } // END s_idxs not empty()
    } // END search for new disulfides
    // For each sulfur that has a disulfide partner, generate a bond command
    // and change residue name.
    for (Iarray::const_iterator idx1 = disulfidePartner.begin(); idx1 != disulfidePartner.end(); ++idx1)
    {
      if (*idx1 != -1) {
        int at1 = cysmask[idx1-disulfidePartner.begin()];
        int at2 = cysmask[*idx1];
        if (at1 < at2) {
          nDisulfides++;
          LeapBond(at1, at2, topIn, outfile);
        }
        ChangeResName(topIn.SetRes(topIn[at1].ResNum()), newcysname);
        resStat[topIn[at1].ResNum()] = VALIDATED;
      }
    }
    mprintf("\tDetected %i disulfide bonds.\n", nDisulfides);
  }
  return 0;
}

// -----------------------------------------------------------------------------
/** Try to determine histidine protonation from existing hydrogens.
  * Change residue names as appropriate.
  * \param topIn Input topology.
  * \param ND1 ND1 atom name.
  * \param NE2 NE2 atom name.
  * \param HisName PDB histidine name.
  * \param HieName Name for epsilon-protonated His.
  * \param HidName Name for delta-protonated His.
  * \param HipName Name for doubly-protonated His.
  */
int Exec_PrepareForLeap::DetermineHisProt( Topology& topIn,
                                           NameType const& ND1,
                                           NameType const& NE2,
                                           NameType const& HisName,
                                           NameType const& HieName,
                                           NameType const& HidName,
                                           NameType const& HipName )
const
{
  mprintf("\tAttempting to determine histidine form from any existing H atoms.\n");
  std::string hisMaskStr = ":" + HisName.Truncated();
  AtomMask mask;
  if (mask.SetMaskString( hisMaskStr )) {
    mprinterr("Error: Invalid His mask string: %s\n", hisMaskStr.c_str());
    return 1;
  }
  if ( topIn.SetupIntegerMask( mask )) return 1;
  mask.MaskInfo();
  Iarray resIdxs = topIn.ResnumsSelectedBy( mask );
  // Loop over selected histidine residues
  unsigned int nchanged = 0;
  for (Iarray::const_iterator rnum = resIdxs.begin();
                              rnum != resIdxs.end(); ++rnum)
  {
    if (debug_ > 1)
      mprintf("DEBUG: %s (%i) (%c)\n", topIn.TruncResNameOnumId(*rnum).c_str(), topIn.Res(*rnum).OriginalResNum(), topIn.Res(*rnum).ChainId());
    int nd1idx = -1;
    int ne2idx = -1;
    Residue const& hisRes = topIn.Res( *rnum );
    for (int at = hisRes.FirstAtom(); at < hisRes.LastAtom(); ++at)
    {
      if ( (topIn[at].Name() == ND1 ) )
        nd1idx = at;
      else if ( (topIn[at].Name() == NE2 ) )
        ne2idx = at;
    }
    if (nd1idx == -1) {
      mprintf("Warning: Atom %s not found for %s; skipping residue.\n", *ND1, topIn.TruncResNameOnumId(*rnum).c_str());
      continue;
    }
    if (ne2idx == -1) {
      mprintf("Warning: Atom %s not found for %s; skipping residue,\n", *NE2, topIn.TruncResNameOnumId(*rnum).c_str());
      continue;
    }
    if (debug_ > 1)
      mprintf("DEBUG: %s nd1idx= %i ne2idx= %i\n",
              topIn.TruncResNameOnumId( *rnum ).c_str(), nd1idx+1, ne2idx+1);
    // Check for H bonded to nd1/ne2
    int nd1h = 0;
    for (Atom::bond_iterator bat = topIn[nd1idx].bondbegin();
                             bat != topIn[nd1idx].bondend();
                           ++bat)
      if ( topIn[*bat].Element() == Atom::HYDROGEN)
        ++nd1h;
    if (nd1h > 1) {
      mprinterr("Error: More than 1 hydrogen bonded to %s\n",
                topIn.ResNameNumAtomNameNum(nd1idx).c_str());
      return 1;
    }
    int ne2h = 0;
    for (Atom::bond_iterator bat = topIn[ne2idx].bondbegin();
                             bat != topIn[ne2idx].bondend();
                           ++bat)
      if ( topIn[*bat].Element() == Atom::HYDROGEN)
        ++ne2h;
    if (ne2h > 1) {
      mprinterr("Error: More than 1 hydrogen bonded to %s\n",
                topIn.ResNameNumAtomNameNum(ne2idx).c_str());
      return 1;
    }
    if (nd1h > 0 && ne2h > 0) {
      mprintf("\t\t%s => %s\n", topIn.TruncResNameOnumId(*rnum).c_str(), *HipName);
      ChangeResName( topIn.SetRes(*rnum), HipName );
      nchanged++;
    } else if (nd1h > 0) {
      mprintf("\t\t%s => %s\n", topIn.TruncResNameOnumId(*rnum).c_str(), *HidName);
      ChangeResName( topIn.SetRes(*rnum), HidName );
      nchanged++;
    } else if (ne2h > 0) {
      mprintf("\t\t%s => %s\n", topIn.TruncResNameOnumId(*rnum).c_str(), *HieName);
      ChangeResName( topIn.SetRes(*rnum), HieName );
      nchanged++;
    }
    //else {
    //  // Default to epsilon
    //  mprintf("\tUsing default name '%s' for %s\n", *HieName, topIn.TruncResNameOnumId(*rnum).c_str());
    //  HisResNames.push_back( HieName );
    //}
  }
  if (nchanged == 0) 
    mprintf("\tNo histidine names were changed.\n");
  else
    mprintf("\t%u histidine names were changed.\n", nchanged);
  return 0;
}

// -----------------------------------------------------------------------------

/** Modify coords according to user wishes. */
int Exec_PrepareForLeap::ModifyCoords( Topology& topIn, Frame& frameIn,
                                       bool remove_water,
                                       std::string const& altLocStr, std::string const& stripMask,
                                       std::string const& waterMask,
                                       Iarray const& resnumsToRemove )
const
{
  // Create a mask denoting which atoms will be kept.
  std::vector<bool> atomsToKeep( topIn.Natom(), true );
  // Previously-determined array of residues to remove
  for (Iarray::const_iterator rnum = resnumsToRemove.begin();
                              rnum != resnumsToRemove.end(); ++rnum)
  {
    Residue const& res = topIn.Res( *rnum );
    mprintf("\tRemoving %s\n", topIn.TruncResNameOnumId( *rnum ).c_str());
    for (int at = res.FirstAtom(); at != res.LastAtom(); at++)
      atomsToKeep[at] = false;
  }
  // User-specified strip mask
  if (!stripMask.empty()) {
    AtomMask mask;
    if (mask.SetMaskString( stripMask )) {
      mprinterr("Error: Invalid mask string '%s'\n", stripMask.c_str());
      return 1;
    }
    if (topIn.SetupIntegerMask( mask )) return 1;
    mask.MaskInfo();
    if (!mask.None()) {
      for (AtomMask::const_iterator atm = mask.begin(); atm != mask.end(); ++atm)
        atomsToKeep[*atm] = false;
    }
  }
  if (remove_water) {
    // Do not use cpptraj definition of solvent in case we have e.g. bound HOH
    AtomMask mask;
    if (mask.SetMaskString( waterMask )) {
      mprinterr("Error: Invalid solvent mask string '%s'\n", waterMask.c_str());
      return 1;
    }
    if (topIn.SetupIntegerMask( mask )) return 1;
    mask.MaskInfo();
    if (!mask.None()) {
      for (AtomMask::const_iterator atm = mask.begin(); atm != mask.end(); ++atm)
        atomsToKeep[*atm] = false;
    }

  }
  // Identify alternate atom location groups.
  if (!altLocStr.empty()) {
    if (topIn.AtomAltLoc().empty()) {
      mprintf("\tNo alternate atom locations.\n");
    } else {
      // Map atom name to atom indices
      typedef std::map<NameType, std::vector<int>> AlocMapType;
      AlocMapType alocMap;
      for (int rnum = 0; rnum != topIn.Nres(); rnum++) {
        alocMap.clear();
        for (int at = topIn.Res(rnum).FirstAtom(); at != topIn.Res(rnum).LastAtom(); at++) {
          if (topIn.AtomAltLoc()[at] != ' ') {
            AlocMapType::iterator it = alocMap.find( topIn[at].Name() );
            if (it == alocMap.end()) {
              alocMap.insert( std::pair<NameType, std::vector<int>>( topIn[at].Name(),
                                                                     std::vector<int>(1, at) ));
            } else {
              it->second.push_back( at );
            }
          }
        } // END loop over atoms in residue
        if (!alocMap.empty()) {
          if (debug_ > 0)
            mprintf("DEBUG: Alternate loc. for %s\n", topIn.TruncResNameOnumId(rnum).c_str());
          // Loop over atoms with alternate locations
          for (AlocMapType::const_iterator it = alocMap.begin(); it != alocMap.end(); ++it) {
            if (debug_ > 0) {
              // Print all alternate atoms
              mprintf("\t'%s'", *(it->first));
              for (std::vector<int>::const_iterator at = it->second.begin();
                                                    at != it->second.end(); ++at)
                mprintf(" %s[%c]", *(topIn[*at].Name()), topIn.AtomAltLoc()[*at]);
              mprintf("\n");
            }
            // For each, choose which location to keep.
            if (altLocStr.size() == 1) {
              // Keep only specified character
              char altLocChar = altLocStr[0];
              for (std::vector<int>::const_iterator at = it->second.begin();
                                                    at != it->second.end(); ++at)
                if (topIn.AtomAltLoc()[*at] != altLocChar)
                  atomsToKeep[*at] = false;
            } else {
              // Keep highest occupancy
              if (topIn.Occupancy().empty()) {
                mprintf("\tNo occupancy.\n"); // TODO error?
              } else {
                int highestOccAt = -1;
                float highestOcc = 0;
                for (std::vector<int>::const_iterator at = it->second.begin();
                                                      at != it->second.end(); ++at)
                {
                  if (highestOccAt == -1) {
                    highestOccAt = *at;
                    highestOcc = topIn.Occupancy()[*at];
                  } else if (topIn.Occupancy()[*at] > highestOcc) {
                    highestOccAt = *at;
                    highestOcc = topIn.Occupancy()[*at];
                  }
                }
                // Set everything beside highest occ to false
                for (std::vector<int>::const_iterator at = it->second.begin();
                                                      at != it->second.end(); ++at)
                  if (*at != highestOccAt)
                    atomsToKeep[*at] = false;
              }
            }
          } // END loop over atoms with alternate locations
        }
      } // END loop over residue numbers
    }
  }

  // Set up mask of only kept atoms.
  AtomMask keptAtoms;
  keptAtoms.SetNatoms( topIn.Natom() );
  for (int idx = 0; idx != topIn.Natom(); idx++) {
    if (atomsToKeep[idx])
      keptAtoms.AddSelectedAtom(idx);
  }
  if (keptAtoms.Nselected() == topIn.Natom())
    // Keeping everything, no modifications
    return 0;
  // Modify top/frame
  Topology* newTop = topIn.modifyStateByMask( keptAtoms );
  if (newTop == 0) {
    mprinterr("Error: Could not create new topology.\n");
    return 1;
  }
  newTop->Brief("After removing atoms:");
  Frame newFrame;
  newFrame.SetupFrameV(newTop->Atoms(), frameIn.CoordsInfo());
  newFrame.SetFrame(frameIn, keptAtoms);

  topIn = *newTop;
  frameIn = newFrame;
  delete newTop;

  return 0;
}

/** Remove any hydrogen atoms. This is done separately from ModifyCoords
  * so that hydrogen atom info can be used to e.g. ID histidine
  * protonation.
  */
int Exec_PrepareForLeap::RemoveHydrogens(Topology& topIn, Frame& frameIn) const {
  AtomMask keptAtoms;
  keptAtoms.SetNatoms( topIn.Natom() );
  for (int idx = 0; idx != topIn.Natom(); idx++)
  {
    if (topIn[idx].Element() != Atom::HYDROGEN)
      keptAtoms.AddSelectedAtom(idx);
  }
  if (keptAtoms.Nselected() == topIn.Natom())
    // Keeping everything, no modification
    return 0;
  // Modify top/frame
  Topology* newTop = topIn.modifyStateByMask( keptAtoms );
  if (newTop == 0) {
    mprinterr("Error: Could not create new topology with no hydrogens.\n");
    return 1;
  }
  newTop->Brief("After removing hydrogen atoms:");
  Frame newFrame;
  newFrame.SetupFrameV(newTop->Atoms(), frameIn.CoordsInfo());
  newFrame.SetFrame(frameIn, keptAtoms);

  topIn = *newTop;
  frameIn = newFrame;
  delete newTop;

  return 0;
}

// -----------------------------------------------------------------------------

/// \return index of oxygen atom bonded to this atom but not in same residue
static inline int getLinkOxygenIdx(Topology const& leaptop, int at, int rnum) {
  int o_idx = -1;
  for (Atom::bond_iterator bat = leaptop[at].bondbegin();
                           bat != leaptop[at].bondend(); ++bat)
  {
    if (leaptop[*bat].Element() == Atom::OXYGEN && leaptop[*bat].ResNum() != rnum) {
      o_idx = *bat;
      break;
    }
  }
  return o_idx;
}

/// \return index of carbon bonded to link oxygen not in same residue
static inline int getLinkCarbonIdx(Topology const& leaptop, int at, int rnum)
{
  int o_idx = getLinkOxygenIdx(leaptop, at, rnum);
  if (o_idx == -1) return o_idx;
  int c_idx = -1;
  for (Atom::bond_iterator bat = leaptop[o_idx].bondbegin();
                           bat != leaptop[o_idx].bondend(); ++bat)
  {
    if (leaptop[*bat].Element() == Atom::CARBON && leaptop[*bat].ResNum() != rnum) {
      c_idx = *bat;
      break;
    }
  }
  return c_idx;
}

/** Run leap to generate topology. Modify the topology if needed. */
int Exec_PrepareForLeap::RunLeap(std::string const& ff_file,
                                 std::string const& leapfilename) const
{
  if (leapfilename.empty()) {
    mprintf("Warning: No leap input file name was specified, not running leap.\n");
    return 0;
  }
  if (ff_file.empty()) {
    mprintf("Warning: No leap input file with force fields was specified, not running leap.\n");
    return 0;
  }
  mprintf("\tExecuting leap.\n");

  std::string topname = leapunitname_ + ".parm7";
  std::string rstname = leapunitname_ + ".rst7";

  Cpptraj::LeapInterface LEAP(debug_);
  LEAP.AddInputFile( ff_file );
  LEAP.AddInputFile( leapfilename );
  LEAP.AddCommand("saveamberparm " + leapunitname_ + " " +
                  topname + " " + rstname);

  if (LEAP.RunLeap()) {
    mprinterr("Error: Leap failed.\n");
    return 1;
  }

  // Load the leap topology;
  Topology leaptop;
  ParmFile parm;
  if (parm.ReadTopology(leaptop, topname, debug_)) return 1;

  bool top_is_modified = false;
  // Go through each residue. Find ones that need to be adjusted.
  // NOTE: If deoxy carbons are ever handled, need to add H1 hydrogen and
  //       add the former -OH charge to the carbon.
  for (int rnum = 0; rnum != leaptop.Nres(); rnum++)
  {
    Residue const& res = leaptop.Res(rnum);
    if (res.Name() == "SO3") {
      int o_idx = -1;
      // Need to adjust the charge on the bonded oxygen by +0.031
      for (int at = res.FirstAtom(); at != res.LastAtom(); at++) {
        if (leaptop[at].Element() == Atom::SULFUR) {
          o_idx = getLinkOxygenIdx( leaptop, at, rnum );
          if (o_idx != -1) break;
        }
      }
      if (o_idx == -1) {
        mprinterr("Error: Could not find oxygen link atom for '%s'\n",
                  leaptop.TruncResNameOnumId(rnum).c_str());
        return 1;
      }
      double newcharge = leaptop[o_idx].Charge() + 0.031;
      mprintf("\tFxn group '%s'; changing charge on %s from %f to %f\n", *(res.Name()),
              leaptop.AtomMaskName(o_idx).c_str(), leaptop[o_idx].Charge(), newcharge);
      leaptop.SetAtom(o_idx).SetCharge( newcharge );
      top_is_modified = true;
    } else if (res.Name() == "MEX") {
      int c_idx = -1;
      // Need to adjust the charge on the carbon bonded to link oxygen by -0.039
      for (int at = res.FirstAtom(); at != res.LastAtom(); at++) {
        if (leaptop[at].Element() == Atom::CARBON) {
          c_idx = getLinkCarbonIdx( leaptop, at, rnum );
          if (c_idx != -1) break;
        }
      }
      if (c_idx == -1) {
        mprinterr("Error: Could not find carbon bonded to oxygen link atom for '%s'\n",
                  leaptop.TruncResNameOnumId(rnum).c_str());
        return 1;
      }
      double newcharge = leaptop[c_idx].Charge() - 0.039;
      mprintf("\tFxn group '%s'; changing charge on %s from %f to %f\n", *(res.Name()),
              leaptop.AtomMaskName(c_idx).c_str(), leaptop[c_idx].Charge(), newcharge);
      leaptop.SetAtom(c_idx).SetCharge( newcharge );
      top_is_modified = true;
    } else if (res.Name() == "ACX") {
      int c_idx = -1;
      // Need to adjust the charge on the carbon bonded to link oxygen by +0.008
      for (int at = res.FirstAtom(); at != res.LastAtom(); at++) {
        if (leaptop[at].Element() == Atom::CARBON) {
          // This needs to be the acetyl carbon, ensure it is bonded to an oxygen
          for (Atom::bond_iterator bat = leaptop[at].bondbegin();
                                   bat != leaptop[at].bondend(); ++bat)
          {
            if (leaptop[*bat].Element() == Atom::OXYGEN) {
              c_idx = getLinkCarbonIdx( leaptop, at, rnum );
              if (c_idx != -1) break;
            }
          }
          if (c_idx != -1) break;
        }
      }
      if (c_idx == -1) {
        mprinterr("Error: Could not find carbon bonded to oxygen link atom for '%s'\n",
                  leaptop.TruncResNameOnumId(rnum).c_str());
        return 1;
      }
      double newcharge = leaptop[c_idx].Charge() + 0.008;
      mprintf("\tFxn group '%s'; changing charge on %s from %f to %f\n", *(res.Name()),
              leaptop.AtomMaskName(c_idx).c_str(), leaptop[c_idx].Charge(), newcharge);
      leaptop.SetAtom(c_idx).SetCharge( newcharge );
      top_is_modified = true;
    }
  }

  // DEBUG: Print out total charge on each residue
  double total_q = 0;
  for (Topology::res_iterator res = leaptop.ResStart(); res != leaptop.ResEnd(); ++res)
  {
    double tcharge = 0;
    for (int at = res->FirstAtom(); at != res->LastAtom(); ++at) {
      total_q += leaptop[at].Charge();
      tcharge += leaptop[at].Charge();
    }
    if (debug_ > 0) {
      mprintf("DEBUG:\tResidue %10s charge= %12.5f\n",
              leaptop.TruncResNameOnumId(res-leaptop.ResStart()).c_str(), tcharge);
    }
  }
  mprintf("\tTotal charge: %16.8f\n", total_q);

  // If topology was modified, write it back out
  if (top_is_modified) {
    mprintf("\tWriting modified topology back to '%s'\n", topname.c_str());
    parm.WriteTopology(leaptop, topname, parm.CurrentFormat(), debug_);
  }

  return 0;
}

/** Print warnings for residues that will need to be modified in leap. */
void Exec_PrepareForLeap::LeapFxnGroupWarning(Topology const& topIn, int rnum) {
  Residue const& res = topIn.Res(rnum);
  if ( res.Name() == "SO3" ) {
    mprintf("Warning: Residue '%s'; after LEaP, will need to adjust the charge on the link oxygen by +0.031.\n",
            topIn.TruncResNameNum(rnum).c_str());
  } else if ( res.Name() == "MEX" ) {
    mprintf("Warning: Residue '%s'; after LEaP, will need to adjust the charge on the carbon bonded to link oxygen by -0.039.\n",
            topIn.TruncResNameNum(rnum).c_str());
  } else if ( res.Name() == "ACX" ) {
    mprintf("Warning: Residue '%s'; after LEaP, will need to adjust the charge on the carbon bonded to link oxygen by +0.008.\n",
            topIn.TruncResNameNum(rnum).c_str());
  }
}

// Exec_PrepareForLeap::Help()
void Exec_PrepareForLeap::Help() const
{
  mprintf("\tcrdset <coords set> [frame <#>] name <out coords set>\n"
          "\t[pdbout <pdbfile> [terbymol]]\n"
          "\t[leapunitname <unit>] [out <leap input file> [runleap <ff file>]]\n"
          "\t[skiperrors]\n"
          "\t[nowat [watermask <watermask>] [noh]\n"
          "\t[keepaltloc {<alt loc ID>|highestocc}]\n"
          "\t[stripmask <stripmask>] [solventresname <solventresname>]\n"
          "\t[molmask <molmask> ...] [determinemolmask <mask>]\n"
          "\t[{nohisdetect |\n"
          "\t  [nd1 <nd1>] [ne2 <ne2] [hisname <his>] [hiename <hie>]\n"
          "\t  [hidname <hid>] [hipname <hip]}]\n"
          "\t[{nodisulfides |\n"
          "\t  existingdisulfides |\n"
          "\t  [cysmask <cysmask>] [disulfidecut <cut>] [newcysname <name>]}]\n"
          "\t[{nosugars |\n"
          "\t  sugarmask <sugarmask> [noc1search] [nosplitres]\n"
          "\t  [resmapfile <file>]\n"
          "\t  [hasglycam] [determinesugarsby {geom|name}]\n"
          "\t }]\n"
          "  Prepare the structure in the given coords set for easier processing\n"
          "  with the LEaP program from AmberTools. Any existing/potential\n"
          "  disulfide bonds will be identified and the residue names changed\n"
          "  to <name> (CYX by default), and if specified any sugars\n"
          "  recognized in the <sugarmask> region will be identified and have\n"
          "  their names changed to Glycam names. Disulfides and sugars will\n"
          "  have any inter-residue bonds removed, and the appropriate LEaP\n"
          "  input to add the bonds back once the structure has been loaded\n"
          "  into LEaP will be written to <leap input file>.\n"
         );
}

// Exec_PrepareForLeap::Execute()
Exec::RetType Exec_PrepareForLeap::Execute(CpptrajState& State, ArgList& argIn)
{
  mprintf("\tPREPAREFORLEAP:\n");
  mprintf("# Citation: Roe, D.R.; Bergonzo, C.; \"PrepareForLeap: An Automated Tool for\n"
          "#           Fast PDB-to-Parameter Generation.\"\n"
          "#           J. Comp. Chem. (2022), V. 43, I. 13, pp 930-935.\n" );
  debug_ = State.Debug();
  errorsAreFatal_ = !argIn.hasKey("skiperrors");
  // Get input coords
  std::string crdset = argIn.GetStringKey("crdset");
  if (crdset.empty()) {
    mprinterr("Error: Must specify input COORDS set with 'crdset'\n");
    return CpptrajState::ERR;
  }
  DataSet* ds = State.DSL().FindSetOfGroup( crdset, DataSet::COORDINATES );
  if (ds == 0) {
    mprinterr("Error: No COORDS set found matching %s\n", crdset.c_str());
    return CpptrajState::ERR;
  }
  DataSet_Coords& coords = static_cast<DataSet_Coords&>( *((DataSet_Coords*)ds) );
  // Get frame from input coords
  int tgtframe = argIn.getKeyInt("frame", 1) - 1;
  mprintf("\tUsing frame %i from COORDS set %s\n", tgtframe+1, coords.legend());
  if (tgtframe < 0 || tgtframe >= (int)coords.Size()) {
    mprinterr("Error: Frame is out of range.\n");
    return CpptrajState::ERR;
  }
  Frame frameIn = coords.AllocateFrame();
  coords.GetFrame(tgtframe, frameIn);


  // Copy input topology, may be modified.
  Topology topIn = coords.Top();

  // Allocate output COORDS data set
  std::string outname = argIn.GetStringKey("name");
  if (outname.empty()) {
    mprinterr("Error: Must specify output COORDS set with 'name'\n");
    return CpptrajState::ERR;
  }
  DataSet_Coords_CRD* outCoords = (DataSet_Coords_CRD*)State.DSL().AddSet( DataSet::COORDS, outname );
  if (outCoords == 0) {
    mprinterr("Error: Could not allocate output COORDS set.\n");
    return CpptrajState::ERR;
  }
  mprintf("\tPrepared system will be saved to COORDS set '%s'\n", outCoords->legend());

  std::string leapffname = argIn.GetStringKey("runleap");
  if (!leapffname.empty()) {
#   ifdef _MSC_VER
    mprinterr("Error: Cannot use LEaP interface on windows.\n");
    return CpptrajState::ERR;
#   else
    mprintf("\tWill attempt to run leap with force fields specified in file '%s'\n",
            leapffname.c_str());
#   endif
  }

  std::string pdbout = argIn.GetStringKey("pdbout");
  if (!pdbout.empty())
    mprintf("\tPDB will be written to %s\n", pdbout.c_str());
  else {
    if (!leapffname.empty()) {
      mprinterr("Error: Must specify PDB file name with 'pdbout' if 'runleap' specified.\n");
      return CpptrajState::ERR;
    }
  }
  std::string pdb_ter_arg;
  if (!argIn.hasKey("terbymol")) {
    mprintf("\tUsing original TER cards where possible.\n");
    pdb_ter_arg.assign("pdbter");
  } else
    mprintf("\tGenerating TER cards based on molecular connectivity.\n");

  std::string leapfilename = argIn.GetStringKey("out");
  if (!leapfilename.empty())
    mprintf("\tWriting leap input to '%s'\n", leapfilename.c_str());
  else {
    if (!leapffname.empty()) {
      mprinterr("Error: Must specify leap input file name with 'out' if 'runleap' specified.\n");
      return CpptrajState::ERR;
    }
  }
  leapunitname_ = argIn.GetStringKey("leapunitname", "m");
  mprintf("\tUsing leap unit name: %s\n", leapunitname_.c_str());
  if (validDouble(leapunitname_))
    mprintf("Warning: LEaP unit name '%s' is a valid number; this may confuse some LEaP commands.\n",
             leapunitname_.c_str());
  solventResName_ = argIn.GetStringKey("solventresname", "HOH");
  mprintf("\tSolvent residue name: %s\n", solventResName_.c_str());
  // TODO functional group stuff should be in a file

  bool prepare_sugars = !argIn.hasKey("nosugars");
  if (!prepare_sugars)
    mprintf("\tNot attempting to prepare sugars.\n");
  else
    mprintf("\tWill attempt to prepare sugars.\n");

  // Load PDB residue names recognized by amber
  if (LoadPdbResNames( argIn.GetStringKey("resnamefile" ) )) {
    mprinterr("Error: PDB residue name file load failed.\n");
    return CpptrajState::ERR;
  }
  mprintf("\t%zu PDB residue names recognized by Amber FFs.\n", pdb_res_names_.size());
  // DEBUG
  if (debug_ > 0) {
    mprintf("\tPDB residue names recognized by Amber FFs:\n");
    for (SetType::const_iterator it = pdb_res_names_.begin(); it != pdb_res_names_.end(); ++it)
      mprintf("\t  %s\n", *(*it));
  }

  // Load PDB to glycam residue name map
  SugarBuilder sugarBuilder;
  if (prepare_sugars) {
    if (sugarBuilder.LoadGlycamPdbResMap( argIn.GetStringKey("resmapfile" ) )) {
      mprinterr("Error: PDB to glycam name map load failed.\n");
      return CpptrajState::ERR;
    }
  }

  Iarray pdbResToRemove;
  std::string removeArg = argIn.GetStringKey("remove");
  if (!removeArg.empty()) {
    if (removeArg == "unrecognized") {
      mprintf("\tRemoving unrecognized PDB residues.\n");
      pdbResToRemove = GetUnrecognizedPdbResidues( topIn, sugarBuilder );
    } else if (removeArg == "isolated") {
      mprintf("\tRemoving unrecognized and isolated PDB residues.\n");
      Iarray unrecognizedPdbRes = GetUnrecognizedPdbResidues( topIn, sugarBuilder );
      pdbResToRemove = GetIsolatedUnrecognizedResidues( topIn, unrecognizedPdbRes );
    } else {
      mprinterr("Error: Unrecognized keyword for 'remove': %s\n", removeArg.c_str());
      return CpptrajState::ERR;
    }
  }

  // Deal with any coordinate modifications
  bool remove_water     = argIn.hasKey("nowat");
  std::string waterMask = argIn.GetStringKey("watermask", ":" + solventResName_);
  bool remove_h         = argIn.hasKey("noh");
  std::string altLocArg = argIn.GetStringKey("keepaltloc");
  if (!altLocArg.empty()) {
    if (altLocArg != "highestocc" &&
        altLocArg.size() > 1)
    {
      mprinterr("Error: Invalid keyword for 'keepaltloc' '%s'; must be 'highestocc' or 1 character.\n",
                altLocArg.c_str());
      return CpptrajState::ERR;
    }
  }
  std::string stripMask = argIn.GetStringKey("stripmask");

  // If keeping highest alt loc, check that alt locs and occupancies are present.
  if (altLocArg == "highestocc") {
    if (topIn.AtomAltLoc().empty()) {
      mprintf("Warning: 'highestocc' specified but no atom alternate location info.\n");
      altLocArg.clear();
    } else if (topIn.Occupancy().empty()) {
      mprintf("Warning: 'highestocc' specified but no atom occupancy info.\n");
      altLocArg.clear();
    }
  }
  // Check if alternate atom location IDs are present 
  if (!topIn.AtomAltLoc().empty()) {
    // For LEaP, must have only 1 atom alternate location
    char firstAltLoc = ' ';
    for (std::vector<char>::const_iterator altLocId = topIn.AtomAltLoc().begin();
                                           altLocId != topIn.AtomAltLoc().end();
                                         ++altLocId)
    {
      if (firstAltLoc == ' ') {
        // Find first non-blank alternate location ID
        if (*altLocId != ' ')
          firstAltLoc = *altLocId;
      } else if (*altLocId != ' ' && *altLocId != firstAltLoc) {
        // Choose a default if necessary
        if (altLocArg.empty()) {
          altLocArg.assign(1, firstAltLoc);
          mprintf("Warning: '%s' has atoms with multiple alternate location IDs, which\n"
                  "Warning:  are not supported by LEaP. Keeping only '%s'.\n"
                  "Warning: To choose a specific location to keep use the 'keepaltloc <char>'\n"
                  "Warning:  keyword.\n", coords.legend(), altLocArg.c_str());
         }
         break;
      }
    }
  }

  if (remove_water)
    mprintf("\tRemoving solvent. Solvent mask= '%s'\n", waterMask.c_str());
  if (remove_h)
    mprintf("\tRemoving hydrogens.\n");
  if (!altLocArg.empty())
    mprintf("\tIf present, keeping only alternate atom locations denoted by '%s'\n", altLocArg.c_str());
  if (!stripMask.empty())
    mprintf("\tRemoving atoms in mask '%s'\n", stripMask.c_str());
  if (ModifyCoords(topIn, frameIn, remove_water, altLocArg, stripMask, waterMask, pdbResToRemove))
  {
    mprinterr("Error: Modification of '%s' failed.\n", coords.legend());
    return CpptrajState::ERR;
  }

  // Do histidine detection before H atoms are removed
  if (!argIn.hasKey("nohisdetect")) {
    std::string nd1name = argIn.GetStringKey("nd1", "ND1");
    std::string ne2name = argIn.GetStringKey("ne2", "NE2");
    std::string hisname = argIn.GetStringKey("hisname", "HIS");
    std::string hiename = argIn.GetStringKey("hiename", "HIE");
    std::string hidname = argIn.GetStringKey("hidname", "HID");
    std::string hipname = argIn.GetStringKey("hipname", "HIP");
    mprintf("\tHistidine protonation detection:\n");
    mprintf("\t\tND1 atom name                   : %s\n", nd1name.c_str());
    mprintf("\t\tNE2 atom name                   : %s\n", ne2name.c_str());
    mprintf("\t\tHistidine original residue name : %s\n", hisname.c_str());
    mprintf("\t\tEpsilon-protonated residue name : %s\n", hiename.c_str());
    mprintf("\t\tDelta-protonated residue name   : %s\n", hidname.c_str());
    mprintf("\t\tDoubly-protonated residue name  : %s\n", hipname.c_str());
    // Add epsilon, delta, and double-protonated names as recognized.
    pdb_res_names_.insert( hiename );
    pdb_res_names_.insert( hidname );
    pdb_res_names_.insert( hipname );
    if (DetermineHisProt( topIn,
                          nd1name, ne2name,
                          hisname, hiename, hidname, hipname)) {
      mprinterr("Error: HIS protonation detection failed.\n");
      return CpptrajState::ERR;
    }
  }

  // Remove hydrogens
  if (remove_h) {
    if (RemoveHydrogens(topIn, frameIn)) return CpptrajState::ERR;
  }

  // See if sugars already have glycam names
  hasGlycam_ = argIn.hasKey("hasglycam");
  if (hasGlycam_)
    mprintf("\tAssuming sugars already have glycam residue names.\n");

  // Get sugar mask or default sugar mask
  std::string sugarmaskstr = argIn.GetStringKey("sugarmask");
  if (!sugarmaskstr.empty()) {
    if (!prepare_sugars) {
      mprinterr("Error: Cannot specify 'nosugars' and 'sugarmask'\n");
      return CpptrajState::ERR;
    }
  } else if (hasGlycam_) {
    sugarmaskstr = GenGlycamResMaskString();
  } else if (prepare_sugars) {
    // No sugar mask specified; create one from names in pdb_to_glycam_ map.
    sugarmaskstr.assign(":");
    for (MapType::const_iterator mit = pdb_to_glycam_.begin(); mit != pdb_to_glycam_.end(); ++mit)
    {
      if (mit != pdb_to_glycam_.begin())
        sugarmaskstr.append(",");
      sugarmaskstr.append( mit->first.Truncated() );
    }
  }

  // Get how sugars should be determined (geometry/name)
  std::string determineSugarsBy = argIn.GetStringKey("determinesugarsby", "geometry");
  if (determineSugarsBy == "geometry") {
    useSugarName_ = false;
    mprintf("\tWill determine sugar anomer type/configuration by geometry.\n");
  } else if (determineSugarsBy == "name") {
    useSugarName_ = true;
    mprintf("\tWill determine sugar anomer type/configuration from residue name.\n");
  } else {
    mprinterr("Error: Invalid argument for 'determinesugarsby': %s\n", determineSugarsBy.c_str());
    return CpptrajState::ERR;
  }

  // If preparing sugars, need to set up an atom map and potentially
  // search for terminal sugars/missing bonds. Do this here after all atom
  // modifications have been done.
  std::vector<Sugar> sugarResidues;
  if (prepare_sugars) {
    // Set up an AtomMap for this residue to help determine stereocenters.
    // This is required by the IdSugarRing() function.
    //myMap_.SetDebug(debug_); // DEBUG
    if (myMap_.Setup(topIn, frameIn)) {
      mprinterr("Error: Atom map setup failed\n");
      return CpptrajState::ERR;
    }
    myMap_.DetermineAtomIDs();

    bool splitres = !argIn.hasKey("nosplitres");
    if (splitres)
      mprintf("\tWill split off recognized sugar functional groups into separate residues.\n");
    else
      mprintf("\tNot splitting recognized sugar functional groups into separate residues.\n");
    bool c1bondsearch = !argIn.hasKey("noc1search");
    if (c1bondsearch)
      mprintf("\tWill search for missing bonds to sugar anomeric atoms.\n");
    else
      mprintf("\tNot searching for missing bonds to sugar anomeric atoms.\n");
    // May need to modify sugar structure/topology, either by splitting
    // C1 hydroxyls of terminal sugars into ROH residues, and/or by
    // adding missing bonds to C1 atoms.
    // This is done before any identification takes place since we want
    // to identify based on the most up-to-date topology.
    if (FixSugarsStructure(sugarResidues, sugarmaskstr, topIn, frameIn,
                           c1bondsearch, splitres))
    {
      mprinterr("Error: Sugar structure modification failed.\n");
      return CpptrajState::ERR;
    }
    // NOTE: If IdSugarRing() is to be used after this point, the map
    //       will need to be recreated.
    // Since FixSugarsStructure() can re-order atoms, need
    // to recreate the map.
    //myMap_.ClearMap();
    //if (myMap_.Setup(topIn, frameIn)) {
    //  mprinterr("Error: Atom map second setup failed\n");
    //  return CpptrajState::ERR;
    //}
    //myMap_.DetermineAtomIDs();
  }

  // ----- Below here, no more removing/reordering atoms. ------------ 

  // Each residue starts out unknown.
  ResStatArray resStat( topIn.Nres() );

  // Get masks for molecules now since bond info in topology may be modified later.
  std::vector<AtomMask> molMasks;
  std::string mstr = argIn.GetStringKey("molmask");
  while (!mstr.empty()) {
    mprintf("\tAll atoms selected by '%s' will be in same molecule.\n", mstr.c_str());
    molMasks.push_back( AtomMask() );
    if (molMasks.back().SetMaskString( mstr )) {
      mprinterr("Error: Invalid mask.\n");
      return CpptrajState::ERR;
    }
    if (topIn.SetupIntegerMask( molMasks.back() )) return CpptrajState::ERR;
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
    if (topIn.SetupCharMask( determineMolMask )) return CpptrajState::ERR;
    determineMolMask.MaskInfo();
    if (determineMolMask.None()) {
      mprinterr("Error: Nothing selected by mask.\n");
      return CpptrajState::ERR;
    }
  }

  //CpptrajFile* outfile = State.DFL().AddCpptrajFile(leapfilename,
  //                                                  "LEaP Input", DataFileList::TEXT, true);
  // NOTE: This needs to contain ONLY leap input, so dont put it on the master file list
  CpptrajFile LEAPOUT;
  if (LEAPOUT.OpenWrite(leapfilename)) return CpptrajState::ERR;
  CpptrajFile* outfile = &LEAPOUT;
  if (outfile == 0) return CpptrajState::ERR;
  mprintf("\tLEaP input containing 'loadpdb' and bond commands for disulfides,\n"
          "\t  sugars, etc will be written to '%s'\n", outfile->Filename().full());
  // Add the loadpdb command if we are writing a PDB file.
  // TODO add 'addPdbResMap { { 1 "NH2" "NHE" } }' to recognize NHE?
  if (!pdbout.empty())
    outfile->Printf("%s = loadpdb %s\n", leapunitname_.c_str(), pdbout.c_str());

  // Disulfide search
  if (!argIn.hasKey("nodisulfides")) {
    if (SearchForDisulfides( argIn.getKeyDouble("disulfidecut", 2.5),
                             argIn.GetStringKey("newcysname", "CYX"),
                             argIn.GetStringKey("cysmask", ":CYS@SG"),
                            !argIn.hasKey("existingdisulfides"),
                             topIn, frameIn, outfile ))
    {
      mprinterr("Error: Disulfide search failed.\n");
      return CpptrajState::ERR;
    }
  } else {
    mprintf("\tNot searching for disulfides.\n");
  }

  // Prepare sugars
  if (prepare_sugars) {
    if (PrepareSugars(sugarmaskstr, sugarResidues, topIn, frameIn, outfile)) {
      mprinterr("Error: Sugar preparation failed.\n");
      return CpptrajState::ERR;
    }
  } else {
    mprintf("\tNot preparing sugars.\n");
  }

  // Count any solvent molecules
  if (!remove_water) {
    NameType solvName(solventResName_);
    unsigned int nsolvent = 0;
    for (Topology::res_iterator res = topIn.ResStart(); res != topIn.ResEnd(); ++res) {
      if ( res->Name() == solvName) {
        nsolvent++;
        resStat[res-topIn.ResStart()] = ResStatArray::VALIDATED;
        // Set as terminal; TODO is this needed? Leap seems ok with not having TER for HOH
        topIn.SetRes(res-topIn.ResStart()).SetTerminal(true);
      }
    }
    if (nsolvent > 0) mprintf("\t%u solvent residues.\n", nsolvent);
  }

  // Residue validation.
  //mprintf("\tResidues with potential problems:\n");
  int fatal_errors = 0;
  static const char* msg1 = "Potential problem : ";
  static const char* msg2 = "Fatal problem     : ";
  for (ResStatArray::iterator it = resStat.begin(); it != resStat.end(); ++it)
  {
    LeapFxnGroupWarning(topIn, it-resStat.begin());
    //if ( *it == VALIDATED )
    //  mprintf("\t\t%s VALIDATED\n", topIn.TruncResNameOnumId(it-resStat_.begin()).c_str());
    //else
    //  mprintf("\t\t%s UNKNOWN\n", topIn.TruncResNameOnumId(it-resStat_.begin()).c_str());
    // ----- Warnings --------
    if ( *it == ResStatArray::UNKNOWN ) {
      SetType::const_iterator pname = pdb_res_names_.find( topIn.Res(it-resStat.begin()).Name() );
      if (pname == pdb_res_names_.end())
        mprintf("\t%s%s is an unrecognized name and may not have parameters.\n",
                msg1, topIn.TruncResNameOnumId(it-resStat.begin()).c_str());
      else
        *it = VALIDATED;
    } else if ( *it == ResStatArray::SUGAR_NAME_MISMATCH ) {
        mprintf("\t%s%s sugar anomer type and/or configuration is not consistent with name.\n",
                msg1, topIn.TruncResNameOnumId(it-resStat.begin()).c_str());
    // ----- Fatal Errors ----
    } else if ( *it == ResStatArray::SUGAR_UNRECOGNIZED_LINK_RES ) {
        mprintf("\t%s%s is linked to a sugar but has no sugar-linkage form.\n",
                msg2, topIn.TruncResNameOnumId(it-resStat.begin()).c_str());
        fatal_errors++;
    } else if ( *it == ResStatArray::SUGAR_UNRECOGNIZED_LINKAGE ) {
        mprintf("\t%s%s is a sugar with an unrecognized linkage.\n",
                msg2, topIn.TruncResNameOnumId(it-resStat.begin()).c_str());
        fatal_errors++;
    } else if ( *it == ResStatArray::SUGAR_NO_LINKAGE ) {
        mprintf("\t%s%s is an incomplete sugar with no linkages.\n",
                msg2, topIn.TruncResNameOnumId(it-resStat.begin()).c_str());
        fatal_errors++;
    } else if ( *it == ResStatArray::SUGAR_NO_CHAIN_FOR_LINK ) {
        mprintf("\t%s%s could not identify chain atoms for determining linkages.\n",
                msg2, topIn.TruncResNameOnumId(it-resStat.begin()).c_str());
        fatal_errors++;
/*    } else if ( *it == SUGAR_MISSING_C1X ) { // TODO should this be a warning
        mprintf("\t%s%s Sugar is missing anomeric carbon substituent.\n",
                msg2, topIn.TruncResNameOnumId(it-resStat_.begin()).c_str());
        fatal_errors++;*/
    } else if ( *it == ResStatArray::SUGAR_SETUP_FAILED ) {
        mprintf("\t%s%s Sugar setup failed and could not be identified.\n",
                msg2, topIn.TruncResNameOnumId(it-resStat.begin()).c_str());
        fatal_errors++;
    }
  }

  // Try to set terminal residues
  if (!molMasks.empty() || determineMolMask.MaskStringSet()) {
    // Reset terminal status
    for (int rnum = 0; rnum != topIn.Nres(); rnum++)
      topIn.SetRes(rnum).SetTerminal(false);
    // The final residue of each molMask is terminal
    for (std::vector<AtomMask>::const_iterator mask = molMasks.begin();
                                               mask != molMasks.end(); ++mask)
    {
      //std::vector<int> Rnums = coords.Top().ResnumsSelectedBy( *mask );
      int lastAtom = mask->back();
      int lastRes = topIn[lastAtom].ResNum();
      mprintf("\tSetting residue %s as terminal.\n",
        topIn.TruncResNameOnumId(lastRes).c_str());
      topIn.SetRes(lastRes).SetTerminal( true );
    }
    // Set ter based on connectivity
    if (determineMolMask.MaskStringSet()) {
      if (FindTerByBonds(topIn, determineMolMask)) {
        mprinterr("Error: Could not set TER by connectivity.\n");
        return CpptrajState::ERR;
      }
    }
  }

  // Setup output COORDS
  outCoords->CoordsSetup( topIn, coords.CoordsInfo() );
  outCoords->AddFrame( frameIn );

  if (!pdbout.empty()) {
    Trajout_Single PDB;
    PDB.SetDebug( debug_ );
    if (PDB.InitTrajWrite( pdbout, "topresnum " + pdb_ter_arg, State.DSL(), TrajectoryFile::PDBFILE)) {
      mprinterr("Error: Could not initialize output PDB\n");
      return CpptrajState::ERR;
    }
    if (PDB.SetupTrajWrite(outCoords->TopPtr(), outCoords->CoordsInfo(), 1)) {
      mprinterr("Error: Could not set up output PDB\n");
      return CpptrajState::ERR;
    }
    PDB.PrintInfo(1);
    PDB.WriteSingle(0, frameIn);
    PDB.EndTraj();
  }

  outfile->CloseFile();

  if (fatal_errors > 0) {
    if (errorsAreFatal_) {
      mprinterr("Error: %i errors were encountered that will prevent LEaP from running successfully.\n", fatal_errors);
      return CpptrajState::ERR;
    } else {
      mprintf("Warning: %i errors were encountered that will prevent LEaP from running successfully.\n", fatal_errors);
      mprintf("Warning: Continuing on anyway, but final structure **NEEDS VALIDATION**.\n");
    }
  }
  // Run leap if needed
  if (!leapffname.empty()) {
    if (RunLeap( leapffname, leapfilename )) {
      mprinterr("Error: Running leap failed.\n");
      return CpptrajState::ERR;
    }
  }

  return CpptrajState::OK;
}
