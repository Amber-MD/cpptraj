#include "Exec_PrepareForLeap.h"
#include "CpptrajStdio.h"
#include "DistRoutines.h"
#include "CharMask.h"
#include "TorsionRoutines.h"
#include "Constants.h"
#include "CpptrajFile.h"
#include "Trajout_Single.h"
#include <stack>
#include <cctype> // tolower
#include <algorithm> // sort

// Exec_PrepareForLeap::Help()
void Exec_PrepareForLeap::Help() const
{
  mprintf("\tcrdset <coords set> [frame <#>] [out <file>]\n"
          "\t[{nodisulfides |\n"
          "\t  existingdisulfides |\n"
          "\t  [cysmask <cysmask>] [disulfidecut <cut>] [newcysname <name>]}]\n"
          "\t[{nosugars | sugarmask <sugarmask>}] [resmapfile <file>]\n"
          "\t[leapunitname <unit>] [pdbout <pdbout>]\n"
          "\t[molmask <molmask> ...] [determinemolmask <mask>]\n"
          "  Prepare the structure in the given coords set for easier processing\n"
          "  with the LEaP program from AmberTools. Any existing/potential\n"
          "  disulfide bonds will be identified and the residue names changed\n"
          "  to <name> (CYX by default), and if specified any sugars\n"
          "  recognized in the <sugarmask> region will be identified and have\n"
          "  their names changed to Glycam names. Disulfides and sugars will\n"
          "  have any inter-residue bonds removed, and the appropriate LEaP\n"
          "  input to add the bonds back once the structure has been loaded\n"
          "  into LEaP will be written to <file>.\n"
         );
}

/// Used to change residue name to nameIn
static inline void ChangeResName(Residue& res, NameType const& nameIn) {
  if (res.Name() != nameIn) {
    mprintf("\t    Changing residue %s to %s\n", *(res.Name()), *nameIn);
    res.SetName( nameIn );
  }
}

/// Used to change atom name to nameIn
static inline void ChangeAtomName(Atom& atm, NameType const& nameIn) {
  if (atm.Name() != nameIn) {
    mprintf("\t    Changing atom %s to %s\n", *(atm.Name()), *nameIn);
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
    case 'S':
      if      (linkstr == "C2") linkcode = "0";
      break;
    default:
      if      (linkstr == "C1") linkcode = "0";
      else if (linkstr == "C1O2") linkcode = "2";
      else if (linkstr == "C1O3") linkcode = "3";
      else if (linkstr == "C1O4") linkcode = "4";
      else if (linkstr == "C1O6") linkcode = "6";
      else if (linkstr == "C1O2O3") linkcode = "Z";
      else if (linkstr == "C1O2O4") linkcode = "Y";
      else if (linkstr == "C1O2O6") linkcode = "X";
      else if (linkstr == "C1O3O4") linkcode = "W";
      else if (linkstr == "C1O3O6") linkcode = "V";
      else if (linkstr == "C1O4O6") linkcode = "U";
      else if (linkstr == "C1O2O3O4") linkcode = "T";
      else if (linkstr == "C1O2O3O6") linkcode = "S";
      else if (linkstr == "C1O2O4O6") linkcode = "R";
      else if (linkstr == "C1O3O4O6") linkcode = "Q";
      else if (linkstr == "C1O2O3O4O6") linkcode = "P";
      break;
  }
  if (linkcode.empty())
    mprintf("Warning: Could not determine link code.\n");
  return linkcode;
}

/** Load reduced interal PDB to Glycam map. */
void Exec_PrepareForLeap::SetGlycamPdbResMap() {
  pdb_to_glycam_.insert( PairType("NAG", 'Y') );
  pdb_to_glycam_.insert( PairType("FUC", 'F') );
  pdb_to_glycam_.insert( PairType("GAL", 'L') );
  pdb_to_glycam_.insert( PairType("BMA", 'M') );
  pdb_to_glycam_.insert( PairType("MAN", 'M') );
}

/** Load PDB to Glycam residue map from file. */
int Exec_PrepareForLeap::LoadGlycamPdbResMap(std::string const& fnameIn)
{
  std::string fname;
  if (fnameIn.empty()) {
    // Check CPPTRAJHOME
    const char* env = getenv("CPPTRAJHOME");
    if (env != 0) {
      fname.assign(env);
      fname += "/dat/Carbohydrate_PDB_Glycam_Names.txt";
    }
    mprintf("Info: Parameter file path from CPPTRAJHOME variable: '%s'\n", fname.c_str());
  }
  if (fname.empty()) {
    mprintf("Warning: No PDB->Glycam file specified and/or CPPTRAJHOME not set.\n"
            "Warning: Using only basic PDB residue name recognition.\n");
    SetGlycamPdbResMap();
    return 0;
  }
  mprintf("\tReading PDB residue name -> Glycam name map from '%s'\n", fname.c_str());

  CpptrajFile infile;
  if (infile.OpenRead(fname)) {
    mprinterr("Error: Could not open Glycam residue map file.\n");
    return 1;
  }
  const char* ptr = 0;
  while ( (ptr = infile.NextLine()) != 0 ) {
    ArgList argline( ptr, " " );
    if (argline[0][0] != '#') {
      if (argline.Nargs() != 3) {
        mprinterr("Error: Expected only 3 data columns in '%s', got %i\n",
                  infile.Filename().full(), argline.Nargs());
        mprinterr("Error: %s\n", ptr);
        return 1;
      }
      ArgList pdbnames( argline[2], "," );
      if (pdbnames.Nargs() < 1) {
        mprinterr("Error: No pdb names found.\n");
        mprinterr("Error: %s\n", ptr);
        return 1;
      }
      // TODO handle glycam res names with > 1 char
      for (int n = 0; n < pdbnames.Nargs(); n++)
        pdb_to_glycam_.insert( PairType(pdbnames[n], argline[1][0]) );
    }
  } // END loop over file
  infile.CloseFile();

  return 0;
}

/// Recursive function for following bonds of an atom to a target atom
static void FollowBonds(int atm, Topology const& topIn, int idx, std::vector<int>& ring_atoms, int tgt_atom, std::vector<bool>& Visited, bool& found)
{
  Visited[atm] = true;
  for (int i = 0; i != idx; i++)
    mprintf("\t");
  mprintf("At atom %s\n", topIn.ResNameNumAtomNameNum(atm).c_str());
  ring_atoms[idx] = atm;
  // Assume we have started at the target atom
  if (idx > 0 && atm == tgt_atom) {
    found = true;
    return;
  }
  // Follow all atoms bonded to this atom
  for (Atom::bond_iterator bat = topIn[atm].bondbegin(); bat != topIn[atm].bondend(); ++bat)
  {
    if (!Visited[*bat]) {
      FollowBonds( *bat, topIn, idx+1, ring_atoms, tgt_atom, Visited, found );
      if (found) return;
    }
  }
}

/** Attempt to identify sugar residue, form, and linkages. */
int Exec_PrepareForLeap::IdentifySugar(int rnum, Topology* topIn,
                                       Frame const& frameIn, CharMask const& cmask,
                                       CpptrajFile* outfile, std::set<BondType>& sugarBondsToRemove)
const
{
  Residue& res = topIn->SetRes(rnum);
  // Try to ID the base sugar type from the input name.
  char resChar = ' ';

  MapType::const_iterator pdb_glycam = pdb_to_glycam_.find( res.Name() );
  if ( pdb_glycam == pdb_to_glycam_.end() ) { 
    mprinterr("Error: Could not identify sugar from residue name '%s'\n", *res.Name());
    return 1;
  }
  resChar = pdb_glycam->second;

  mprintf("\tSugar %s %i glycam name: %c\n", *res.Name(), rnum+1, resChar);

  // Change PDB names to Glycam ones
  for (int at = res.FirstAtom(); at != res.LastAtom(); at++)
  {
    if (resChar == 'Y') {
      if ( (*topIn)[at].Name() == "C7" )
        ChangeAtomName(topIn->SetAtom(at), "C2N");
      else if ( (*topIn)[at].Name() == "O7" )
        ChangeAtomName(topIn->SetAtom(at), "O2N");
      else if ( (*topIn)[at].Name() == "C8" )
        ChangeAtomName(topIn->SetAtom(at), "CME");
    } else if (resChar == 'S') {
      if ( (*topIn)[at].Name() == "C10" )
        ChangeAtomName(topIn->SetAtom(at), "C5N");
      else if ( (*topIn)[at].Name() == "O10" )
        ChangeAtomName(topIn->SetAtom(at), "O5N");
      else if ( (*topIn)[at].Name() == "C11" )
        ChangeAtomName(topIn->SetAtom(at), "CME");
    }
  }

  // Try to identify the sugar ring. Potential starting atoms are oxygens
  // bonded to two carbon atoms.
  std::vector<int> potentialRingStartAtoms;
  for (int at = res.FirstAtom(); at != res.LastAtom(); at++)
  {
    Atom const& currentAtom = (*topIn)[at];
    if (currentAtom.Element() == Atom::OXYGEN)
    {
      if (currentAtom.Nbonds() == 2) {
        if ( (*topIn)[currentAtom.Bond(0)].Element() == Atom::CARBON &&
             (*topIn)[currentAtom.Bond(0)].ResNum() == rnum &&
             (*topIn)[currentAtom.Bond(1)].Element() == Atom::CARBON &&
             (*topIn)[currentAtom.Bond(1)].ResNum() == rnum )
        {
          potentialRingStartAtoms.push_back( at );
        }
      }
    }
  }
  // TODO handle case where multiple potential ring start atoms exist
  if (potentialRingStartAtoms.empty()) {
    mprinterr("Error: Ring oxygen could not be identified\n");
    return 1;
  } else if (potentialRingStartAtoms.size() > 1) {
    mprinterr("Error: Multiple potential ring start atoms:\n");
    for (std::vector<int>::const_iterator it = potentialRingStartAtoms.begin();
                                          it != potentialRingStartAtoms.end();
                                         ++it)
      mprinterr("Error:   %s\n", topIn->ResNameNumAtomNameNum(*it).c_str());
    return 1;
  }

  // Out of the potential ring start atoms, see which ones are actually
  // part of a ring. Potential ring start atoms only have 2 bonds,
  // each one to a carbon.
  int n_ring_atoms = 0;
  int ring_oxygen_atom = -1; // e.g. O5
  std::vector<bool> Visited;
  Visited.reserve( topIn->Natom() );
  for (std::vector<int>::const_iterator ringat = potentialRingStartAtoms.begin();
                                        ringat != potentialRingStartAtoms.end();
                                      ++ringat)
  {
    Visited.assign( topIn->Natom(), true );
    for (int at = res.FirstAtom(); at != res.LastAtom(); at++)
      if (at != *ringat)
        Visited[at] = false;
    std::vector<int> ring_atoms( topIn->Res(rnum).NumAtoms(), -1 );
    bool ring_complete = false;
    // Since we have already established that *ringat is an oxygen bonded
    // to two carbons, just start at the first carbon to see if we can
    // get to the second carbon.
    FollowBonds((*topIn)[*ringat].Bond(0), *topIn, 0, ring_atoms,
                (*topIn)[*ringat].Bond(1), Visited, ring_complete);
    mprintf("DEBUG: Potential ring start atom %s, Ring complete = %i",
            topIn->ResNameNumAtomNameNum(*ringat).c_str(), (int)ring_complete);
    // TODO handle the case where multiple potential ring start atoms exist
    if (ring_complete) {
      ring_oxygen_atom = *ringat;
      // Use Visited as a mask with ring atoms
      Visited.assign( topIn->Natom(), false );
      Visited[ring_oxygen_atom] = true;
      n_ring_atoms = 1;
      mprintf(" :"); // DEBUG
      for (std::vector<int>::const_iterator it = ring_atoms.begin(); it != ring_atoms.end(); ++it)
      {
        mprintf(" %i", *it + 1);
        if (*it == -1) break;
        Visited[*it] = true;
        ++n_ring_atoms;
      }
      mprintf("\n"); // DEBUG
    }
  }
  mprintf("\t  Number of ring atoms= %i\n", n_ring_atoms);
  if (n_ring_atoms == 0 || ring_oxygen_atom == -1) {
    mprinterr("Error: Sugar ring atoms could not be identified.\n");
    return 1;
  }



  // Try to identify the ring oxygen. Should be bonded to 2 carbons
  // in the same residue.
  int ring_c_beg       = -1; // e.g. C1, "lowest" ring C, anomeric carbon
  int ring_c_beg_C     = -1; // e.g. C2, next C in ring bonded to anomeric carbon
  int ring_c_beg_X     = -1; // from anomeric C to non ring atom, non hydrogen 
  int ring_c_end       = -1; // e.g. C5, "highest" ring C, a chiral center
  int ring_c_end_C     = -1; // e.g. C4, from c_end to previous C in ring
  int ring_c_end_X     = -1; // e.g. C6, from c_end to non ring atom in same residue

  // TODO this is temporary
  if (ring_oxygen_atom != -1) {
    Atom const& currentAtom = (*topIn)[ring_oxygen_atom];
    if (currentAtom.Bond(0) < currentAtom.Bond(1)) {
      ring_c_beg = currentAtom.Bond(0);
      ring_c_end = currentAtom.Bond(1);
    } else {
      ring_c_beg = currentAtom.Bond(1);
      ring_c_end = currentAtom.Bond(0);
    }
  }

/*  for (int at = res.FirstAtom(); at != res.LastAtom(); at++)
  {
    Atom const& currentAtom = (*topIn)[at];
    if (currentAtom.Element() == Atom::OXYGEN)
    {
      if (currentAtom.Nbonds() == 2) {
        if ( (*topIn)[currentAtom.Bond(0)].Element() == Atom::CARBON &&
             (*topIn)[currentAtom.Bond(0)].ResNum() == rnum &&
             (*topIn)[currentAtom.Bond(1)].Element() == Atom::CARBON &&
             (*topIn)[currentAtom.Bond(1)].ResNum() == rnum )
        {
          if (ring_oxygen_atom != -1) {
            mprinterr("Error: Two potential ring oxygen atoms: %i and %i\n",
                      at+1, ring_oxygen_atom+1);
            return 1;
          }
          ring_oxygen_atom = at;
          if (currentAtom.Bond(0) < currentAtom.Bond(1)) {
            ring_c_beg = currentAtom.Bond(0);
            ring_c_end = currentAtom.Bond(1);
          } else {
            ring_c_beg = currentAtom.Bond(1);
            ring_c_end = currentAtom.Bond(0);
          }
        }
      }
    }
  }*/
  if (ring_oxygen_atom == -1) {
    mprinterr("Error: Ring oxygen atom could not be identified.\n");
    return 1;
  }
  if (ring_c_beg == -1) {
    mprinterr("Error: Lowest ring C could not be identified.\n");
    return 1;
  }
  if (ring_c_end == -1) {
    mprinterr("Error: Highest ring C could not be identified.\n");
    return 1;
  }
  mprintf("\t  Ring Cbeg-O-Cend atoms: %s-%s-%s\n",
          topIn->ResNameNumAtomNameNum(ring_c_beg).c_str(),
          topIn->ResNameNumAtomNameNum(ring_oxygen_atom).c_str(),
          topIn->ResNameNumAtomNameNum(ring_c_end).c_str());

  // Try to identify the ring atoms. Start from ring_c_beg, get to ring_c_end
/*  std::vector<bool> Visited( topIn->Natom(), true );
  for (int at = res.FirstAtom(); at != res.LastAtom(); at++)
    if (at != ring_oxygen_atom)
      Visited[at] = false;
  std::vector<int> ring_atoms( topIn->Res(rnum).NumAtoms(), -1 );
  bool ring_complete = false;
  FollowBonds( ring_c_beg, *topIn, 0, ring_atoms, ring_c_end, Visited, ring_complete );
  mprintf("DEBUG: Ring %i:", (int)ring_complete);
  // Use Visited as a mask with ring atoms
  Visited.assign( topIn->Natom(), false );
  Visited[ring_oxygen_atom] = true;
  int n_ring_atoms = 1;
  for (std::vector<int>::const_iterator it = ring_atoms.begin(); it != ring_atoms.end(); ++it)
  {
    mprintf(" %i", *it + 1);
    if (*it == -1) break;
    Visited[*it] = true;
    ++n_ring_atoms;
  }
  mprintf("\n");
  mprintf("\t  Number of ring atoms= %i\n", n_ring_atoms);
  if (!ring_complete) {
    mprinterr("Error: Sugar ring atoms could not be identified.\n");
    return 1;
  }*/

  // Get the substituent of first ring C (e.g. C1) that is to a non-ring atom, non hydrogen 
  // Get the substituent of first ring C (e.g. C1) that is part of the ring (e.g. C2)
  for ( Atom::bond_iterator bat = (*topIn)[ring_c_beg].bondbegin();
                            bat != (*topIn)[ring_c_beg].bondend();
                          ++bat )
  {
    if ( (*topIn)[*bat].Element() != Atom::HYDROGEN &&
         !Visited[*bat] )
    {
      if (ring_c_beg_X != -1) {
        // If there are two non-ring, non-hydrogen substituents, prioritize
        // the one that is part of this residue.
        bool bat_in_res = ((*topIn)[*bat].ResNum() == rnum);
        bool X_in_res   = ((*topIn)[ring_c_beg_X].ResNum() == rnum);
        if ( (bat_in_res && X_in_res) || (!bat_in_res && !X_in_res) ) {
          mprinterr("Error: Two potential substituents for anomeric carbon: %s and %s\n",
                    topIn->ResNameNumAtomNameNum(*bat).c_str(),
                    topIn->ResNameNumAtomNameNum(ring_c_beg_X).c_str());
          return 1;
        } else if (bat_in_res) {
          ring_c_beg_X = *bat;
        }
      } else
        ring_c_beg_X = *bat;
    } else if ( (*topIn)[*bat].Element() == Atom::CARBON &&
                Visited[*bat] )
    {
      if (ring_c_beg_C != -1) {
        mprinterr("Error: Two potential ring carbons bonded to anomeric carbon: %s and %s\n",
                  topIn->ResNameNumAtomNameNum(*bat).c_str(),
                  topIn->ResNameNumAtomNameNum(ring_c_beg_C).c_str());
        return 1;
      }
      ring_c_beg_C = *bat;
    }
  }
  if (ring_c_beg_X == -1) {
    // If the Cx (C1 substituent, usually a different residue) index is
    // not found this usually means missing inter-residue bond.
    // Alternatively, this could be an isolated sugar missing an -OH
    // group, so make this non-fatal.
    mprintf("Warning: Ring C1 substituent could not be identified.\n"
            "Warning:   If '%s' is from a topology without complete bonding information\n"
            "Warning:   (e.g. a PDB file), try loading the topology with the\n"
            "Warning:   'searchtype grid' keywords instead.\n", topIn->c_str());
    mprintf("Warning: This can also happen for isolated sugars missing e.g. a -OH\n"
            "Warning:   group. In that case coordinates for the missing sugar atoms\n"
            "Warning:   may need to be generated.\n");
    return 0;
  }
  mprintf("\t  C1 X substituent: %s\n",
          topIn->ResNameNumAtomNameNum(ring_c_beg_X).c_str());
  if (ring_c_beg_C == -1) {
    mprinterr("Error: Next ring atom after C1 could not be identified.\n");
    return 1;
  }
  mprintf("\t  C1 C ring substituent: %s\n",
          topIn->ResNameNumAtomNameNum(ring_c_beg_C).c_str());

  // Get substituent of last ring C (e.g. C5) that is not part of the ring (e.g. C6)
  // Get substituent of last ring C (e.g. C5) that is part of the ring (e.g. C4)
  for ( Atom::bond_iterator bat = (*topIn)[ring_c_end].bondbegin();
                            bat != (*topIn)[ring_c_end].bondend();
                          ++bat )
  {
    if ( (*topIn)[*bat].Element() == Atom::CARBON &&
         !Visited[*bat] )
    {
      if (ring_c_end_X != -1) {
        mprinterr("Error: Two potential substituents for ring end C: %s and %s\n",
                  topIn->ResNameNumAtomNameNum(*bat).c_str(),
                  topIn->ResNameNumAtomNameNum(ring_c_end_X).c_str());
        return 1;
      }
      ring_c_end_X = *bat;
    } else if ( (*topIn)[*bat].Element() == Atom::CARBON &&
      Visited[*bat] )
    {
      if (ring_c_end_C != -1) {
        mprinterr("Error: Two potential ring atoms bonded to ring end C: %s and %s\n",
                  topIn->ResNameNumAtomNameNum(*bat).c_str(),
                  topIn->ResNameNumAtomNameNum(ring_c_end_C).c_str());
        return 1;
      }
      ring_c_end_C = *bat;
    }
  }
  if (ring_c_end_X == -1) {
    mprinterr("Error: Highest ring C substituent could not be identified.\n");
    return 1;
  }
  mprintf("\t  C5 non-ring substituent: %s\n",
          topIn->ResNameNumAtomNameNum(ring_c_end_X).c_str());
  if (ring_c_end_C == -1) {
    mprinterr("Error: Highest ring C previous ring atom could not be identified.\n");
    return 1;
  }
  mprintf("\t  C5 previous ring atom: %s\n",
          topIn->ResNameNumAtomNameNum(ring_c_end_C).c_str());

  // Try to identify alpha/beta.
  // The alpha form has the CH2OH substituent (C5-C6 etc in Glycam) on the 
  // opposite side of the OH on the anomeric carbon (C1 in Glycam), while
  // in the beta form it is on the same side, i.e.
  // Alpha - C1 and C5 substituents are on opposite sides.
  // Beta  - C1 and C5 substituents are on the same side.
  std::string formStr;
  // Determine alpha/beta
  double t_c5 = Torsion( frameIn.XYZ(ring_c_end_C),     frameIn.XYZ(ring_c_end),
                         frameIn.XYZ(ring_oxygen_atom), frameIn.XYZ(ring_c_end_X) );
  double t_c1 = Torsion( frameIn.XYZ(ring_oxygen_atom), frameIn.XYZ(ring_c_beg),
                         frameIn.XYZ(ring_c_beg_C),     frameIn.XYZ(ring_c_beg_X) );
  mprintf("\t  A/B torsion around C5 = %f deg\n", t_c5 * Constants::RADDEG);
  mprintf("\t  A/B torsion around C1 = %f deg\n", t_c1 * Constants::RADDEG);
  bool c5up = (t_c5 > 0);
  bool c1up = (t_c1 > 0);
  if (c1up == c5up) {
    mprintf("\t  Beta form\n");
    formStr = "B";
  } else {
    mprintf("\t  Alpha form\n");
    formStr = "A";
  }
  // Determine D/L
  // Check the chirality around the C5 atom.
  bool isDform = !c5up;
  /*double torsion = Torsion( frameIn.XYZ(ring_c_end_C), frameIn.XYZ(ring_c_end),
                            frameIn.XYZ(ring_c_end_X), frameIn.XYZ(ring_oxygen_atom) );
  mprintf("\t  D/L Torsion around C5= %f deg\n", torsion * Constants::RADDEG);*/
  if (isDform) {
    mprintf("\t  D form\n");
  } else {
    mprintf("\t  L form\n");
    //isDform = false;
  }

  // Identify linkages to other residues.
  // Use a set to store linkages so they are in alphabetical order for easier identification.
  std::set<NameType> linkages;
  // Bonds to non sugars to be removed since these will confuse tleap
  BondArray bondsToRemove;
  // Loop over sugar atoms
  for (int at = res.FirstAtom(); at != res.LastAtom(); at++)
  {
    // Check for bonds to other residues
    for (Atom::bond_iterator bat = (*topIn)[at].bondbegin();
                             bat != (*topIn)[at].bondend(); ++bat)
    {
      if ((*topIn)[*bat].ResNum() != rnum) {
        if (!cmask.AtomInCharMask(*bat)) {
          mprintf("\t  Sugar %s bonded to non-sugar %s\n",
                  topIn->ResNameNumAtomNameNum(at).c_str(),
                  topIn->ResNameNumAtomNameNum(*bat).c_str());
          linkages.insert( (*topIn)[at].Name() );
          bondsToRemove.push_back( BondType(at, *bat, -1) );
          // Check if this is a recognized linkage to non-sugar TODO put in another file?
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
          mprintf("\t  Sugar %s bonded to sugar %s\n",
                  topIn->ResNameNumAtomNameNum(at).c_str(),
                  topIn->ResNameNumAtomNameNum(*bat).c_str());
          linkages.insert( (*topIn)[at].Name() );
          // Also remove inter-sugar bonds since leap cant handle branching
          if (at < *bat)
            sugarBondsToRemove.insert( BondType(at, *bat, -1) );
          else
            sugarBondsToRemove.insert( BondType(*bat, at, -1) );
        }
      }
    } // END loop over bonded atoms
  } // END loop over residue atoms

  // Determine linkage
  mprintf("\t  Link atoms:");
  for (std::set<NameType>::const_iterator it = linkages.begin();
                                          it != linkages.end(); ++it)
    mprintf(" %4s", *(*it));
  mprintf("\n");
  std::string linkcode = LinkageCode(resChar, linkages);
  mprintf("\t  Linkage code: %s\n", linkcode.c_str());
  if (linkcode.empty()) {
    mprinterr("Error: Unrecognized sugar linkage.\n");
    return 1;
  }
  // Modify residue char to indicate D form if necessary.
  // We do this here and not above so as not to mess with the
  // linkage determination.
  if (!isDform)
    resChar = tolower( resChar );
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

/** Prepare sugars for leap. */
int Exec_PrepareForLeap::PrepareSugars(AtomMask& sugarMask, DataSet_Coords& coords,
                                       Frame const& frameIn, CpptrajFile* outfile)
const
{
  std::set<BondType> sugarBondsToRemove;
  mprintf("\tPreparing sugars selected by '%s'\n", sugarMask.MaskString());
  if (coords.Top().SetupIntegerMask( sugarMask )) return 1;
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
      if (IdentifySugar(*rnum, coords.TopPtr(), frameIn, cmask, outfile, sugarBondsToRemove))
        return 1;
    } // END loop over sugar residues
    // Remove bonds between sugars
    for (std::set<BondType>::const_iterator bnd = sugarBondsToRemove.begin();
                                            bnd != sugarBondsToRemove.end(); ++bnd)
    {
      LeapBond(bnd->A1(), bnd->A2(), *(coords.TopPtr()), outfile);
      coords.TopPtr()->RemoveBond(bnd->A1(), bnd->A2());
    }
  }
  // Bonds to sugars have been removed, so regenerate molecule info
  coords.TopPtr()->DetermineMolecules();
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

/** Search for disulfide bonds. */
int Exec_PrepareForLeap::SearchForDisulfides(double disulfidecut, std::string const& newcysnamestr,
                                             std::string const& cysmaskstr, bool searchForNewDisulfides,
                                             DataSet_Coords& coords, Frame const& frameIn,
                                             CpptrajFile* outfile)
const
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
  if (coords.Top().SetupIntegerMask( cysmask )) return 1; 
  cysmask.MaskInfo();
  if (cysmask.None())
    mprintf("Warning: No cysteine sulfur atoms selected by %s\n", cysmaskstr.c_str());
  else {
    int nDisulfides = 0;
    double cut2 = disulfidecut * disulfidecut;
    // Try to find potential disulfide sites.
    // Keep track of which atoms will be part of disulfide bonds.
    std::vector<int> disulfidePartner( cysmask.Nselected(), -1 );
    // First, check for existing disulfides.
    for (AtomMask::const_iterator at1 = cysmask.begin(); at1 != cysmask.end(); ++at1)
    {
      for (AtomMask::const_iterator at2 = at1 + 1; at2 != cysmask.end(); ++at2)
      {
        if (coords.Top()[*at1].IsBondedTo(*at2)) {
          mprintf("\tExisting disulfide: %s to %s\n",
                  coords.Top().ResNameNumAtomNameNum(*at1).c_str(),
                  coords.Top().ResNameNumAtomNameNum(*at2).c_str());
          int idx1 = (int)(at1 - cysmask.begin());
          int idx2 = (int)(at2 - cysmask.begin());
          disulfidePartner[idx1] = idx2;
          disulfidePartner[idx2] = idx1;
        }
      }
    }
    // DEBUG - Print current array
    mprintf("DEBUG: Disulfide partner array after existing:\n");
    for (std::vector<int>::const_iterator it = disulfidePartner.begin(); it != disulfidePartner.end(); ++it)
    {
      mprintf("  S %i [%li]", cysmask[it-disulfidePartner.begin()]+1, it-disulfidePartner.begin());
      if (*it == -1)
        mprintf(" None.\n");
      else
        mprintf(" to S %i [%i]\n", cysmask[*it]+1, *it);
    }
    // Second, search for new disulfides from remaining sulfurs.
    if (searchForNewDisulfides) {
      // Only search with atoms that do not have an existing partner.
      std::vector<int> s_idxs; // Indices into cysmask/disulfidePartner
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

        for (std::vector<int>::const_iterator it1 = s_idxs.begin(); it1 != s_idxs.end(); ++it1)
        {
          int at1 = cysmask[*it1];
          for (std::vector<int>::const_iterator it2 = it1 + 1; it2 != s_idxs.end(); ++it2)
          {
            int at2 = cysmask[*it2];
            double r2 = DIST2_NoImage(frameIn.XYZ(at1), frameIn.XYZ(at2));
            if (r2 < cut2)
              D2.push_back( D2Pair(r2, IdxPair(*it1, *it2)) );
          }
        }
        std::sort(D2.begin(), D2.end());
        mprintf("DEBUG: Sorted S-S array:\n");
        for (D2Array::const_iterator it = D2.begin(); it != D2.end(); ++it)
        {
          int at1 = cysmask[it->second.first];
          int at2 = cysmask[it->second.second];
          mprintf("  %8i - %8i = %g Ang.\n", at1+1, at2+2, sqrt(it->first));
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
                    coords.Top().ResNameNumAtomNameNum(at1).c_str(),
                    coords.Top().ResNameNumAtomNameNum(at2).c_str(), sqrt(it->first));
            disulfidePartner[it->second.first ] = it->second.second;
            disulfidePartner[it->second.second] = it->second.first;
          }
        } // END loop over sorted distances
      } // END s_idxs not empty()
    } // END search for new disulfides
    // For each sulfur that has a disulfide partner, generate a bond command
    // and change residue name.
    for (std::vector<int>::const_iterator idx1 = disulfidePartner.begin(); idx1 != disulfidePartner.end(); ++idx1)
    {
      if (*idx1 != -1) {
        int at1 = cysmask[idx1-disulfidePartner.begin()];
        int at2 = cysmask[*idx1];
        if (at1 < at2) {
          nDisulfides++;
          LeapBond(at1, at2, coords.Top(), outfile);
        }
        ChangeResName(coords.TopPtr()->SetRes(coords.Top()[at1].ResNum()), newcysname);
      }
    }
    mprintf("\tDetected %i disulfide bonds.\n", nDisulfides);
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

  // Check that coords has no issues
  if (!coords.Top().AtomAltLoc().empty()) {
    // Must have only 1 atom alternate location
    char firstAltLoc = ' ';
    for (std::vector<char>::const_iterator altLocId = coords.Top().AtomAltLoc().begin();
                                           altLocId != coords.Top().AtomAltLoc().end();
                                         ++altLocId)
    {
      if (firstAltLoc == ' ') {
        // Find first non-blank alternate location ID
        if (*altLocId != ' ')
          firstAltLoc = *altLocId;
      } else if (*altLocId != ' ' && *altLocId != firstAltLoc) {
        mprinterr("Error: '%s' has atoms with multiple alternate location IDs, which\n"
                  "Error:  are not supported by LEaP. Use the 'keepaltloc <char>'\n"
                  "Error:  keyword for parm/trajin etc. to select which alternate\n"
                  "Error:  atom locations to keep.\n", coords.legend());
        return CpptrajState::ERR;
      }
    }
    mprintf("\t'%s' only contains atoms from alternate location ID '%c'\n", coords.legend(), firstAltLoc);
  }

  int tgtframe = argIn.getKeyInt("frame", 1) - 1;
  mprintf("\tUsing frame %i from COORDS set %s\n", tgtframe+1, coords.legend());
  if (tgtframe < 0 || tgtframe >= (int)coords.Size()) {
    mprinterr("Error: Frame is out of range.\n");
    return CpptrajState::ERR;
  }
  Frame frameIn = coords.AllocateFrame();
  coords.GetFrame(tgtframe, frameIn);

  std::string pdbout = argIn.GetStringKey("pdbout");
  if (!pdbout.empty())
    mprintf("\tPDB will be written to %s\n", pdbout.c_str());

  // Load PDB to glycam residue name map
  if (LoadGlycamPdbResMap( argIn.GetStringKey("resmapfile" ) )) {
    mprinterr("Error: PDB to glycam name map load failed.\n");
    return CpptrajState::ERR;
  }
  // DEBUG - print residue name map
  mprintf("\tResidue name map:\n");
  for (MapType::const_iterator mit = pdb_to_glycam_.begin(); mit != pdb_to_glycam_.end(); ++mit)
    mprintf("\t  %4s -> %c\n", *(mit->first), mit->second);

  leapunitname_ = argIn.GetStringKey("leapunitname", "m");
  mprintf("\tUsing leap unit name: %s\n", leapunitname_.c_str());

  bool prepare_sugars = !argIn.hasKey("nosugars");
  if (!prepare_sugars)
    mprintf("\tNot attempting to prepare sugars.\n");
  else
    mprintf("\tWill attempt to prepare sugars.\n");

  AtomMask sugarMask;
  std::string sugarmaskstr = argIn.GetStringKey("sugarmask");
  if (!sugarmaskstr.empty()) {
    if (!prepare_sugars) {
      mprinterr("Error: Cannot specify 'nosugars' and 'sugarmask'\n");
      return CpptrajState::ERR;
    }
    if (sugarMask.SetMaskString(sugarmaskstr)) {
      mprinterr("Error: Setting sugar mask string.\n");
      return CpptrajState::ERR;
    }
  } else if (prepare_sugars) {
    // No sugar mask specified; create one from names in pdb_to_glycam_ map.
    sugarmaskstr.assign(":");
    for (MapType::const_iterator mit = pdb_to_glycam_.begin(); mit != pdb_to_glycam_.end(); ++mit)
    {
      if (mit != pdb_to_glycam_.begin())
        sugarmaskstr.append(",");
      sugarmaskstr.append( mit->first.Truncated() );
    }
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
  mprintf("\tLEaP input containing 'loadpdb' and bond commands for disulfides,\n"
          "\t  sugars, etc will be written to '%s'\n", outfile->Filename().full());
  // Add the loadpdb command if we are writing a PDB file.
  if (!pdbout.empty())
    outfile->Printf("%s = loadpdb %s\n", leapunitname_.c_str(), pdbout.c_str());

  // Disulfide search
  if (!argIn.hasKey("nodisulfides")) {
    if (SearchForDisulfides( argIn.getKeyDouble("disulfidecut", 2.5),
                             argIn.GetStringKey("newcysname", "CYX"),
                             argIn.GetStringKey("cysmask", ":CYS@SG"),
                            !argIn.hasKey("existingdisulfides"),
                             coords, frameIn, outfile ))
    {
      mprinterr("Error: Disulfide search failed.\n");
      return CpptrajState::ERR;
    }
  } else {
    mprintf("\tNot searching for disulfides.\n");
  }

  // Prepare sugars
  if (prepare_sugars) {
    if (PrepareSugars(sugarMask, coords, frameIn, outfile)) {
      mprinterr("Error: Sugar preparation failed.\n");
      return CpptrajState::ERR;
    }
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

  if (!pdbout.empty()) {
    Trajout_Single PDB;
    if (PDB.InitTrajWrite( pdbout, "topresnum", State.DSL(), TrajectoryFile::PDBFILE)) {
      mprinterr("Error: Could not initialize output PDB\n");
      return CpptrajState::ERR;
    }
    if (PDB.SetupTrajWrite(coords.TopPtr(), coords.CoordsInfo(), 1)) {
      mprinterr("Error: Could not set up output PDB\n");
      return CpptrajState::ERR;
    }
    PDB.PrintInfo(1);
    PDB.WriteSingle(0, frameIn);
    PDB.EndTraj();
  }

  return CpptrajState::OK;
}
