#include "Exec_PrepareForLeap.h"
#include "CpptrajStdio.h"
#include "DistRoutines.h"
#include "CharMask.h"
#include "TorsionRoutines.h"
#include "Constants.h"
#include "CpptrajFile.h"
#include "Trajout_Single.h"
#include "DataSet_Coords_CRD.h"
#include <stack>
#include <cctype> // tolower
#include <algorithm> // sort

// ----- Sugar Class -----------------------------------------------------------
Exec_PrepareForLeap::Sugar::Sugar(int rn) :
  rnum_(rn),
  ring_oxygen_atom_(-1),
  anomeric_atom_(-1),
  ano_ref_atom_(-1)
{}

Exec_PrepareForLeap::Sugar::Sugar(int rn, int roa, int aa, int ara, std::vector<int> const& RA) :
  rnum_(rn),
  ring_oxygen_atom_(roa),
  anomeric_atom_(aa),
  ano_ref_atom_(ara),
  ring_atoms_(RA)
{}

void Exec_PrepareForLeap::Sugar::PrintInfo(Topology const& topIn) const {
  if (NotSet()) {
    mprintf("\t%s : Not Set.\n", topIn.TruncResNameNum(rnum_).c_str());
    //mprintf("DEBUG:\tres=%i %c\n", topIn.Res(rnum_).OriginalResNum(), topIn.Res(rnum_).ChainId());
  } else {
    mprintf("\t%s :\n", topIn.TruncResNameNum(rnum_).c_str());
    mprintf("\t\tRing O           : %s\n", topIn.TruncAtomNameNum(ring_oxygen_atom_).c_str());
    mprintf("\t\tAnomeric C       : %s\n", topIn.TruncAtomNameNum(anomeric_atom_).c_str());
    mprintf("\t\tAnomeric ref. C  : %s\n", topIn.TruncAtomNameNum(ano_ref_atom_).c_str());
    mprintf("\t\tNum ring atoms   : %u\n", NumRingAtoms());
    mprintf("\t\tNon-O Ring atoms :");
    for (std::vector<int>::const_iterator it = ring_atoms_.begin(); it != ring_atoms_.end(); ++it)
      mprintf(" %s", topIn.TruncAtomNameNum(*it).c_str());
    mprintf("\n");
  }
}

unsigned int Exec_PrepareForLeap::Sugar::NumRingAtoms() const {
  if (NotSet()) return 0;
  return ring_atoms_.size() + 1;
}

// -----------------------------------------------------------------------------

/** CONSTRUCTOR */
Exec_PrepareForLeap::Exec_PrepareForLeap() : Exec(COORDS),
  errorsAreFatal_(true),
  debug_(0)
{
  SetHidden(true);
}

// Exec_PrepareForLeap::Help()
void Exec_PrepareForLeap::Help() const
{
  mprintf("\tcrdset <coords set> [frame <#>] name <out coords set> [pdbout <pdbfile>]\n"
          "\t[leapunitname <unit>] [out <leap input file> [skiperrors]\n"
          "\t[nowat [watername <watername>] [noh] [keepaltloc <alt loc ID>]\n"
          "\t[stripmask <stripmask>] [solventresname <solventresname>]\n"
          "\t[{nohisdetect |\n"
          "\t  [nd1 <nd1>] [ne2 <ne2] [hisname <his>] [hiename <hie>]\n"
          "\t  [hidname <hid>] [hipname <hip]}]\n"
          "\t[{nodisulfides |\n"
          "\t  existingdisulfides |\n"
          "\t  [cysmask <cysmask>] [disulfidecut <cut>] [newcysname <name>]}]\n"
          "\t[{nosugars | sugarmask <sugarmask> [noc1search]}] [resmapfile <file>]\n"
          "\t[molmask <molmask> ...] [determinemolmask <mask>]\n"
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

/// Used to change residue name to nameIn
void Exec_PrepareForLeap::ChangeResName(Residue& res, NameType const& nameIn) const {
  if (res.Name() != nameIn) {
    if (debug_ > 0) mprintf("\t    Changing residue %s to %s\n", *(res.Name()), *nameIn);
    res.SetName( nameIn );
  }
}

/// Used to change atom name to nameIn
void Exec_PrepareForLeap::ChangeAtomName(Atom& atm, NameType const& nameIn) const {
  if (atm.Name() != nameIn) {
    if (debug_ > 0) mprintf("\t    Changing atom %s to %s\n", *(atm.Name()), *nameIn);
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
  //mprintf("\t  linkstr= '%s'\n", linkstr.c_str());
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
    mprintf("Warning: Could not determine link code for link atoms '%s'.\n", linkstr.c_str());
  return linkcode;
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
  std::string fname;
  if (fnameIn.empty()) {
    // Check CPPTRAJHOME
    const char* env = getenv("CPPTRAJHOME");
    if (env != 0) {
      fname.assign(env);
      fname += "/dat/PDB_ResidueNames.txt";
    }
    mprintf("Info: Parameter file path from CPPTRAJHOME variable: '%s'\n", fname.c_str());
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

/** Load reduced interal PDB to Glycam map. */
void Exec_PrepareForLeap::SetGlycamPdbResMap() {
  pdb_to_glycam_.insert( PairType("NAG", 'Y') );
  pdb_to_glycam_.insert( PairType("FUC", 'F') );
  pdb_to_glycam_.insert( PairType("GAL", 'L') );
  pdb_to_glycam_.insert( PairType("BMA", 'M') );
  pdb_to_glycam_.insert( PairType("MAN", 'M') );
  // TODO internal atom name map
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
  // Describe which section of the file we are in
  enum SectionType { PDB_RESMAP_SECTION = 0, PDB_ATOMMAP_SECTION, PDB_LINKAGE_RES_SECTION };
  SectionType section = PDB_RESMAP_SECTION;
  while ( (ptr = infile.NextLine()) != 0 ) {
    ArgList argline( ptr, " " );
    // Check for section change first
    if (argline.Nargs() < 1) {
      if (section == PDB_RESMAP_SECTION) {
        //mprintf("DEBUG: Section change.\n");
        section = PDB_ATOMMAP_SECTION;
      } else if (section == PDB_ATOMMAP_SECTION) {
        section = PDB_LINKAGE_RES_SECTION;
      }
    } else if (argline[0][0] != '#') {
      // Skipping comments, read sections
      if (section == PDB_RESMAP_SECTION) {
        // "<Name>" <glycam reschar> <pdb resname list>
        if (argline.Nargs() != 3) {
          mprinterr("Error: Expected only 3 columns in '%s' res map section, got %i\n",
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
      } else if (section == PDB_ATOMMAP_SECTION) {
        // <glycam reschar list> <PDB atomname to glycam atomname pair> ...
        if (argline.Nargs() < 2) {
          mprinterr("Error: Expected at least 2 columns in '%s' atom map section, got %i\n",
                    infile.Filename().full(), argline.Nargs());
          mprinterr("Error: %s\n", ptr);
          return 1;
        }
        // TODO handle glycam res names with > 1 char
        ArgList glycamnames( argline[0], "," );
        if (glycamnames.Nargs() < 1) {
          mprinterr("Error: No Glycam names found.\n");
          mprinterr("Error: %s\n", ptr);
          return 1;
        }
        int glycam_map_idx = (int)pdb_glycam_name_maps_.size();
        pdb_glycam_name_maps_.push_back(NameMapType());
        NameMapType& currentMap = pdb_glycam_name_maps_.back();
        for (int col = 1; col < argline.Nargs(); col++) {
          ArgList namepair( argline[col], "," );
          if (namepair.Nargs() != 2) {
            mprinterr("Error: Expected only 2 names for name pair, got %i\n", namepair.Nargs());
            mprinterr("Error: %s\n", ptr);
          }
          currentMap.insert( NamePairType(NameType(namepair[0]), NameType(namepair[1])) );
        } // END loop over name pair columns
        // Map will be for each glycam res
        for (ArgList::const_iterator gres = glycamnames.begin(); gres != glycamnames.end(); ++gres)
          glycam_res_idx_map_.insert( ResIdxPairType( (*gres)[0], glycam_map_idx ) );
      } else if (section == PDB_LINKAGE_RES_SECTION) {
        // <pdb linkage res name> <glycam linkage res name>
        if (argline.Nargs() != 2) {
          mprinterr("Error: Expected only 2 columns in '%s' linkage res map section, got %i\n",
                    infile.Filename().full(), argline.Nargs());
          mprinterr("Error: %s\n", ptr);
        }
        pdb_glycam_linkageRes_map_.insert( NamePairType(NameType(argline[0]),
                                                        NameType(argline[1])) );
      }
    } // END not comment
  } // END loop over file
  infile.CloseFile();

  return 0;
}

/** Determine torsion around non-ring stereocenter. */
int Exec_PrepareForLeap::CalcStereocenterTorsion(double& torsion,
                                                 int idx,
                                                 Topology const& topIn, Frame const& frameIn)
const
{
  // Assume part of a carbon chain, e.g.
  //      X
  //      |
  // C0 - C - C1
  //      |
  //      Y
  // where X is heavy atom and Y may or may not exist if it is hydrogen.
  // Try to order like so:
  // C0 - C - C1 - X
  int i0 = -1;
  int i1 = -1;
  int ix = -1;
  //int iy = -1;

  for (Atom::bond_iterator bat = topIn[idx].bondbegin(); bat != topIn[idx].bondend(); ++bat)
  {
    if ( topIn[*bat].Element() == Atom::CARBON ) {
      if (i0 == -1) {
        i0 = *bat;
      } else if (i1 == -1) {
        i1 = *bat;
      } else {
        mprinterr("Error: Too many carbons around stereocenter %s\n",
                  topIn.ResNameNumAtomNameNum(idx).c_str());
        return 1;
      }
    } else {
      if (ix == -1) {
        ix = *bat;
      } else {
        if ( topIn[*bat].Element() > topIn[ix].Element() ) {
          //iy = ix;
          ix = *bat;
        } //else
          //iy = ix;
      }
    }
  }
  if (i0 == -1) {mprinterr("Error: CalcStereocenterTorsion: Lower C is empty.\n"); return 1;}
  if (i1 == -1) {mprinterr("Error: CalcStereocenterTorsion: Higher C is empty.\n"); return 1;}
  if (ix == -1) {mprinterr("Error: CalcStereocenterTorsion: X substituent is empty.\n"); return 1;}
  torsion = Torsion(frameIn.XYZ(i0), frameIn.XYZ(idx), frameIn.XYZ(i1), frameIn.XYZ(ix));
  mprintf("\t  Stereocenter torsion: %s-%s-%s-%s = %f\n",
          *(topIn[i0].Name()), *(topIn[idx].Name()), *(topIn[i1].Name()), *(topIn[ix].Name()),
          torsion * Constants::RADDEG);
  return 0;
}
         

/** Determine torsion around anomeric reference carbon. */
int Exec_PrepareForLeap::CalcAnomericRefTorsion(double& torsion,
                                                int ano_ref_atom, int ring_oxygen_atom,
                                                Topology const& topIn, Frame const& frameIn,
                                                std::vector<bool> const& IsRingAtom)
const
{
  int ano_ref_atom_Y = -1;
  int ano_ref_atom_C = -1;
  // Get substituent of last ring C (e.g. C5) that is not part of the ring (e.g. C6)
  // Get substituent of last ring C (e.g. C5) that is part of the ring (e.g. C4)
  for ( Atom::bond_iterator bat = topIn[ano_ref_atom].bondbegin();
                            bat != topIn[ano_ref_atom].bondend();
                          ++bat )
  {
    if ( topIn[*bat].Element() == Atom::CARBON && !IsRingAtom[*bat] ) {
      if (ano_ref_atom_Y != -1) {
        mprinterr("Error: Two potential substituents for anomeric ref: %s and %s\n",
                  topIn.ResNameNumAtomNameNum(*bat).c_str(),
                  topIn.ResNameNumAtomNameNum(ano_ref_atom_Y).c_str());
        return 1;
      }
      ano_ref_atom_Y = *bat;
    } else if ( topIn[*bat].Element() == Atom::CARBON && IsRingAtom[*bat] ) {
      if (ano_ref_atom_C != -1) {
        mprinterr("Error: Two potential ring atoms bonded to anomeric ref: %s and %s\n",
                  topIn.ResNameNumAtomNameNum(*bat).c_str(),
                  topIn.ResNameNumAtomNameNum(ano_ref_atom_C).c_str());
        return 1;
      }
      ano_ref_atom_C = *bat;
    }
  }
  if (debug_ > 0)
    mprintf("\t  Anomeric ref carbon                   : %s\n",
            topIn.ResNameNumAtomNameNum(ano_ref_atom).c_str());
  if (ano_ref_atom_Y == -1) {
    mprinterr("Error: Anomeric reference Y substituent could not be identified.\n");
    return 1;
  }
  if (debug_ > 0)
    mprintf("\t  Anomeric reference substituent        : %s\n",
            topIn.ResNameNumAtomNameNum(ano_ref_atom_Y).c_str());
  if (ano_ref_atom_C == -1) {
    mprinterr("Error: Anomeric reference ring C previous ring atom could not be identified.\n");
    return 1;
  }
  if (debug_ > 0)
    mprintf("\t  Anomeric reference previous ring atom : %s\n",
            topIn.ResNameNumAtomNameNum(ano_ref_atom_C).c_str());
  torsion = Torsion( frameIn.XYZ(ano_ref_atom_C),   frameIn.XYZ(ano_ref_atom),
                     frameIn.XYZ(ano_ref_atom_Y),   frameIn.XYZ(ring_oxygen_atom) );
  if (debug_ > 0)
    mprintf("\t  Anomeric reference torsion            = %f\n", torsion * Constants::RADDEG);
  return 0;
}

/** Determine torsion around the anomeric carbon. */
int Exec_PrepareForLeap::CalcAnomericTorsion(double& torsion,
                                             int ring_oxygen_atom, int anomeric_atom,
                                             int rnum,
                                             Topology const& topIn, Frame const& frameIn,
                                             std::vector<bool> const& IsRingAtom)
const
{
  int anomeric_atom_X = -1;
  int anomeric_atom_C = -1;
  // Get the substituent of first ring C (e.g. C1) that is to a non-ring atom, non hydrogen 
  // Get the substituent of first ring C (e.g. C1) that is part of the ring (e.g. C2)
  for ( Atom::bond_iterator bat = topIn[anomeric_atom].bondbegin();
                            bat != topIn[anomeric_atom].bondend();
                          ++bat )
  {
    if ( topIn[*bat].Element() != Atom::HYDROGEN && !IsRingAtom[*bat] ) {
      if (anomeric_atom_X != -1) {
        // If there are two non-ring, non-hydrogen substituents, prioritize
        // the one that is part of this residue.
        bool bat_in_res = (topIn[*bat].ResNum() == rnum);
        bool X_in_res   = (topIn[anomeric_atom_X].ResNum() == rnum);
        if ( (bat_in_res && X_in_res) || (!bat_in_res && !X_in_res) ) {
          mprinterr("Error: Two potential substituents for anomeric carbon: %s and %s\n",
                    topIn.ResNameNumAtomNameNum(*bat).c_str(),
                    topIn.ResNameNumAtomNameNum(anomeric_atom_X).c_str());
          return 1;
        } else if (bat_in_res) {
          anomeric_atom_X = *bat;
        }
      } else
        anomeric_atom_X = *bat;
    } else if ( topIn[*bat].Element() == Atom::CARBON && IsRingAtom[*bat] ) {
      if (anomeric_atom_C != -1) {
        mprinterr("Error: Two potential ring carbons bonded to anomeric carbon: %s and %s\n",
                  topIn.ResNameNumAtomNameNum(*bat).c_str(),
                  topIn.ResNameNumAtomNameNum(anomeric_atom_C).c_str());
        return 1;
      }
      anomeric_atom_C = *bat;
    }
  }
  if (debug_ > 0)
    mprintf("\t  Anomeric carbon             : %s\n", topIn.ResNameNumAtomNameNum(anomeric_atom).c_str());
  if (anomeric_atom_X == -1) {
    // If the Cx (C1 substituent, usually a different residue) index is
    // not found this usually means missing inter-residue bond.
    // Alternatively, this could be an isolated sugar missing an -OH
    // group, so make this non-fatal.
    mprintf("Warning: Anomeric C non-ring substituent could not be identified.\n"
            "Warning: This can happen if the sugar is bonded to something that\n"
            "Warning:  is missing, e.g. a -OH group. In that case the coordinates\n"
            "Warning   for the missing atoms may need to be generated.\n");
    return -1;
  }
  if (debug_ > 0)
    mprintf("\t  Anomeric X substituent      : %s\n",
            topIn.ResNameNumAtomNameNum(anomeric_atom_X).c_str());
  if (anomeric_atom_C == -1) {
    mprinterr("Error: Next ring atom after anomeric C could not be identified.\n");
    return 1;
  }
  if (debug_ > 0)
    mprintf("\t  Anomeric C ring substituent : %s\n",
            topIn.ResNameNumAtomNameNum(anomeric_atom_C).c_str());
  torsion = Torsion( frameIn.XYZ(anomeric_atom_X), frameIn.XYZ(anomeric_atom),
                     frameIn.XYZ(anomeric_atom_C), frameIn.XYZ(ring_oxygen_atom) );
  if (debug_ > 0)
    mprintf("\t  Anomeric torsion            = %f\n", torsion * Constants::RADDEG);
  return 0;
}

/// Recursive function for finding and recording all carbons
static void Find_Carbons(int atm, Topology const& topIn, std::vector<bool>& Visited,
                         std::vector<int>& remainingChainCarbons)
{
  remainingChainCarbons.push_back( atm );
  Visited[atm] = true;
  // Follow all carbons bonded to this atom
  for (Atom::bond_iterator bat = topIn[atm].bondbegin(); bat != topIn[atm].bondend(); ++bat)
  {
    if (topIn[*bat].Element() == Atom::CARBON && !Visited[*bat]) {
      Find_Carbons( *bat, topIn, Visited, remainingChainCarbons );
    }
  }
}

/** Find remaining non-ring carbons in chain starting from ring end atom. */
int Exec_PrepareForLeap::FindRemainingChainCarbons(std::vector<int>& remainingChainCarbons,
                                                   int start_c, Topology const& topIn, int rnum,
                                                   std::vector<bool> const& IsRingAtom)
const
{
  Residue const& res = topIn.Res(rnum);
  std::vector<bool> Visited(topIn.Natom(), true);
  for (int at = res.FirstAtom(); at != res.LastAtom(); at++)
    if (!IsRingAtom[at])
      Visited[at] = false;

  for (Atom::bond_iterator bat = topIn[start_c].bondbegin();
                           bat != topIn[start_c].bondend();
                         ++bat)
  {
    if ( !Visited[*bat] && topIn[*bat].Element() == Atom::CARBON )
      Find_Carbons(*bat, topIn, Visited, remainingChainCarbons);
  }
  return 0;
}

/// Recursive function for following bonds of an atom to a target atom
static void FollowBonds(int atm, Topology const& topIn, int idx, std::vector<int>& ring_atoms, int tgt_atom, std::vector<bool>& Visited, bool& found)
{
  Visited[atm] = true;
  //for (int i = 0; i != idx; i++)
  //  mprintf("\t");
  //mprintf("At atom %s\n", topIn.ResNameNumAtomNameNum(atm).c_str());
  ring_atoms[idx] = atm;
  // Assume we have started at the target atom
  if (idx > 0 && atm == tgt_atom) {
    found = true;
    return;
  }
  // Follow all atoms bonded to this atom
  for (Atom::bond_iterator bat = topIn[atm].bondbegin(); bat != topIn[atm].bondend(); ++bat)
  {
    if (topIn[*bat].Element() == Atom::CARBON && !Visited[*bat]) {
      FollowBonds( *bat, topIn, idx+1, ring_atoms, tgt_atom, Visited, found );
      if (found) return;
    }
  }
}

/** Attempt to find any missing linkages to the anomeric carbon in sugar. */
int Exec_PrepareForLeap::FindSugarC1Linkages(Sugar const& sugar, Topology& topIn, Frame const& frameIn)
const
{
  if (sugar.NotSet()) {
    mprintf("Warning: Sugar %s is not set up. Skipping linkage detection.\n",
            topIn.TruncResNameNum(sugar.ResNum()).c_str());
    return 0; // TODO return 1?
  }
  int rnum1 = sugar.ResNum();
  int c_beg = sugar.AnomericAtom();
  Residue const& res1 = topIn.SetRes(rnum1);

  // residue first atom to residue first atom cutoff^2
  const double rescut2 = 64.0;
  // bond cutoff offset
  const double offset = 0.2;

  Atom::AtomicElementType a1Elt = topIn[c_beg].Element(); // Should always be C
  if (debug_ > 0)
    mprintf("DEBUG: Anomeric ring carbon: %s\n", topIn.ResNameNumAtomNameNum(c_beg).c_str());
  // Loop over other residues
  for (int rnum2 = 0; rnum2 < topIn.Nres(); rnum2++)
  {
    if (rnum2 != rnum1) {
      Residue const& res2 = topIn.Res(rnum2);
      // Ignore solvent residues
      if (res2.Name() != solventResName_) {
        int at1 = res1.FirstAtom();
        int at2 = res2.FirstAtom();
        // Initial residue-residue distance based on first atoms in each residue
        double dist2_1 = DIST2_NoImage( frameIn.XYZ(at1), frameIn.XYZ(at2) );
        if (dist2_1 < rescut2) {
          if (debug_ > 1)
            mprintf("DEBUG: %s to %s = %f\n",
                    topIn.TruncResNameNum(rnum1).c_str(), topIn.TruncResNameNum(rnum2).c_str(),
                    sqrt(dist2_1));
          // Do the rest of the atoms in res2 to the anomeric carbon
          for (; at2 != res2.LastAtom(); ++at2)
          {
            if (!topIn[c_beg].IsBondedTo(at2)) {
              double D2 = DIST2_NoImage( frameIn.XYZ(c_beg), frameIn.XYZ(at2) );
              Atom::AtomicElementType a2Elt = topIn[at2].Element();
              double cutoff2 = Atom::GetBondLength(a1Elt, a2Elt) + offset;
              cutoff2 *= cutoff2;
              if (D2 < cutoff2) {
                mprintf("\t  Adding bond between %s and %s\n",
                        topIn.ResNameNumAtomNameNum(c_beg).c_str(),
                        topIn.ResNameNumAtomNameNum(at2).c_str());
                topIn.AddBond(c_beg, at2);
              }
            }
          } // END loop over res2 atoms
        } // END res1-res2 distance cutoff
      } // END res2 is not solvent
    } // END res1 != res2
  } // END res2 loop over other residues

  return 0;
}

/** Identify sugar oxygen, anomeric and ref carbons, and ring atoms. */
Exec_PrepareForLeap::Sugar Exec_PrepareForLeap::IdSugarRing(int rnum, Topology const& topIn,
                                                            int& err)
{
  err = 0;
  Residue const& res = topIn.Res(rnum);
  // Try to identify the sugar ring. Potential starting atoms are oxygens
  // bonded to two carbon atoms. Also save potential stereocenter indices
  // (i.e. carbons bonded to 4 other atoms). Since input structure may not
  // have any hydrogens, count bonds to heavy atoms only, must make at
  // least 3 bonds (otherwise e.g. likely 2 hydrogens).
  std::vector<int> potentialRingStartAtoms;
  for (int at = res.FirstAtom(); at != res.LastAtom(); at++)
  {
    Atom const& currentAtom = topIn[at];
    if (currentAtom.Element() == Atom::OXYGEN) {
      if (currentAtom.Nbonds() == 2) {
        if ( topIn[currentAtom.Bond(0)].Element() == Atom::CARBON &&
             topIn[currentAtom.Bond(0)].ResNum() == rnum &&
             topIn[currentAtom.Bond(1)].Element() == Atom::CARBON &&
             topIn[currentAtom.Bond(1)].ResNum() == rnum )
        {
          potentialRingStartAtoms.push_back( at );
        }
      }
    }
  }
  // TODO handle case where multiple potential ring start atoms exist
  if (potentialRingStartAtoms.empty()) {
    mprintf("Warning: Ring oxygen could not be identified for %s\n",
            topIn.TruncResNameNum(rnum).c_str());
    resStat_[rnum] = SUGAR_MISSING_RING_O;
    return Sugar(rnum);
  } else if (potentialRingStartAtoms.size() > 1) {
    mprinterr("Error: Multiple potential ring start atoms:\n");
    for (std::vector<int>::const_iterator it = potentialRingStartAtoms.begin();
                                          it != potentialRingStartAtoms.end();
                                         ++it)
      mprinterr("Error:   %s\n", topIn.ResNameNumAtomNameNum(*it).c_str());
    err = 1;
    return Sugar(rnum);
  }

  // The anomeric carbon is the carbon that was part of the carbonyl group
  // in the straight chain. It is therefore typically the carbon with fewer
  // bonds to other carbons.
  int anomeric_atom = -1;    // e.g. C1
  // The anomeric reference carbon is the stereocenter farthest from the
  // anomeric carbon in the ring.
  int ano_ref_atom = -1;     // e.g. C5
  // Out of the potential ring start atoms, see which ones are actually
  // part of a ring. Potential ring start atoms only have 2 bonds,
  // each one to a carbon.
  int ring_oxygen_atom = -1; // e.g. O5
  // This will hold ring atoms, not including the ring oxygen
  std::vector<int> RA;
  for (std::vector<int>::const_iterator ringat = potentialRingStartAtoms.begin();
                                        ringat != potentialRingStartAtoms.end();
                                      ++ringat)
  {
    std::vector<bool> Visited( topIn.Natom(), true );
    for (int at = res.FirstAtom(); at != res.LastAtom(); at++)
      if (at != *ringat)
        Visited[at] = false;
    std::vector<int> ring_atoms( topIn.Res(rnum).NumAtoms(), -1 );
    bool ring_complete = false;
    // Since we have already established that *ringat is an oxygen bonded
    // to two carbons, just start at the first carbon to see if we can
    // get to the second carbon.
    int c_beg, c_end;
    if (topIn[*ringat].Bond(0) < topIn[*ringat].Bond(1)) {
      c_beg = topIn[*ringat].Bond(0);
      c_end = topIn[*ringat].Bond(1);
    } else {
      c_beg = topIn[*ringat].Bond(1);
      c_end = topIn[*ringat].Bond(0);
    }
    FollowBonds(c_beg, topIn, 0, ring_atoms,
                c_end, Visited, ring_complete);
    if (debug_ > 0)
      mprintf("DEBUG: Potential ring start atom %s, Ring complete = %i",
              topIn.ResNameNumAtomNameNum(*ringat).c_str(), (int)ring_complete);
    // TODO handle the case where multiple potential ring start atoms exist
    if (ring_complete) {
      ring_oxygen_atom = *ringat;
      // Try to ascertain which carbon might be the anomeric carbon (i.e. the
      // carbon that originally started the chain). Tie goes to lower index.
      int c_beg_bonds_to_C = 0;
      for (Atom::bond_iterator bat = topIn[c_beg].bondbegin(); bat != topIn[c_beg].bondend(); ++bat)
        if (topIn[*bat].Element() == Atom::CARBON)
          c_beg_bonds_to_C++;
      int c_end_bonds_to_C = 0;
      for (Atom::bond_iterator bat = topIn[c_end].bondbegin(); bat != topIn[c_end].bondend(); ++bat)
        if (topIn[*bat].Element() == Atom::CARBON)
          c_end_bonds_to_C++;
      if (debug_ > 0)
        mprintf("(%s bonds to C= %i, %s bonds to C = %i)", // DEBUG
                topIn.ResNameNumAtomNameNum(c_beg).c_str(), c_beg_bonds_to_C,
                topIn.ResNameNumAtomNameNum(c_end).c_str(), c_end_bonds_to_C);
      if (c_beg_bonds_to_C <= c_end_bonds_to_C) {
        anomeric_atom = c_beg;
        ano_ref_atom  = c_end;
      } else {
        anomeric_atom = c_end;
        ano_ref_atom = c_beg;
      }
      RA.clear();
      if (debug_ > 0) mprintf(" :"); // DEBUG
      for (std::vector<int>::const_iterator it = ring_atoms.begin(); it != ring_atoms.end(); ++it)
      {
        if (debug_ > 0) mprintf(" %i", *it + 1);
        if (*it == -1) break;
        RA.push_back( *it );
      }
      if (debug_ > 0) mprintf("\n"); // DEBUG
    }
  }
  if (RA.empty() || ring_oxygen_atom == -1) {
    mprinterr("Error: Sugar ring atoms could not be identified.\n");
    err = 1;
    return Sugar(rnum);
  }
  if (debug_ > 0)
    mprintf("\t  Ring oxygen         : %s\n", topIn.ResNameNumAtomNameNum(ring_oxygen_atom).c_str());
  return Sugar(rnum, ring_oxygen_atom, anomeric_atom, ano_ref_atom, RA);
}

/** Change PDB atom names in residue to glycam ones. */
int Exec_PrepareForLeap::ChangePdbAtomNamesToGlycam(char resChar, Residue const& res, Topology& topIn)
const
{
  // Get the appropriate map
  ResIdxMapType::const_iterator resIdxPair = glycam_res_idx_map_.find( resChar );
  if (resIdxPair == glycam_res_idx_map_.end()) {
    // No map needed for this residue
    //mprintf("DEBUG: No atom map for residue '%c'.\n", resChar);
    return 0;
  }
  NameMapType const& currentMap = pdb_glycam_name_maps_[resIdxPair->second];
  // Change PDB names to Glycam ones
  for (int at = res.FirstAtom(); at != res.LastAtom(); at++)
  {
    NameMapType::const_iterator namePair = currentMap.find( topIn[at].Name() );
    if (namePair != currentMap.end())
      ChangeAtomName( topIn.SetAtom(at), namePair->second );
  }
  return 0;
}

/** Attempt to identify sugar residue, form, and linkages. */
int Exec_PrepareForLeap::IdentifySugar(Sugar const& sugar, Topology& topIn,
                                       Frame const& frameIn, CharMask const& cmask,
                                       CpptrajFile* outfile, std::set<BondType>& sugarBondsToRemove)
{
  if (sugar.NotSet()) {
    mprintf("Warning: Sugar %s is not set up. Skipping sugar identification.\n",
            topIn.TruncResNameNum(sugar.ResNum()).c_str());
    return 0; // TODO return 1?
  }

  int rnum = sugar.ResNum();
  Residue& res = topIn.SetRes(rnum);
  // Try to ID the base sugar type from the input name.
  char resChar = ' ';

  MapType::const_iterator pdb_glycam = pdb_to_glycam_.find( res.Name() );
  if ( pdb_glycam == pdb_to_glycam_.end() ) { 
    mprinterr("Error: Could not identify sugar from residue name '%s'\n", *res.Name());
    return 1;
  }
  resChar = pdb_glycam->second;

  mprintf("\tSugar %s %i glycam name: %c\n", *res.Name(), rnum+1, resChar);
  if (debug_ > 0)
    mprintf("DEBUG:\tOriginal #= %i chain %c\n", res.OriginalResNum(), res.ChainId());

  // Change PDB names to Glycam ones
  if (ChangePdbAtomNamesToGlycam(resChar, res, topIn)) {
    mprinterr("Error: Changing PDB atom names to Glycam failed.\n");
    return 1;
  }

  // Create an array with all ring atoms set to true
  std::vector<bool> IsRingAtom;
  IsRingAtom.assign( topIn.Natom(), false );
  IsRingAtom[sugar.RingOxygenAtom()] = true;
  for (Sugar::const_iterator it = sugar.ringbegin(); it != sugar.ringend(); ++it)
    IsRingAtom[ *it ] = true;

  double t_c1 = 0;
  double t_c5 = 0;

  // Calculate torsion around anomeric carbon:
  int ret = CalcAnomericTorsion(t_c1, sugar.RingOxygenAtom(), sugar.AnomericAtom(),
                                rnum, topIn, frameIn, IsRingAtom);
  if (ret < 0) {
    // This means C1 X substituent missing; non-fatal.
    resStat_[rnum] = SUGAR_MISSING_C1X;
    return 0;
  } else if (ret > 0) {
    // Error
    return 1;
  }

  // Calculate torsion around anomeric reference:
  if (CalcAnomericRefTorsion(t_c5, sugar.AnomericRefAtom(), sugar.RingOxygenAtom(),
                             topIn, frameIn, IsRingAtom))
    return 1;

  // Determine alpha/beta 
  std::string formStr;
  bool c5up = (t_c5 > 0);
  bool c1up = (t_c1 > 0);
  if (c1up == c5up) {
    mprintf("\t  Beta form\n");
    formStr = "B";
  } else {
    mprintf("\t  Alpha form\n");
    formStr = "A";
  }

  // Find the rest of the carbons in the chain in order to find the
  // stereocenter with the highest index. Start from final ring carbon.
  int highest_stereocenter = -1;
  std::vector<int> remainingChainCarbons;
  if (FindRemainingChainCarbons(remainingChainCarbons, sugar.AnomericRefAtom(),
                                topIn, rnum, IsRingAtom))
    return 1;
  if (debug_ > 0) mprintf("\t  Remaining chain carbons:\n");
  for (std::vector<int>::const_iterator it = remainingChainCarbons.begin();
                                        it != remainingChainCarbons.end();
                                      ++it)
  {
    if (debug_ > 0) mprintf("\t\t%s", topIn.ResNameNumAtomNameNum(*it).c_str());
    // Count number of bonds to heavy atoms.
    Atom const& currentAtom = topIn[*it];
    int n_heavyat_bonds = 0;
    for (Atom::bond_iterator bat = currentAtom.bondbegin(); bat != currentAtom.bondend(); ++bat)
    {
      if ( topIn[*bat].Element() != Atom::HYDROGEN )
        ++n_heavyat_bonds;
    }
    if (n_heavyat_bonds == 3) {
      if (debug_ > 0) mprintf(" Potential stereocenter");
      if (*it > highest_stereocenter) // Is absolute index the best way to do this?
        highest_stereocenter = *it;
    }
    if (debug_ > 0) mprintf("\n");
  }
  if (highest_stereocenter == -1) {
    // This means that ano_ref_atom is the highest stereocenter.
    highest_stereocenter = sugar.AnomericRefAtom();
  }
  mprintf("\t  Highest stereocenter: %s\n", topIn.ResNameNumAtomNameNum(highest_stereocenter).c_str());
  // Determine D/L
  bool isDform;
  if (highest_stereocenter != sugar.AnomericRefAtom()) {
    double stereo_t = 0;
    if (CalcStereocenterTorsion(stereo_t, highest_stereocenter, topIn, frameIn))
      return 1;
    isDform = (stereo_t > 0);
  } else {
    isDform = c5up;
  }
  if (isDform) {
    mprintf("\t  D form\n");
  } else {
    mprintf("\t  L form\n");
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
    for (Atom::bond_iterator bat = topIn[at].bondbegin();
                             bat != topIn[at].bondend(); ++bat)
    {
      if (topIn[*bat].ResNum() != rnum) {
        if (!cmask.AtomInCharMask(*bat)) {
          mprintf("\t  Sugar %s bonded to non-sugar %s\n",
                  topIn.ResNameNumAtomNameNum(at).c_str(),
                  topIn.ResNameNumAtomNameNum(*bat).c_str());
          linkages.insert( topIn[at].Name() );
          bondsToRemove.push_back( BondType(at, *bat, -1) );
          // Check if this is a recognized linkage to non-sugar
          Residue& pres = topIn.SetRes( topIn[*bat].ResNum() );
          NameMapType::const_iterator lname = pdb_glycam_linkageRes_map_.find( pres.Name() );
          if (lname != pdb_glycam_linkageRes_map_.end()) {
            mprintf("DEBUG: Link residue name for %s found: %s\n", *(lname->first), *(lname->second));
            ChangeResName( pres, lname->second );
            resStat_[topIn[*bat].ResNum()] = VALIDATED;
          } else {
            mprintf("Warning: Unrecognized link residue %s, not modifying name.\n", *pres.Name());
            resStat_[topIn[*bat].ResNum()] = UNRECOGNIZED_SUGAR_LINKAGE;
          }
        } else {
          mprintf("\t  Sugar %s bonded to sugar %s\n",
                  topIn.ResNameNumAtomNameNum(at).c_str(),
                  topIn.ResNameNumAtomNameNum(*bat).c_str());
          linkages.insert( topIn[at].Name() );
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
  if (debug_ > 0) {
    mprintf("\t  Link atoms:");
    for (std::set<NameType>::const_iterator it = linkages.begin();
                                            it != linkages.end(); ++it)
      mprintf(" %4s", *(*it));
    mprintf("\n");
  }
  std::string linkcode = LinkageCode(resChar, linkages);
  if (debug_ > 0) mprintf("\t  Linkage code: %s\n", linkcode.c_str());
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
    LeapBond(bnd->A1(), bnd->A2(), topIn, outfile);
    topIn.RemoveBond(bnd->A1(), bnd->A2());
  }
  // Set new residue name
  NameType newResName( linkcode + std::string(1,resChar) + formStr );
  mprintf("\t  Changing %s to Glycam resname: %s\n", topIn.TruncResNameNum(rnum).c_str(), *newResName);
  ChangeResName(res, newResName);
  resStat_[rnum] = VALIDATED;
  return 0;
}

/** Prepare sugars for leap. */
int Exec_PrepareForLeap::PrepareSugars(AtomMask& sugarMask, Topology& topIn,
                                       Frame const& frameIn, CpptrajFile* outfile,
                                       bool findC1linkages)
{
  std::set<BondType> sugarBondsToRemove;
  mprintf("\tPreparing sugars selected by '%s'\n", sugarMask.MaskString());
  if (topIn.SetupIntegerMask( sugarMask )) return 1;
  sugarMask.MaskInfo();
  if (sugarMask.None())
    mprintf("Warning: No sugar atoms selected by %s\n", sugarMask.MaskString());
  else {
    CharMask cmask( sugarMask.ConvertToCharMask(), sugarMask.Nselected() );
    // Get sugar residue numbers
    std::vector<int> sugarResNums = topIn.ResnumsSelectedBy( sugarMask );
    // Try to identify sugar rings
    std::vector<Sugar> Sugars;
    Sugars.reserve( sugarResNums.size() );
    for (std::vector<int>::const_iterator rnum = sugarResNums.begin();
                                          rnum != sugarResNums.end(); ++rnum)
    {
      int err = 0;
      Sugars.push_back( IdSugarRing(*rnum, topIn, err) );
      if (err != 0) {
        if (errorsAreFatal_) {
          mprinterr("Error: Problem identifying sugar ring for %s\n",
                    topIn.TruncResNameNum(*rnum).c_str());
          return 1;
        } else
          mprintf("Warning: Problem identifying sugar ring for %s\n",
                    topIn.TruncResNameNum(*rnum).c_str());
      }
      Sugars.back().PrintInfo( topIn );
    }
    // For each sugar residue, try to fill in missing linkages
    if (findC1linkages) {
      mprintf("\tAttempting to identify missing linkages to sugar anomeric carbons.\n");
      for (std::vector<Sugar>::const_iterator sugar = Sugars.begin();
                                              sugar != Sugars.end(); ++sugar)
      {
        if (FindSugarC1Linkages(*sugar, topIn, frameIn)) {
          if (errorsAreFatal_)
            return 1;
          else
            mprintf("Warning: Finding anomeric C linkages for %s failed, skipping.\n",
                    topIn.TruncResNameNum( sugar->ResNum() ).c_str());
        }
      }
    } else {
      mprintf("\tNot attempting to identify missing linkages to sugar anomeric carbons.\n");
    }
    // For each sugar residue, see if it is bonded to a non-sugar residue.
    // If it is, remove that bond but record it.
    for (std::vector<Sugar>::const_iterator sugar = Sugars.begin();
                                            sugar != Sugars.end(); ++sugar)
    {
      //Residue const& Res = coords.Top().Res(*rnum);
      // See if we recognize this sugar.
      if (IdentifySugar(*sugar, topIn, frameIn, cmask, outfile, sugarBondsToRemove))
      {
        if (errorsAreFatal_)
          return 1;
        else
          mprintf("Warning: Preparation of sugar %s failed, skipping.\n",
                  topIn.TruncResNameNum( sugar->ResNum() ).c_str());
      }
    } // END loop over sugar residues
    // Remove bonds between sugars
    for (std::set<BondType>::const_iterator bnd = sugarBondsToRemove.begin();
                                            bnd != sugarBondsToRemove.end(); ++bnd)
    {
      LeapBond(bnd->A1(), bnd->A2(), topIn, outfile);
      topIn.RemoveBond(bnd->A1(), bnd->A2());
    }
  }
  // Bonds to sugars have been removed, so regenerate molecule info
  topIn.DetermineMolecules();
  return 0;
}

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
  std::vector<int> atomMolNum( topIn.Natom(), -1 );
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
              topIn.TruncResNameNum(lastRes).c_str());
      topIn.SetRes(lastRes).SetTerminal( true );
    }
  }
  return 0;
}

/** Search for disulfide bonds. */
int Exec_PrepareForLeap::SearchForDisulfides(double disulfidecut, std::string const& newcysnamestr,
                                             std::string const& cysmaskstr, bool searchForNewDisulfides,
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
        if (topIn[*at1].IsBondedTo(*at2)) {
          mprintf("\tExisting disulfide: %s to %s\n",
                  topIn.ResNameNumAtomNameNum(*at1).c_str(),
                  topIn.ResNameNumAtomNameNum(*at2).c_str());
          int idx1 = (int)(at1 - cysmask.begin());
          int idx2 = (int)(at2 - cysmask.begin());
          disulfidePartner[idx1] = idx2;
          disulfidePartner[idx2] = idx1;
        }
      }
    }
    // DEBUG - Print current array
    if (debug_ > 1) {
      mprintf("DEBUG: Disulfide partner array after existing:\n");
      for (std::vector<int>::const_iterator it = disulfidePartner.begin(); it != disulfidePartner.end(); ++it)
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
    for (std::vector<int>::const_iterator idx1 = disulfidePartner.begin(); idx1 != disulfidePartner.end(); ++idx1)
    {
      if (*idx1 != -1) {
        int at1 = cysmask[idx1-disulfidePartner.begin()];
        int at2 = cysmask[*idx1];
        if (at1 < at2) {
          nDisulfides++;
          LeapBond(at1, at2, topIn, outfile);
        }
        ChangeResName(topIn.SetRes(topIn[at1].ResNum()), newcysname);
        resStat_[topIn[at1].ResNum()] = VALIDATED;
      }
    }
    mprintf("\tDetected %i disulfide bonds.\n", nDisulfides);
  }
  return 0;
}

/** Modify coords according to user wishes. */
int Exec_PrepareForLeap::ModifyCoords( Topology& topIn, Frame& frameIn,
                                       bool remove_water,
                                       char altLocChar, std::string const& stripMask,
                                       std::string const& waterMask )
const
{
  // Create a mask denoting which atoms will be kept.
  std::vector<bool> atomsToKeep( topIn.Natom(), true );
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
/*
    unsigned int nRemoved = 0;
    for (Topology::mol_iterator mol = topIn.MolStart(); mol != topIn.MolEnd(); ++mol) {
      if (mol->IsSolvent()) {
        for (Unit::const_iterator seg = mol->MolUnit().segBegin();
                                  seg != mol->MolUnit().segEnd();
                                ++seg)
          for (int satm = seg->Begin(); satm < seg->End(); ++satm)
          {
            atomsToKeep[satm] = false;
            ++nRemoved;
          }
      }
    }
    if (nRemoved == 0)
      mprintf("\tNo solvent to remove.\n");
    else
      mprintf("\t# solvent removed: %u\n", nRemoved);
*/
  }
  //if (remove_h) {
  //  for (Topology::atom_iterator atom = topIn.begin(); atom != topIn.end(); ++atom) {
  //    if (atom->Element() == Atom::HYDROGEN)
  //      atomsToKeep[atom - topIn.begin()] = false;
  //  }
  //}
  if (altLocChar != '\0') {
    if (topIn.AtomAltLoc().empty()) {
      mprintf("\tNo alternate atom locations.\n");
    } else {
      for (int idx = 0; idx != topIn.Natom(); idx++) {
        if (topIn.AtomAltLoc()[idx] != ' ' &&
            topIn.AtomAltLoc()[idx] != altLocChar)
          atomsToKeep[idx] = false;
      }
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

/** Try to determine histidine protonation.
  * \param HisResNames Array containing final histidine residue names.
  * \param HisResIdxs Array containing residue indices of histidine residues (correspond to HisResNames).
  * \param topIn Input topology.
  * \param ND1 ND1 atom name.
  * \param NE2 NE2 atom name.
  * \param HisName PDB histidine name.
  * \param HieName Name for epsilon-protonated His.
  * \param HidName Name for delta-protonated His.
  * \param HipName Name for doubly-protonated His.
  */
int Exec_PrepareForLeap::DetermineHisProt( std::vector<NameType>& HisResNames,
                                           std::vector<int>& HisResIdxs,
                                           Topology const& topIn,
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
  std::vector<int> resIdxs = topIn.ResnumsSelectedBy( mask );

  HisResIdxs.clear();
  HisResIdxs.reserve( resIdxs.size() );
  HisResNames.clear();
  HisResNames.reserve( resIdxs.size() );
  for (std::vector<int>::const_iterator rnum = resIdxs.begin();
                                        rnum != resIdxs.end(); ++rnum)
  {
    if (debug_ > 1)
      mprintf("DEBUG: %s (%i) (%c)\n", topIn.TruncResNameNum(*rnum).c_str(), topIn.Res(*rnum).OriginalResNum(), topIn.Res(*rnum).ChainId());
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
      mprintf("Warning: Atom %s not found for %s; skipping residue.\n", *ND1, topIn.TruncResNameNum(*rnum).c_str());
      continue;
    }
    if (ne2idx == -1) {
      mprintf("Warning: Atom %s not found for %s; skipping residue,\n", *NE2, topIn.TruncResNameNum(*rnum).c_str());
      continue;
    }
    if (debug_ > 1)
      mprintf("DEBUG: %s nd1idx= %i ne2idx= %i\n",
              topIn.TruncResNameNum( *rnum ).c_str(), nd1idx+1, ne2idx+1);
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
      HisResNames.push_back( HipName );
      HisResIdxs.push_back( *rnum );
    } else if (nd1h > 0) {
      HisResNames.push_back( HidName );
      HisResIdxs.push_back( *rnum );
    } else if (ne2h > 0) {
      HisResNames.push_back( HieName );
      HisResIdxs.push_back( *rnum );
    }
    //else {
    //  // Default to epsilon
    //  mprintf("\tUsing default name '%s' for %s\n", *HieName, topIn.TruncResNameNum(*rnum).c_str());
    //  HisResNames.push_back( HieName );
    //}
  }
  if (!HisResIdxs.empty()) {
    mprintf("\tChanged histidine names:\n");
    for (unsigned int idx = 0; idx < HisResIdxs.size(); idx++)
      mprintf("\t\t%i %s\n", HisResIdxs[idx]+1, *HisResNames[idx]);
  } else
    mprintf("\tNo histidine names were changed.\n");
  return 0;
}

// Exec_PrepareForLeap::Execute()
Exec::RetType Exec_PrepareForLeap::Execute(CpptrajState& State, ArgList& argIn)
{
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

  std::string pdbout = argIn.GetStringKey("pdbout");
  if (!pdbout.empty())
    mprintf("\tPDB will be written to %s\n", pdbout.c_str());

  leapunitname_ = argIn.GetStringKey("leapunitname", "m");
  mprintf("\tUsing leap unit name: %s\n", leapunitname_.c_str());
  solventResName_ = argIn.GetStringKey("solventresname", "HOH");
  mprintf("\tSolvent residue name: %s\n", solventResName_.c_str());

  bool prepare_sugars = !argIn.hasKey("nosugars");
  if (!prepare_sugars)
    mprintf("\tNot attempting to prepare sugars.\n");
  else
    mprintf("\tWill attempt to prepare sugars.\n");

  // Deal with any coordinate modifications
  bool remove_water     = argIn.hasKey("nowat");
  std::string waterMask = argIn.GetStringKey("watername", ":" + solventResName_);
  bool remove_h         = argIn.hasKey("noh");
  std::string altLocArg = argIn.GetStringKey("keepaltloc");
  char altLocChar = '\0';
  if (!altLocArg.empty())
    altLocChar = altLocArg[0];
  std::string stripMask = argIn.GetStringKey("stripmask");

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
        if (altLocChar == '\0') {
          altLocChar = firstAltLoc;
          mprintf("Warning: '%s' has atoms with multiple alternate location IDs, which\n"
                  "Warning:  are not supported by LEaP. Keeping only '%c'.\n"
                  "Warning: To choose a specific location to keep use the 'keepaltloc <char>'\n"
                  "Warning:  keyword.\n", coords.legend(), altLocChar);
         }
         break;
      }
    }
  }

  if (remove_water)
    mprintf("\tRemoving solvent. Solvent mask= '%s'\n", waterMask.c_str());
  if (remove_h)
    mprintf("\tRemoving hydrogens.\n");
  if (altLocChar != '\0')
    mprintf("\tIf present, keeping only alternate atom locations denoted by '%c'\n", altLocChar);
  if (!stripMask.empty())
    mprintf("\tRemoving atoms in mask '%s'\n", stripMask.c_str());
  if (ModifyCoords( topIn, frameIn, remove_water, altLocChar, stripMask, waterMask ))
  {
    mprinterr("Error: Modification of '%s' failed.\n", coords.legend());
    return CpptrajState::ERR;
  }

  // Do histidine detection before H atoms are removed
  std::vector<int> HisResIdxs;
  std::vector<NameType> HisResNames;
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
    if (DetermineHisProt( HisResNames, HisResIdxs, topIn,
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

  // Each residue starts out unknown.
  resStat_.assign( topIn.Nres(), UNKNOWN );

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
  if (prepare_sugars) {
    if (LoadGlycamPdbResMap( argIn.GetStringKey("resmapfile" ) )) {
      mprinterr("Error: PDB to glycam name map load failed.\n");
      return CpptrajState::ERR;
    }
    mprintf("\t%zu entries in PDB to glycam name map.\n", pdb_to_glycam_.size());
    if (debug_ > 0) {
      // DEBUG - print residue name map
      mprintf("\tResidue name map:\n");
      for (MapType::const_iterator mit = pdb_to_glycam_.begin(); mit != pdb_to_glycam_.end(); ++mit)
        mprintf("\t  %4s -> %c\n", *(mit->first), mit->second);
      // DEBUG - print atom name maps
      mprintf("\tRes char to atom map index map:\n");
      for (ResIdxMapType::const_iterator mit = glycam_res_idx_map_.begin(); mit != glycam_res_idx_map_.end(); ++mit)
        mprintf("\t  %c -> %i\n", mit->first, mit->second);
      mprintf("\tAtom name maps:\n");
      for (std::vector<NameMapType>::const_iterator it = pdb_glycam_name_maps_.begin(); it != pdb_glycam_name_maps_.end(); ++it)
      {
        mprintf("\t  %li)", it - pdb_glycam_name_maps_.begin());
        for (NameMapType::const_iterator mit = it->begin(); mit != it->end(); ++mit)
          mprintf(" %s:%s", *(mit->first), *(mit->second));
        mprintf("\n");
      }
      // DEBUG - print linkage res map
      mprintf("\tLinkage res name map:\n");
      for (NameMapType::const_iterator mit = pdb_glycam_linkageRes_map_.begin(); mit != pdb_glycam_linkageRes_map_.end(); ++mit)
        mprintf("\t  %s -> %s\n", *(mit->first), *(mit->second));
    }
  }

  // Get sugar mask or default sugar mask
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
    if (PrepareSugars(sugarMask, topIn, frameIn, outfile, !argIn.hasKey("noc1search"))) {
      mprinterr("Error: Sugar preparation failed.\n");
      return CpptrajState::ERR;
    }
  }

  // Count any solvent molecules
  if (!remove_water) {
    NameType solvName(solventResName_);
    unsigned int nsolvent = 0;
    for (Topology::res_iterator res = topIn.ResStart(); res != topIn.ResEnd(); ++res) {
      if ( res->Name() == solvName) {
        nsolvent++;
        resStat_[res-topIn.ResStart()] = VALIDATED;
      }
    }
    if (nsolvent > 0) mprintf("\t%u solvent residues.\n", nsolvent);
  }

  // Residue validation.
  //mprintf("\tResidues with potential problems:\n");
  static const char* msg = "Potential problem: ";
  for (ResStatArray::iterator it = resStat_.begin(); it != resStat_.end(); ++it)
  {
    //if ( *it == VALIDATED )
    //  mprintf("\t\t%s VALIDATED\n", topIn.TruncResNameNum(it-resStat_.begin()).c_str());
    //else
    //  mprintf("\t\t%s UNKNOWN\n", topIn.TruncResNameNum(it-resStat_.begin()).c_str());
    if ( *it == UNKNOWN ) {
      SetType::const_iterator pname = pdb_res_names_.find( topIn.Res(it-resStat_.begin()).Name() );
      if (pname == pdb_res_names_.end())
        mprintf("\t%s%s is an unrecognized name and may not have parameters.\n",
                msg, topIn.TruncResNameNum(it-resStat_.begin()).c_str());
      else
        *it = VALIDATED;
    } else if ( *it == UNRECOGNIZED_SUGAR_LINKAGE ) {
        mprintf("\t%s%s is linked to a sugar but has no sugar-linkage form.\n",
                msg, topIn.TruncResNameNum(it-resStat_.begin()).c_str());
    } else if ( *it == SUGAR_MISSING_C1X ) {
        mprintf("\t%s%s Sugar is missing anomeric carbon substituent and cannot be identified.\n",
                msg, topIn.TruncResNameNum(it-resStat_.begin()).c_str());
    } else if ( *it == SUGAR_MISSING_RING_O ) {
        mprintf("\t%s%s Sugar is missing ring oxygen and cannot be identified.\n",
                msg, topIn.TruncResNameNum(it-resStat_.begin()).c_str());
    }
  }

  // Change HIS res names
  for (unsigned int idx = 0; idx != HisResIdxs.size(); idx++)
    ChangeResName( topIn.SetRes(HisResIdxs[idx]), HisResNames[idx] );

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
        topIn.TruncResNameNum(lastRes).c_str());
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
    if (PDB.InitTrajWrite( pdbout, "topresnum", State.DSL(), TrajectoryFile::PDBFILE)) {
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

  return CpptrajState::OK;
}
