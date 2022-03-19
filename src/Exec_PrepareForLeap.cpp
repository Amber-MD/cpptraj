#include "Exec_PrepareForLeap.h"
#include "CharMask.h"
#include "Chirality.h"
#include "Constants.h"
#include "CpptrajFile.h"
#include "CpptrajStdio.h"
#include "DataSet_Coords_CRD.h"
#include "DistRoutines.h"
#include "LeapInterface.h"
#include "ParmFile.h"
#include "TorsionRoutines.h"
#include "Trajout_Single.h"
#include "StringRoutines.h" // integerToString
#include <stack>
#include <cctype> // tolower
#include <algorithm> // sort

const char* Exec_PrepareForLeap::ringstr_[] = {"pyranose", "furanose", "?"};
const char* Exec_PrepareForLeap::formstr_[] = {"alpha", "beta", "?"};
const char* Exec_PrepareForLeap::chirstr_[] = {"D", "L", "?"};

// ===== Sugar Class ===========================================================
/** CONSTRUCTOR - Incomplete setup; set anomeric atom as residue first atom
  *               so that ResNum() works.
  */
Exec_PrepareForLeap::Sugar::Sugar(StatType status, int firstat) :
  stat_(status),
  ring_oxygen_atom_(-1),
  anomeric_atom_(firstat),
  ano_ref_atom_(-1),
  highest_stereocenter_(-1),
  ringType_(UNKNOWN_RING),
  isMissingAtoms_(false)
{}

/** CONSTRUCTOR - Partial setup; ring O, anomeric atom, ring and chain atoms. */
Exec_PrepareForLeap::Sugar::Sugar(StatType status, int roa, int aa,
                                  Iarray const& RA, Iarray const& CA) :
  stat_(status),
  ring_oxygen_atom_(roa),
  anomeric_atom_(aa),
  ano_ref_atom_(-1),
  highest_stereocenter_(-1),
  ringType_(UNKNOWN_RING),
  isMissingAtoms_(true),
  ring_atoms_(RA),
  chain_atoms_(CA)
{
  if (ring_oxygen_atom_ != -1) {
    if (RA.size() == 5)
      ringType_ = PYRANOSE;
    else if (RA.size() == 4)
      ringType_ = FURANOSE;
  }
}

/** CONSTRUCTOR - Set ring atom indices and ring type. */
Exec_PrepareForLeap::Sugar::Sugar(int roa, int aa, int ara, int hs,
                                  Iarray const& RA, Iarray const& CA, bool isMissing) :
  stat_(SETUP_OK),
  ring_oxygen_atom_(roa),
  anomeric_atom_(aa),
  ano_ref_atom_(ara),
  highest_stereocenter_(hs),
  ringType_(UNKNOWN_RING),
  isMissingAtoms_(isMissing),
  ring_atoms_(RA),
  chain_atoms_(CA)
{
  if (ring_oxygen_atom_ != -1) {
    if (RA.size() == 5)
      ringType_ = PYRANOSE;
    else if (RA.size() == 4)
      ringType_ = FURANOSE;
  }
}

/** \return Residue number based on the anomeric carbon atom index. */
int Exec_PrepareForLeap::Sugar::ResNum(Topology const& topIn) const {
  return topIn[anomeric_atom_].ResNum();
}

/** Correspond to Sugar::StatType */
const char* Exec_PrepareForLeap::Sugar::StatTypeStr_[] = {
  "OK", "Missing ring oxygen", "Multiple ring oxygens", "Missing chain atoms",
  "Missing anomeric reference", "Missing configurational carbon" };

/** Print info about the sugar to STDOUT. */
void Exec_PrepareForLeap::Sugar::PrintInfo(Topology const& topIn) const {
  if (NotSet()) {
    mprintf("\t%s : Not Set. %s\n", topIn.TruncResNameOnumId(ResNum(topIn)).c_str(),
            StatTypeStr_[stat_]);
  } else {
    mprintf("\t%s : %s\n", topIn.TruncResNameOnumId(ResNum(topIn)).c_str(),
            StatTypeStr_[stat_]);
  }
  //mprintf("\t%s :\n", topIn.TruncResNameOnumId(ResNum(topIn)).c_str());
  if (isMissingAtoms_)
    mprintf("\t\tIs missing atoms.\n");
  if (ring_oxygen_atom_ != -1)
    mprintf("\t\tRing O           : %s\n", topIn.TruncAtomNameNum(ring_oxygen_atom_).c_str());
  if (anomeric_atom_ != -1)
    mprintf("\t\tAnomeric C       : %s\n", topIn.TruncAtomNameNum(anomeric_atom_).c_str());
  if (ano_ref_atom_ != -1)
    mprintf("\t\tAnomeric ref. C  : %s\n", topIn.TruncAtomNameNum(ano_ref_atom_).c_str());
  if (highest_stereocenter_ != -1)
    mprintf("\t\tConfig. C        : %s\n", topIn.TruncAtomNameNum(highest_stereocenter_).c_str());
  mprintf("\t\tNum ring atoms   : %u\n", NumRingAtoms());
  static const char* RingPrefixStr[] = { "Hexo", "Pento", "Other" };
  int rpsidx = 2;
  if (chain_atoms_.size() == 6)
    rpsidx = 0;
  else if (chain_atoms_.size() == 5)
    rpsidx = 1;
  static const char* RingTypeStr[] = { "Pyranose", "Furanose", "Unknown" };
  mprintf("\t\tRing Type        : %s%s\n", RingPrefixStr[rpsidx], RingTypeStr[ringType_]);
  mprintf("\t\tNon-O Ring atoms :");
  for (Iarray::const_iterator it = ring_atoms_.begin(); it != ring_atoms_.end(); ++it)
    mprintf(" %s", topIn.TruncAtomNameNum(*it).c_str());
  mprintf("\n\t\tChain atoms      :");
  for (Iarray::const_iterator it = chain_atoms_.begin(); it != chain_atoms_.end(); ++it)
    mprintf(" %s", topIn.TruncAtomNameNum(*it).c_str());
  mprintf("\n");
}

/** \return Number of ring atoms (including oxygen) */
unsigned int Exec_PrepareForLeap::Sugar::NumRingAtoms() const {
  if (NotSet()) return 0;
  return ring_atoms_.size() + 1;
}

/// \return what oldidx should be according to atomMap
inline static int find_new_idx(int oldidx, std::vector<int> const& atomMap, int at0, int at1) {
  for (int newidx = at0; newidx != at1; newidx++)
    if (atomMap[newidx] == oldidx)
      return newidx;
  return -1;
}

/** Remap internal indices according to given map. */
void Exec_PrepareForLeap::Sugar::RemapIndices(Iarray const& atomMap, int at0, int at1) {
  // Always try the anomeric atom
  anomeric_atom_ = find_new_idx(anomeric_atom_, atomMap, at0, at1);
  //mprintf("DEBUG: Ring O old = %i ring O new %i\n", ring_oxygen_atom_+1, atomMap[ring_oxygen_atom_]+1);
  if (ring_oxygen_atom_ != -1)
    ring_oxygen_atom_     = find_new_idx(ring_oxygen_atom_, atomMap, at0, at1);
  if (ano_ref_atom_ != -1)
    ano_ref_atom_         = find_new_idx(ano_ref_atom_, atomMap, at0, at1);
  if (highest_stereocenter_ != -1)
    highest_stereocenter_ = find_new_idx(highest_stereocenter_, atomMap, at0, at1);
  for (Iarray::iterator it = ring_atoms_.begin(); it != ring_atoms_.end(); ++it)
    *it = find_new_idx(*it, atomMap, at0, at1);
  for (Iarray::iterator it = chain_atoms_.begin(); it != chain_atoms_.end(); ++it)
    *it = find_new_idx(*it, atomMap, at0, at1);
}

// ===== SugarToken class ======================================================
/** CONSTRUCTOR */
Exec_PrepareForLeap::SugarToken::SugarToken() :
  form_(UNKNOWN_FORM),
  chir_(UNKNOWN_CHIR),
  ring_(UNKNOWN_RING)
{}

/** CONSTRUCTOR - name, glycam code, form, chirality, ring type */
Exec_PrepareForLeap::SugarToken::SugarToken(std::string const& fn, std::string const& gc,
                                            FormTypeEnum ft, ChirTypeEnum ct, RingTypeEnum rt) :
  name_(fn),
  glycamCode_(gc),
  form_(ft),
  chir_(ct),
  ring_(rt)
{}

/** CONSTRUCTOR - ring type */
Exec_PrepareForLeap::SugarToken::SugarToken(RingTypeEnum rt) :
  form_(UNKNOWN_FORM),
  chir_(UNKNOWN_CHIR),
  ring_(rt)
{}

/** \return String containing name, glycam code, and form-chirality-ring type. */
std::string Exec_PrepareForLeap::SugarToken::InfoStr() const {
  return std::string("\"" + name_ + "\" " + glycamCode_ + " " +
                     std::string(formstr_[form_]) + "-" +
                     std::string(chirstr_[chir_]) + "-" +
                     std::string(ringstr_[ring_]));
}

/** Set up from line: <res> <code> <form> <chir> <ring> <name>
  *                   0     1      2      3      4      5
  * \return Residue name <res>
  */
std::string Exec_PrepareForLeap::SugarToken::SetFromLine(ArgList const& line) {
  const char* lineIn = line.ArgLine();
  if (line.Nargs() != 6) {
    mprinterr("Error: Expected 6 columns, got %i\n"
              "Error: %s\n", line.Nargs(), lineIn);
    return std::string("");
  }
  name_ = line[5];
  glycamCode_ = line[1];
  if (line[2] == "A")
    form_ = ALPHA;
  else if (line[2] == "B")
    form_ = BETA;
  else {
    mprinterr("Error: Unrecognized anomer type: %s\n"
              "Error: Line: %s\n", line[2].c_str(), lineIn);
    return std::string("");
  }
  if (line[3] == "D")
    chir_ = IS_D;
  else if (line[3] == "L")
    chir_ = IS_L;
  else {
    mprinterr("Error: Unrecognized configuration: %s\n"
              "Error: Line: %s\n", line[3].c_str(), lineIn);
    return std::string("");
  }
  if (line[4] == "P")
    ring_ = PYRANOSE;
  else if (line[4] == "F")
    ring_ = FURANOSE;
  else {
    mprinterr("Error: Unrecognized ring: %s\n"
              "Error: Line: %s\n", line[4].c_str(), lineIn);
    return std::string("");
  }

  return line[0];
}
// =============================================================================

/** CONSTRUCTOR */
Exec_PrepareForLeap::Exec_PrepareForLeap() : Exec(COORDS),
  errorsAreFatal_(true),
  hasGlycam_(false),
  useSugarName_(false),
  debug_(0)
{
  SetHidden(false);
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

/** Load reduced interal PDB to Glycam map. */
void Exec_PrepareForLeap::SetGlycamPdbResMap() {
  pdb_to_glycam_.insert( PairType("64K",
    SugarToken("alpha-D-arabinopyranose", "A", ALPHA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("AHR",
    SugarToken("alpha-L-arabinofuranose", "A", ALPHA, IS_L, FURANOSE)) );
  pdb_to_glycam_.insert( PairType("ARA",
    SugarToken("alpha-L-arabinopyranose", "A", ALPHA, IS_L, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("ARB",
    SugarToken("beta-L-arabinopyranose", "A", BETA, IS_L, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("AXR",
    SugarToken("methyl alpha-D-arabinofuranoside", "A", ALPHA, IS_D, FURANOSE)) );
  pdb_to_glycam_.insert( PairType("BXX",
    SugarToken("beta-D-arabinofuranose", "A", BETA, IS_D, FURANOSE)) );
  pdb_to_glycam_.insert( PairType("BXY",
    SugarToken("alpha-D-arabinofuranose", "A", ALPHA, IS_D, FURANOSE)) );
  pdb_to_glycam_.insert( PairType("FUB",
    SugarToken("beta-L-arabinofuranose", "A", BETA, IS_L, FURANOSE)) );
  pdb_to_glycam_.insert( PairType("SEJ",
    SugarToken("beta-D-arabinopyranose", "A", BETA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("LDY",
    SugarToken("alpha-D-lyxopyranose", "D", ALPHA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("Z4W",
    SugarToken("beta-D-lyxopyranose", "D", BETA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("0MK",
    SugarToken("beta-L-ribopyranose", "R", BETA, IS_L, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("32O",
    SugarToken("beta-L-ribofuranose", "R", BETA, IS_L, FURANOSE)) );
  pdb_to_glycam_.insert( PairType("BDR",
    SugarToken("beta-D-ribofuranose", "R", BETA, IS_D, FURANOSE)) );
  pdb_to_glycam_.insert( PairType("RIB",
    SugarToken("alpha-D-ribofuranose", "R", ALPHA, IS_D, FURANOSE)) );
  pdb_to_glycam_.insert( PairType("RIP",
    SugarToken("beta-D-ribopyranose", "R", BETA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("YYM",
    SugarToken("alpha-D-ribopyranose", "R", ALPHA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("Z6J",
    SugarToken("alpha-L-ribofuranose", "R", ALPHA, IS_L, FURANOSE)) );
  pdb_to_glycam_.insert( PairType("HSY",
    SugarToken("alpha-L-xylopyranose", "X", ALPHA, IS_L, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("HSZ",
    SugarToken("beta-D-xylopyranose", "X", BETA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("LXC",
    SugarToken("beta-L-xylopyranose", "X", BETA, IS_L, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("XYP",
    SugarToken("beta-D-xylopyranose", "X", BETA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("XYS",
    SugarToken("alpha-D-xylopyranose", "X", ALPHA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("XYZ",
    SugarToken("beta-D-xylofuranose", "X", BETA, IS_D, FURANOSE)) );
  pdb_to_glycam_.insert( PairType("AFD",
    SugarToken("alpha-D-allopyranose", "N", ALPHA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("ALL",
    SugarToken("beta-D-allopyranose", "N", BETA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("VDS",
    SugarToken("beta-D-allofuranose", "N", BETA, IS_D, FURANOSE)) );
  pdb_to_glycam_.insert( PairType("VDV",
    SugarToken("alpha-D-allofuranose", "N", ALPHA, IS_D, FURANOSE)) );
  pdb_to_glycam_.insert( PairType("WOO",
    SugarToken("beta-L-allopyranose", "N", BETA, IS_L, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("Z2D",
    SugarToken("alpha-L-allopyranose", "N", ALPHA, IS_L, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("3MK",
    SugarToken("beta-L-altropyranose", "E", BETA, IS_L, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("SHD",
    SugarToken("alpha-D-altropyranose", "E", ALPHA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("Z6H",
    SugarToken("alpha-L-altropyranose", "E", ALPHA, IS_L, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("GAL",
    SugarToken("beta-D-galactopyranose", "L", BETA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("GIV",
    SugarToken("beta-L-galactopyranose", "L", BETA, IS_L, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("GLA",
    SugarToken("alpha-D-galactopyranose", "L", ALPHA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("GXL",
    SugarToken("alpha-L-galactopyranose", "L", ALPHA, IS_L, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("GZL",
    SugarToken("beta-D-galactofuranose", "L", BETA, IS_D, FURANOSE)) );
  pdb_to_glycam_.insert( PairType("BGC",
    SugarToken("beta-D-glucopyranose", "G", BETA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("GLC",
    SugarToken("alpha-D-glucopyranose", "G", ALPHA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("GU4",
    SugarToken("2,3,4,6-tetra-O-sulfonato-alpha-D-glucopyranose", "G", ALPHA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("MGL",
    SugarToken("methyl beta-D-glucopyranoside", "G", BETA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("Z8T",
    SugarToken("beta-L-glucopyranose", "G", BETA, IS_L, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("4GL",
    SugarToken("alpha-D-gulopyranose", "K", ALPHA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("GL0",
    SugarToken("beta-D-gulopyranose", "K", BETA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("GUP",
    SugarToken("alpha-L-gulopyranose", "K", ALPHA, IS_L, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("Z8H",
    SugarToken("beta-L-gulopyranose", "K", BETA, IS_L, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("BMA",
    SugarToken("beta-D-mannopyranose", "M", BETA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("MAN",
    SugarToken("alpha-D-mannopyranose", "M", ALPHA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("A5C",
    SugarToken("alpha-L-talofuranose", "T", ALPHA, IS_L, FURANOSE)) );
  pdb_to_glycam_.insert( PairType("SDY",
    SugarToken("beta-D-talopyranose", "T", BETA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("ZEE",
    SugarToken("beta-L-talopyranose", "T", BETA, IS_L, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("BDF",
    SugarToken("beta-D-fructopyranose", "C", BETA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("FRU",
    SugarToken("beta-D-fructofuranose", "C", BETA, IS_D, FURANOSE)) );
  pdb_to_glycam_.insert( PairType("LFR",
    SugarToken("beta-L-fructofuranose", "C", BETA, IS_L, FURANOSE)) );
  pdb_to_glycam_.insert( PairType("YYJ",
    SugarToken("1,3,4,6-tetra-O-sulfo-beta-D-fructofuranose", "C", BETA, IS_D, FURANOSE)) );
  pdb_to_glycam_.insert( PairType("Z9N",
    SugarToken("alpha-D-fructofuranose", "C", ALPHA, IS_D, FURANOSE)) );
  pdb_to_glycam_.insert( PairType("PSV",
    SugarToken("alpha-D-psicofuranose", "P", ALPHA, IS_D, FURANOSE)) );
  pdb_to_glycam_.insert( PairType("SF6",
    SugarToken("alpha-L-psicofuranose", "P", ALPHA, IS_L, FURANOSE)) );
  pdb_to_glycam_.insert( PairType("SF9",
    SugarToken("beta-L-psicofuranose", "P", BETA, IS_L, FURANOSE)) );
  pdb_to_glycam_.insert( PairType("TTV",
    SugarToken("beta-D-psicofuranose", "P", BETA, IS_D, FURANOSE)) );
  pdb_to_glycam_.insert( PairType("SOE",
    SugarToken("alpha-L-sorbopyranose", "B", ALPHA, IS_L, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("UEA",
    SugarToken("beta-D-sorbofuranose", "B", BETA, IS_D, FURANOSE)) );
  pdb_to_glycam_.insert( PairType("T6T",
    SugarToken("alpha-D-tagatopyranose", "J", ALPHA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("FCA",
    SugarToken("alpha-D-fucopyranose", "F", ALPHA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("FCB",
    SugarToken("beta-D-fucopyranose", "F", BETA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("FUC",
    SugarToken("alpha-L-fucopyranose", "F", ALPHA, IS_L, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("FUL",
    SugarToken("beta-L-fucopyranose", "F", BETA, IS_L, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("GYE",
    SugarToken("beta-D-fucofuranose", "F", BETA, IS_D, FURANOSE)) );
  pdb_to_glycam_.insert( PairType("MXY",
    SugarToken("2-O-methyl-beta-L-fucopyranose", "F", BETA, IS_L, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("MXZ",
    SugarToken("2-O-methyl-alpha-L-fucopyranose", "F", ALPHA, IS_L, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("G6D",
    SugarToken("alpha-D-quinovopyranose", "Q", ALPHA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("YYK",
    SugarToken("beta-D-quinovopyranose", "Q", BETA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("RAM",
    SugarToken("alpha-L-rhamnopyranose", "H", ALPHA, IS_L, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("RM4",
    SugarToken("beta-L-rhamnopyranose", "H", BETA, IS_L, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("XXR",
    SugarToken("alpha-D-rhamnopyranose", "H", ALPHA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("ADA",
    SugarToken("alpha-D-galactopyranuronic acid", "O", ALPHA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("GTR",
    SugarToken("beta-D-galactopyranuronic acid", "O", BETA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("BDP",
    SugarToken("beta-D-glucopyranuronic acid", "Z", BETA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("GCU",
    SugarToken("alpha-D-glucopyranuronic acid", "Z", ALPHA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("GCV",
    SugarToken("4-O-methyl-alpha-D-glucopyranuronic acid", "Z", ALPHA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("IDR",
    SugarToken("alpha-L-idopyranuronic acid", "U", ALPHA, IS_L, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("IDS",
    SugarToken("2-O-sulfo-alpha-L-idopyranuronic acid", "U", ALPHA, IS_L, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("A2G",
    SugarToken("2-acetamido-2-deoxy-alpha-D-galactopyranose", "V", ALPHA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("ASG",
    SugarToken("2-acetamido-2-deoxy-4-O-sulfo-beta-D-galactopyranose", "V", BETA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("NG6",
    SugarToken("2-acetamido-2-deoxy-6-O-sulfo-beta-D-galactopyranose", "V", BETA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("NGA",
    SugarToken("2-acetamido-2-deoxy-beta-D-galactopyranose", "V", BETA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("YYQ",
    SugarToken("2-acetamido-2-deoxy-alpha-L-galactopyranose", "V", ALPHA, IS_L, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("NAG",
    SugarToken("2-acetamido-2-deoxy-beta-D-glucopyranose", "Y", BETA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("NDG",
    SugarToken("2-acetamido-2-deoxy-alpha-D-glucopyranose", "Y", ALPHA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("NGZ",
    SugarToken("2-acetamido-2-deoxy-alpha-L-glucopyranose", "Y", ALPHA, IS_L, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("BM3",
    SugarToken("2-acetamido-2-deoxy-alpha-D-mannopyranose", "W", ALPHA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("BM7",
    SugarToken("2-acetamido-2-deoxy-beta-D-mannopyranose", "W", BETA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("SIA",
    SugarToken("N-acetyl-alpha-neuraminic acid", "S", ALPHA, IS_D, PYRANOSE)) );
  pdb_to_glycam_.insert( PairType("SLB",
    SugarToken("N-acetyl-beta-neuraminic acid", "S", BETA, IS_D, PYRANOSE)) );
  // PDB to glycam atom name maps
  // 0 - V,W,Y C7,C2N  O7,O2N  C8,CME
  pdb_glycam_name_maps_.push_back(NameMapType());
  pdb_glycam_name_maps_A_.push_back(NameMapType());
  pdb_glycam_name_maps_B_.push_back(NameMapType());
  pdb_glycam_name_maps_.back().insert(NamePairType("C7","C2N"));
  pdb_glycam_name_maps_.back().insert(NamePairType("O7","O2N"));
  pdb_glycam_name_maps_.back().insert(NamePairType("C8","CME"));
  glycam_res_idx_map_.insert( ResIdxPairType( "V", 0 ) );
  glycam_res_idx_map_.insert( ResIdxPairType( "W", 0 ) );
  glycam_res_idx_map_.insert( ResIdxPairType( "Y", 0 ) );
  // 1 - S     C10,C5N O10,O5N C11,CME
  pdb_glycam_name_maps_.push_back(NameMapType());
  pdb_glycam_name_maps_A_.push_back(NameMapType());
  pdb_glycam_name_maps_B_.push_back(NameMapType());
  pdb_glycam_name_maps_.back().insert(NamePairType("C10","C5N"));
  pdb_glycam_name_maps_.back().insert(NamePairType("O10","O5N"));
  pdb_glycam_name_maps_.back().insert(NamePairType("C11","CME"));
  glycam_res_idx_map_.insert( ResIdxPairType( "S", 1 ) );
  // 2 - H     C6,C6M,B
  pdb_glycam_name_maps_.push_back(NameMapType());
  pdb_glycam_name_maps_A_.push_back(NameMapType());
  pdb_glycam_name_maps_B_.push_back(NameMapType());
  pdb_glycam_name_maps_B_.back().insert(NamePairType("C6","C6M"));
  glycam_res_idx_map_.insert( ResIdxPairType( "H", 2 ) );
  // PDB to glycam linkage residue name maps
  pdb_glycam_linkageRes_map_.insert( NamePairType("SER", "OLS") );
  pdb_glycam_linkageRes_map_.insert( NamePairType("THR", "OLT") );
  pdb_glycam_linkageRes_map_.insert( NamePairType("HYP", "OLP") );
  pdb_glycam_linkageRes_map_.insert( NamePairType("ASN", "NLN") );
}

/** Load PDB to Glycam residue map from file. */
int Exec_PrepareForLeap::LoadGlycamPdbResMap(std::string const& fnameIn)
{
  std::string fname = fnameIn;
  if (fnameIn.empty()) {
    // Check CPPTRAJHOME
    const char* env = getenv("CPPTRAJHOME");
    if (env != 0) {
      fname.assign(env);
      fname += "/dat/Carbohydrate_PDB_Glycam_Names.txt";
      mprintf("Info: Parameter file path from CPPTRAJHOME variable: '%s'\n", fname.c_str());
    } else {
      // Check AMBERHOME
      env = getenv("AMBERHOME");
      if (env != 0) {
        fname.assign(env);
        fname += "/AmberTools/src/cpptraj/dat/Carbohydrate_PDB_Glycam_Names.txt";
        mprintf("Info: Parameter file path from AMBERHOME variable: '%s'\n", fname.c_str());
      }
    }
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
        // OLD: "<Name>" <glycam reschar> <pdb resname list>
        // NEW: <res> <glycam> <form> <chirality> <ring type> <full name>
        SugarToken sToken;
        std::string sResName = sToken.SetFromLine(argline);
        if (sResName.empty()) {
          mprinterr("Error: Could not parse residue map section of '%s'\n",
                    infile.Filename().full());
          return 1;
        }
        std::pair<MapType::iterator,bool> ret =
          pdb_to_glycam_.insert( PairType(sResName, sToken) );
        if (!ret.second) {
          mprinterr("Error: Duplicate residue name '%s' in residue map section of '%s'\n",
                    sResName.c_str(), infile.Filename().full());
          return 1;
        }
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
        pdb_glycam_name_maps_A_.push_back(NameMapType());
        pdb_glycam_name_maps_B_.push_back(NameMapType());
        NameMapType& currentMap  = pdb_glycam_name_maps_.back();
        NameMapType& currentMapA = pdb_glycam_name_maps_A_.back();
        NameMapType& currentMapB = pdb_glycam_name_maps_B_.back();
        for (int col = 1; col < argline.Nargs(); col++) {
          ArgList namepair( argline[col], "," );
          NameMapType* currentMapPtr = &currentMap;
          if (namepair.Nargs() == 3) {
            // This name mapping is for a particular anomeric form
            if (namepair[2] == "A")
              currentMapPtr = &currentMapA;
            else if (namepair[2] == "B")
              currentMapPtr = &currentMapB;
            else {
              mprinterr("Error: For name pair, third arg should only be A or B: %s\n", ptr);
              return 1;
            }
          } else if (namepair.Nargs() != 2) {
            mprinterr("Error: Expected only 2 names for name pair, got %i\n", namepair.Nargs());
            mprinterr("Error: %s\n", ptr);
            return 1;
          }
          currentMapPtr->insert( NamePairType(NameType(namepair[0]), NameType(namepair[1])) );
        } // END loop over name pair columns
        // Map will be for each glycam res
        for (ArgList::const_iterator gres = glycamnames.begin(); gres != glycamnames.end(); ++gres)
          glycam_res_idx_map_.insert( ResIdxPairType( *gres, glycam_map_idx ) );
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

// -----------------------------------------------------------------------------
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

/// \return True if the given tgt atom is in the given array
static inline bool AtomIsInArray(std::vector<int> const& RingAtoms, int tgt)
{
  for (std::vector<int>::const_iterator it = RingAtoms.begin(); it != RingAtoms.end(); ++it)
    if (*it == tgt) return true;
  return false;
}

/// \return Position of given tgt atom in the array (if it is in the given array)
static inline int AtomIdxInArray(std::vector<int> const& ChainAtoms, int tgt)
{
  for (std::vector<int>::const_iterator it = ChainAtoms.begin(); it != ChainAtoms.end(); ++it)
    if (*it == tgt)
      return (int)(it - ChainAtoms.begin());
  return -1;
}

/** Find remaining non-ring carbons in chain starting from ring end atom. */
int Exec_PrepareForLeap::FindRemainingChainCarbons(Iarray& remainingChainCarbons,
                                                   int start_c, Topology const& topIn, int rnum,
                                                   Iarray const& RingAtoms)
const
{
  Residue const& res = topIn.Res(rnum);
  std::vector<bool> Visited(topIn.Natom(), true);
  for (int at = res.FirstAtom(); at != res.LastAtom(); at++)
    if (!AtomIsInArray(RingAtoms, at))
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
  int rnum = topIn[atm].ResNum();
  //for (int i = 0; i != idx; i++) // DEBUG
  //  mprintf("\t"); // DBEUG
  //mprintf("At atom %s\n", topIn.ResNameNumAtomNameNum(atm).c_str()); // DEBUG
  ring_atoms[idx] = atm;
  // Assume we have started at the target atom
  if (idx > 0 && atm == tgt_atom) {
    found = true;
    return;
  }
  // Follow all atoms bonded to this atom
  for (Atom::bond_iterator bat = topIn[atm].bondbegin(); bat != topIn[atm].bondend(); ++bat)
  {
    if (topIn[*bat].ResNum() == rnum &&
        topIn[*bat].Element() == Atom::CARBON &&
        !Visited[*bat])
    {
      FollowBonds( *bat, topIn, idx+1, ring_atoms, tgt_atom, Visited, found );
      if (found) return;
    }
  }
}

// -----------------------------------------------
/** Determine torsion around the anomeric carbon. */
int Exec_PrepareForLeap::CalcAnomericTorsion(double& torsion,
                                             int anomeric_atom, int ring_oxygen_atom,
                                             int rnum,
                                             Iarray const& RingAtoms,
                                             Topology const& topIn, Frame const& frameIn)
const
{
  if (debug_ > 0) {
    mprintf("\t  Anomeric carbon             : %s\n", topIn.ResNameNumAtomNameNum(anomeric_atom).c_str());
    mprintf("\t  Ring oxygen atom            : %s\n", topIn.ResNameNumAtomNameNum(ring_oxygen_atom).c_str());
  }
  int anomeric_atom_X = -1;
  int anomeric_atom_C = -1;
  // By definition the anomeric atom should be the first ring atom TODO catch size==1?
  anomeric_atom_C = RingAtoms[1];
  if (anomeric_atom_C == -1) {
    mprinterr("Error: Next ring atom after anomeric C could not be identified.\n");
    return 1;
  }
  if (debug_ > 0)
    mprintf("\t  Anomeric C ring substituent : %s\n",
            topIn.ResNameNumAtomNameNum(anomeric_atom_C).c_str());
  // Get the substituent of the anomeric C (e.g. C1) that is a non-ring atom, non hydrogen 
  for ( Atom::bond_iterator bat = topIn[anomeric_atom].bondbegin();
                            bat != topIn[anomeric_atom].bondend();
                          ++bat )
  {
    if ( *bat != ring_oxygen_atom &&
         topIn[*bat].Element() != Atom::HYDROGEN &&
         !AtomIsInArray(RingAtoms, *bat) )
    {
      if (anomeric_atom_X != -1) {
        // If there are two non-ring, non-hydrogen substituents, prioritize
        // the one that is part of this residue.
        bool bat_in_res = (topIn[*bat].ResNum() == rnum);
        bool X_in_res   = (topIn[anomeric_atom_X].ResNum() == rnum);
        if ( (bat_in_res && X_in_res) || (!bat_in_res && !X_in_res) ) {
          // Both in or both out of residue. Choose based on atomic number.
          if (topIn[*bat].AtomicNumber() == topIn[anomeric_atom_X].AtomicNumber()) {
            mprinterr("Error: Two potential substituents for anomeric carbon: %s and %s\n",
                      topIn.ResNameNumAtomNameNum(*bat).c_str(),
                      topIn.ResNameNumAtomNameNum(anomeric_atom_X).c_str());
            return 1;
          } else if (topIn[*bat].AtomicNumber() > topIn[anomeric_atom_X].AtomicNumber()) {
            anomeric_atom_X = *bat;
          }
        } else if (bat_in_res) {
          anomeric_atom_X = *bat;
        }
      } else
        anomeric_atom_X = *bat;
    }
  }

  if (anomeric_atom_X == -1) {
    // If the Cx (C1 substituent, usually a different residue) index is
    // not found this usually means missing inter-residue bond.
    // Alternatively, this could be an isolated sugar missing an -OH
    // group, so make this non-fatal.
    mprintf("Warning: Anomeric C non-ring substituent could not be identified.\n");
    //        "Warning: This can happen if the sugar is bonded to something that\n"
    //        "Warning:  is missing, e.g. a -OH group. In that case the coordinates\n"
    //        "Warning   for the missing atoms may need to be generated.\n");
    return -1;
  }
  if (debug_ > 0)
    mprintf("\t  Anomeric X substituent      : %s\n",
            topIn.ResNameNumAtomNameNum(anomeric_atom_X).c_str());

  torsion = Torsion( frameIn.XYZ(ring_oxygen_atom), frameIn.XYZ(anomeric_atom),
                     frameIn.XYZ(anomeric_atom_C), frameIn.XYZ(anomeric_atom_X) );
  if (debug_ > 0)
    mprintf("DEBUG: Anomeric torsion %s-%s-%s-%s= %f\n",
            *(topIn[ring_oxygen_atom].Name()),
            *(topIn[anomeric_atom].Name()),
            *(topIn[anomeric_atom_C].Name()),
            *(topIn[anomeric_atom_X].Name()),
            torsion * Constants::RADDEG);
  return 0;
}

/** Determine torsion around anomeric reference carbon. */
int Exec_PrepareForLeap::CalcAnomericRefTorsion(double& torsion,
                                                int ano_ref_atom, int ring_oxygen_atom,
                                                int ring_end_atom, Iarray const& RingAtoms,
                                                Topology const& topIn, Frame const& frameIn)
const
{
  if (debug_ > 0)
    mprintf("\t  Anomeric ref carbon                   : %s\n",
            topIn.ResNameNumAtomNameNum(ano_ref_atom).c_str());
  //      ano_ref_atom_Y
  //           |
  //      ano_ref_atom
  //       |        |
  // ano_ref_atom_0 ano_ref_atom_1
  int ano_ref_atom_Y = -1;
  int ano_ref_atom_0 = -1;
  int ano_ref_atom_1 = -1;
  // This will be the index of the anomeric atom in the RingAtoms array
  int ar_index = -1;
  // Find ring atom that precedes the anomeric reference atom TODO catch size==1?
  for (unsigned int idx = 1; idx != RingAtoms.size(); idx++) {
    if (RingAtoms[idx] == ano_ref_atom) {
        ar_index = (int)idx;
        ano_ref_atom_0 = RingAtoms[idx-1];
        break;
    }
  }
  if (ano_ref_atom_0 == -1) {
    mprinterr("Error: Anomeric reference ring C previous ring atom could not be identified.\n");
    return 1;
  }
  if (debug_ > 0)
    mprintf("\t  Anomeric reference previous ring atom : %s\n",
            topIn.ResNameNumAtomNameNum(ano_ref_atom_0).c_str());
  // If the anomeric reference atom is the ring end atom then ano_ref_atom_1
  // is the ring oxygen.
  if (ano_ref_atom == ring_end_atom) {
    ano_ref_atom_1 = ring_oxygen_atom;
  } else {
    // Anomeric reference atom is somewhere before the ring end atom.
    ano_ref_atom_1 = RingAtoms[ar_index+1];
  }
  if (debug_ > 0)
    mprintf("\t  Anomeric reference next atom          : %s\n",
            topIn.ResNameNumAtomNameNum(ano_ref_atom_1).c_str());
  // Get non-hydrogen substituent of anomeric ref (e.g. C5) that 
  // is not part of the ring (e.g. C6).
  for ( Atom::bond_iterator bat = topIn[ano_ref_atom].bondbegin();
                            bat != topIn[ano_ref_atom].bondend();
                          ++bat )
  {
    if ( *bat != ring_oxygen_atom &&
         topIn[*bat].Element() != Atom::HYDROGEN &&
         !AtomIsInArray(RingAtoms, *bat) )
    {
      if (ano_ref_atom_Y != -1) {
        mprinterr("Error: Two potential non-ring substituents for anomeric ref: %s and %s\n",
                  topIn.ResNameNumAtomNameNum(*bat).c_str(),
                  topIn.ResNameNumAtomNameNum(ano_ref_atom_Y).c_str());
        return 1;
      }
      ano_ref_atom_Y = *bat;
    }
  }
  if (ano_ref_atom_Y == -1) {
    mprinterr("Error: Anomeric reference Y substituent could not be identified.\n");
    return 1;
  }
  if (debug_ > 0)
    mprintf("\t  Anomeric reference substituent        : %s\n",
            topIn.ResNameNumAtomNameNum(ano_ref_atom_Y).c_str());

  torsion = Torsion( frameIn.XYZ(ano_ref_atom_0),   frameIn.XYZ(ano_ref_atom),
                     frameIn.XYZ(ano_ref_atom_1),   frameIn.XYZ(ano_ref_atom_Y) );
  if (debug_ > 0)
    mprintf("DEBUG: Anomeric reference torsion %s-%s-%s-%s= %f\n",
            *(topIn[ano_ref_atom_0].Name()),
            *(topIn[ano_ref_atom].Name()),
            *(topIn[ano_ref_atom_1].Name()),
            *(topIn[ano_ref_atom_Y].Name()),
            torsion * Constants::RADDEG);
  return 0;
}

/** Determine torsion around the configurational carbon. 
  * Calculate torsion around config. carbon C as:
  *   C0-C-Z-C1
  * where C0 is the carbon preceding C in the chain, C1 is the carbon
  * after C in the chain, and Z is the non-hydrogen substituent of C
  * with highest priority. Do it this way to be consistent with how
  * CalcAnomericRefTorsion orders the atoms.
  */
int Exec_PrepareForLeap::CalcConfigCarbonTorsion(double& torsion, int config_carbon,
                                                 Iarray const& carbon_chain,
                                                 Topology const& topIn,
                                                 Frame const& frameIn)
const
{
  int atom_c0 = -1;
  int atom_c1 = -1;
  int atom_z  = -1;
  // Get c0 and c1
  int c_idx = AtomIdxInArray(carbon_chain, config_carbon);
  if (c_idx < 1) {
    mprinterr("Error: Could not determine carbon before config. C '%s'\n",
              topIn.ResNameNumAtomNameNum(config_carbon).c_str());
    return 1;
  }
  atom_c0 = carbon_chain[c_idx-1];
  if ((unsigned int)c_idx+1 >= carbon_chain.size()) {
    mprinterr("Error: Could not determine carbon after config. C '%s'\n",
              topIn.ResNameNumAtomNameNum(config_carbon).c_str());
    return 1;
  }
  atom_c1 = carbon_chain[c_idx+1];

  for (Atom::bond_iterator bat = topIn[config_carbon].bondbegin();
                           bat != topIn[config_carbon].bondend(); ++bat)
  {
    if (topIn[*bat].Element() != Atom::HYDROGEN &&
        !AtomIsInArray(carbon_chain, *bat))
    {
      if (atom_z == -1)
        atom_z = *bat;
      else if (topIn[*bat].AtomicNumber() > topIn[atom_z].AtomicNumber())
        atom_z = *bat;
    }
  }
  if (atom_z == -1) {
    mprinterr("Error: Could not determine substituent for config. C '%s'\n",
              topIn.ResNameNumAtomNameNum(config_carbon).c_str());
    return 1;
  }

  torsion = Torsion( frameIn.XYZ(atom_c0),
                     frameIn.XYZ(config_carbon),
                     frameIn.XYZ(atom_z),
                     frameIn.XYZ(atom_c1) );
  if (debug_ > 0)
    mprintf("DEBUG: Config. C torsion %s-%s-%s-%s= %f\n",
            *(topIn[atom_c0].Name()),
            *(topIn[config_carbon].Name()),
            *(topIn[atom_z].Name()),
            *(topIn[atom_c1].Name()),
            torsion*Constants::RADDEG);
  return 0;
}

// -----------------------------------------------
/** \return true if residue is on the list of residues with missing heteroatoms. */
static bool res_missing_het(int rnum, Topology const& topIn) {
  Residue const& res = topIn.Res(rnum);
  //mprintf("DEBUG: TGT %s %i '%c' '%c'\n", *(res.Name()), res.OriginalResNum(), res.Icode(), res.ChainId());
  for (std::vector<Residue>::const_iterator het = topIn.MissingHet().begin();
                                            het != topIn.MissingHet().end(); ++het)
  {
    //mprintf("DEBUG: HET %s %i '%c' '%c'\n", *(het->Name()), het->OriginalResNum(), het->Icode(), het->ChainId());
    if ( het->OriginalResNum() == res.OriginalResNum() &&
         het->Icode()          == res.Icode() &&
         het->ChainId()        == res.ChainId() &&
         het->Name()           == res.Name() )
      return true;
  }
  return false;
}

/** Identify sugar oxygen, anomeric and ref carbons, and ring atoms. */
Exec_PrepareForLeap::Sugar Exec_PrepareForLeap::IdSugarRing(int rnum, Topology const& topIn)
const
{
  Residue const& res = topIn.Res(rnum);
  bool residue_missing_atoms = res_missing_het(rnum, topIn);
  if (residue_missing_atoms)
    mprintf("\tResidue '%s' is missing atoms.\n", topIn.TruncResNameOnumId(rnum).c_str());

  // Determine candidates for ring oxygen atoms. 
  Iarray potentialRingStartAtoms;
  // This array will hold the 2 flanking carbons for potential ring oxygen atoms.
  std::vector<std::pair<int,int>> ringAtomCarbons;
  for (int at = res.FirstAtom(); at != res.LastAtom(); at++)
  {
    Atom const& currentAtom = topIn[at];
    // Try to identify the sugar ring oxygen. Candidate atoms are oxygens
    // bonded to two carbon atoms in the same residue.
    if (currentAtom.Element() == Atom::OXYGEN) {
      if (currentAtom.Nbonds() > 1) {
        Iarray c_atoms;
        for (Atom::bond_iterator bat = currentAtom.bondbegin();
                                 bat != currentAtom.bondend(); ++bat)
        {
          if (topIn[*bat].Element() == Atom::CARBON &&
              topIn[*bat].ResNum() == rnum)
            c_atoms.push_back( *bat );
        }
        if ( c_atoms.size() == 2 ) {
          potentialRingStartAtoms.push_back( at );
          ringAtomCarbons.push_back( std::pair<int,int>(c_atoms[0], c_atoms[1]) );
        }
      }
    }
  }

  if (potentialRingStartAtoms.empty()) {
    mprintf("Warning: Ring oxygen could not be identified for %s\n",
            topIn.TruncResNameOnumId(rnum).c_str());
    return Sugar(Sugar::MISSING_O, res.FirstAtom());
  }

  // Use the previously-set up AtomMap to help determine stereocenters
  std::vector<bool> atomIsChiral;
  atomIsChiral.reserve( res.NumAtoms() );
  // Since we cannot be certain there will be hydrogens, cannot rely
  // on the AtomMap chiral designations (which assumes hydrogens).
  // Make a chiral center a carbon with at least 3 bonds, all must be
  // to different kinds of atoms.
  int resat = res.FirstAtom();
  for (int iat = 0; iat != res.NumAtoms(); iat++, resat++)
  {
    bool chiral = false;
    if (topIn[resat].Element() == Atom::CARBON && topIn[resat].Nbonds() > 2) {
      if (debug_ > 0)
        mprintf("DEBUG: Atom '%s' potential chiral\n", topIn.TruncResAtomNameNum(resat).c_str());
      chiral = true;
      for (Atom::bond_iterator bat1 = topIn[resat].bondbegin();
                               bat1 != topIn[resat].bondend(); ++bat1)
      {
        std::string unique1 = myMap_[*bat1].Unique();
        for (Atom::bond_iterator bat2 = bat1 + 1; bat2 != topIn[resat].bondend(); ++bat2)
        {
          std::string unique2 = myMap_[*bat2].Unique();
          if (unique1 == unique2) {
            // At least two of the atoms bonded to this atom look the same. Not chiral.
            if (debug_ > 1)
              mprintf("DEBUG: unique strings match %s='%s' %s='%s'\n", *(topIn[*bat1].Name()), unique1.c_str(), *(topIn[*bat2].Name()), unique2.c_str());
            chiral = false;
            break;
          }
        } // END inner loop over bonded atoms
        if (!chiral) break;
      } // END outer loop over bonded atoms
    } // END atom is carbon with > 2 bonds
    
    atomIsChiral.push_back( chiral );
    if (debug_ > 0)
      mprintf("DEBUG: Atom '%s' isChiral= %i\n",
              topIn.TruncResAtomNameNum(resat).c_str(),
              (int)atomIsChiral.back());
  }

  // Using the potential ring start atoms, see if we can actually complete
  // a ring. Identify important atoms as well.
  // This will indicate ring direction
//  int ring_direction = 0;
  // Ring end atom is the last atom in the ring
  int ring_end_atom = -1;
  // The anomeric carbon is the carbon that was part of the carbonyl group
  // in the straight chain. It is therefore typically the carbon with fewer
  // bonds to other carbons.
  int anomeric_atom = -1;    // e.g. C1
  // The anomeric reference carbon is the stereocenter farthest from the
  // anomeric carbon in the ring.
  int ano_ref_atom = -1;     // e.g. C5
  // Ring oxygen atom
  int ring_oxygen_atom = -1; // e.g. O5
  // This will hold the index of the highest stereocenter, e.g. C5
  int highest_stereocenter = -1;
  // This will hold ring atoms, not including the ring oxygen.
  std::vector<Iarray> Ring_Atoms;
  // This will hold carbon chain atoms starting from the anomeric carbon
  Iarray carbon_chain;

  // Out of the potential ring start atoms, see which ones are actually
  // part of a ring. Potential ring start atoms only have 2 bonds,
  // each one to a carbon.
  int ring_atom_idx = -1;
  std::vector<std::pair<int,int>>::iterator catoms = ringAtomCarbons.begin();
  for (Iarray::const_iterator ringat = potentialRingStartAtoms.begin();
                              ringat != potentialRingStartAtoms.end();
                            ++ringat, ++catoms)
  {
    if (debug_ > 0)
      mprintf("DEBUG: Ring start '%s'\n", topIn.ResNameNumAtomNameNum(*ringat).c_str());
    // Mark all atoms as visited except this residue (minus the ring start).
    std::vector<bool> Visited( topIn.Natom(), true );
    for (int at = res.FirstAtom(); at != res.LastAtom(); at++)
      if (at != *ringat)
        Visited[at] = false;
    Iarray ring_atoms( topIn.Res(rnum).NumAtoms(), -1 );
    // Will be set true if complete ring can be found
    bool ring_complete = false;

    // Since we have already established that *ringat is an oxygen bonded
    // to two carbons, just start at the first carbon to see if we can
    // get to the second carbon.
    int c_beg, c_end;
    if (catoms->first < catoms->second) {
      c_beg = catoms->first;
      c_end = catoms->second;
    } else {
      c_beg = catoms->second;
      c_end = catoms->first;
    }
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
      mprintf("DEBUG:\t(%s bonds to C= %i, %s bonds to C = %i)\n", // DEBUG
              topIn.ResNameNumAtomNameNum(c_beg).c_str(), c_beg_bonds_to_C,
              topIn.ResNameNumAtomNameNum(c_end).c_str(), c_end_bonds_to_C);
    if (c_beg_bonds_to_C <= c_end_bonds_to_C) {
      anomeric_atom = c_beg;
      ring_end_atom = c_end;
//      ring_direction = 1;
    } else {
      anomeric_atom = c_end;
      ring_end_atom = c_beg;
//      ring_direction = -1;
    }
//    mprintf("DEBUG: Potential Ring direction= %i\n", ring_direction);
    catoms->first = anomeric_atom;
    catoms->second = ring_end_atom;

    FollowBonds(anomeric_atom, topIn, 0, ring_atoms,
                ring_end_atom, Visited, ring_complete);
    if (debug_ > 0)
      mprintf("DEBUG: Potential ring start atom %s, Ring complete = %i",
              topIn.ResNameNumAtomNameNum(*ringat).c_str(), (int)ring_complete);
    // Create empty array for ring atoms
    Ring_Atoms.push_back(Iarray());
    if (ring_complete) {
      // Able to complete the cycle.
      if (ring_atom_idx == -1)
        ring_atom_idx = (int)(ringat - potentialRingStartAtoms.begin());
      else {
        mprinterr("Error: Multiple potential ring atoms: %s and %s\n",
                  topIn.ResNameNumAtomNameNum(potentialRingStartAtoms[ring_atom_idx]).c_str(),
                  topIn.ResNameNumAtomNameNum(*ringat).c_str());
        return Sugar(Sugar::MULTIPLE_O, res.FirstAtom());
      }
      // Place the ring atoms into an array without the terminating -1
      Iarray& RA = Ring_Atoms.back();
      if (debug_ > 0) mprintf(" :"); // DEBUG
      for (Iarray::const_iterator it = ring_atoms.begin(); it != ring_atoms.end(); ++it)
      {
        if (debug_ > 0) mprintf(" %i", *it + 1);
        if (*it == -1) break;
        RA.push_back( *it );
      }
    }
    if (debug_ > 0) mprintf("\n"); // DEBUG
  } // END loop over potential ring atoms
  if (ring_atom_idx == -1) {
    mprinterr("Error: Sugar ring oxygen could not be identified.\n");
    return Sugar(Sugar::MISSING_O, res.FirstAtom());
  }

  ring_oxygen_atom = potentialRingStartAtoms[ring_atom_idx];
  anomeric_atom = ringAtomCarbons[ring_atom_idx].first;
  ring_end_atom = ringAtomCarbons[ring_atom_idx].second;
  Iarray const& RA = Ring_Atoms[ring_atom_idx];

  // Find anomeric reference atom. Start at ring end and work down to anomeric atom
  for (Iarray::const_iterator arat = RA.end() - 1; arat != RA.begin(); --arat)
    if (atomIsChiral[*arat - topIn.Res(rnum).FirstAtom()]) {
      ano_ref_atom = *arat;
      break;
    }

  // Get complete chain starting from the anomeric carbon
  carbon_chain = RA;
  if (FindRemainingChainCarbons(carbon_chain, ring_end_atom, topIn, rnum, RA)) {
    mprinterr("Error: Could not find remaining chain carbons.\n");
    return Sugar(Sugar::MISSING_CHAIN, res.FirstAtom());
  }
  if (debug_ > 0) {
    mprintf("DEBUG: Complete carbon chain (from anomeric carbon):\n");
    for (Iarray::const_iterator it = carbon_chain.begin(); it != carbon_chain.end(); ++it)
      mprintf("\t\t%s\n", topIn.ResNameNumAtomNameNum(*it).c_str());
  }
  // See if there is chain prior to anomeric carbon
  Iarray previous_chain;
  if (FindRemainingChainCarbons(previous_chain, anomeric_atom, topIn, rnum, RA)) {
    mprinterr("Error: Could not find previous chain carbons.\n");
    return Sugar(Sugar::MISSING_CHAIN, res.FirstAtom());
  }
  if (!previous_chain.empty()) {
    if (debug_ > 0) {
      mprintf("DEBUG: Previous carbon chain (from anomeric carbon):\n");
      for (Iarray::const_iterator it = previous_chain.begin(); it != previous_chain.end(); ++it)
        mprintf("\t\t%s\n", topIn.ResNameNumAtomNameNum(*it).c_str());
    }
    for (Iarray::const_iterator it = carbon_chain.begin(); it != carbon_chain.end(); ++it)
      previous_chain.push_back( *it );
    carbon_chain = previous_chain;
    if (debug_ > 0) {
      mprintf("DEBUG: Complete carbon chain:\n");
      for (Iarray::const_iterator it = carbon_chain.begin(); it != carbon_chain.end(); ++it)
        mprintf("\t\t%s\n", topIn.ResNameNumAtomNameNum(*it).c_str());
    }
  } 
  // Get the index of the highest stereocenter
  for (Iarray::const_iterator it = carbon_chain.begin(); it != carbon_chain.end(); ++it)
  {
    if (atomIsChiral[*it - topIn.Res(rnum).FirstAtom()])
      highest_stereocenter = *it;
  }
  if (debug_ > 0)
    mprintf("DEBUG: Index of highest stereocenter: %s\n",
            topIn.ResNameNumAtomNameNum(highest_stereocenter).c_str());

  if (ano_ref_atom == -1) {
    mprinterr("Error: Anomeric reference atom could not be identified.\n");
    return Sugar(Sugar::MISSING_ANO_REF, ring_oxygen_atom, anomeric_atom, RA, carbon_chain);
  }
  if (highest_stereocenter == -1) {
    mprinterr("Error: Highest stereocenter atom could not be identified.\n");
    return Sugar(Sugar::MISSING_CONFIG, ring_oxygen_atom, anomeric_atom, RA, carbon_chain);
  }

//  if (!ring_complete || RA.empty() || ring_oxygen_atom == -1) {
//    mprinterr("Error: Sugar ring atoms could not be identified.\n");
//    stat = ID_ERR;
//    return Sugar(res.FirstAtom());
//  }
  if (debug_ > 0)
    mprintf("\t  Ring oxygen         : %s\n", topIn.ResNameNumAtomNameNum(ring_oxygen_atom).c_str());
  return Sugar(ring_oxygen_atom, anomeric_atom, ano_ref_atom, highest_stereocenter,
               RA, carbon_chain, residue_missing_atoms);
}

/** Change PDB atom names in residue to glycam ones. */
int Exec_PrepareForLeap::ChangePdbAtomNamesToGlycam(std::string const& resCode, Residue const& res,
                                                    Topology& topIn, FormTypeEnum form)
const
{
  // Get the appropriate map
  ResIdxMapType::const_iterator resIdxPair = glycam_res_idx_map_.find( resCode );
  if (resIdxPair == glycam_res_idx_map_.end()) {
    // No map needed for this residue
    //mprintf("DEBUG: No atom map for residue '%s'.\n", resCode.c_str());
    return 0;
  }
  NameMapType const& currentMap = pdb_glycam_name_maps_[resIdxPair->second];
  NameMapType const* currentMapAB;
  if (form == ALPHA)
    currentMapAB = &(pdb_glycam_name_maps_A_[resIdxPair->second]);
  else
    currentMapAB = &(pdb_glycam_name_maps_B_[resIdxPair->second]);
  // Change PDB names to Glycam ones
  for (int at = res.FirstAtom(); at != res.LastAtom(); at++)
  {
    NameMapType::const_iterator namePair = currentMapAB->find( topIn[at].Name() );
    if (namePair != currentMapAB->end())
      ChangeAtomName( topIn.SetAtom(at), namePair->second );
    else {
      namePair = currentMap.find( topIn[at].Name() );
      if (namePair != currentMap.end())
        ChangeAtomName( topIn.SetAtom(at), namePair->second );
    }
  }
  return 0;
}

/** Determine if anomeric carbon of furanose is up or down. */
int Exec_PrepareForLeap::DetermineUpOrDown(SugarToken& stoken,
                                         Sugar const& sugar,
                                         Topology const& topIn, Frame const& frameIn)
const
{
  int cdebug;
  if (debug_ > 1)
    cdebug = 1;
  else
    cdebug = 0;
  Cpptraj::Chirality::ChiralType ctypeR = Cpptraj::Chirality::
                                          DetermineChirality(sugar.HighestStereocenter(),
                                                             topIn, frameIn, cdebug);
  if (ctypeR == Cpptraj::Chirality::ERR) {
    mprinterr("Error: Could not determine configuration for furanose.\n"); // TODO warn?
    return 1;
  }
  if (ctypeR == Cpptraj::Chirality::IS_R)
    stoken.SetChirality(IS_D);
  else
    stoken.SetChirality(IS_L);

  Cpptraj::Chirality::ChiralType ctypeA = Cpptraj::Chirality::
                                          DetermineChirality(sugar.AnomericAtom(),
                                                             topIn, frameIn, cdebug);
  if (ctypeA == Cpptraj::Chirality::ERR) {
    mprinterr("Error: Could not determine chirality around anomeric atom for furanose.\n"); // TODO warn?
    return 1;
  }

  if (ctypeR == ctypeA) {
    // Up, beta
    stoken.SetForm(BETA);
  } else {
    // Down, alpha
    stoken.SetForm(ALPHA);
  }
  return 0;
}

/** Determine anomeric form of the sugar. */
int Exec_PrepareForLeap::DetermineAnomericForm(SugarToken& stoken,
                                             Sugar& sugarIn,
                                             Topology const& topIn, Frame const& frameIn)
const
{
  Sugar const& sugar = sugarIn;
  // For determining orientation around anomeric carbon need ring
  // oxygen atom and next carbon in the ring.
  double t_an;
//  ChiralRetType ac_chirality = CalcChiralAtomTorsion(t_an, anomeric_atom, topIn, frameIn);
//  mprintf("DEBUG: Based on t_an %s chirality is %s\n",
//          topIn.TruncResNameOnumId(rnum).c_str(), chiralStr[ac_chirality]);
  int ret = CalcAnomericTorsion(t_an, sugar.AnomericAtom(), sugar.RingOxygenAtom(),
                                sugar.ResNum(topIn),
                                sugar.RingAtoms(), topIn, frameIn);
  if (ret < 0) {
    // This means C1 X substituent missing; non-fatal.
    sugarIn.SetStatus( Sugar::MISSING_C1X );
    return 1; // TODO return 0?
  } else if (ret > 0) {
    // Error
    return 1; 
  }
  bool t_an_up = (t_an > 0);

  // For determining orientation around anomeric reference carbon need
  // previous carbon in the chain and either next carbon or ring oxygen.
  double t_ar;
//    ChiralRetType ar_chirality = CalcChiralAtomTorsion(t_ar, ano_ref_atom, topIn, frameIn);
//    mprintf("DEBUG: Based on t_ar %s chirality is %s\n",
//            topIn.TruncResNameOnumId(rnum).c_str(), chiralStr[ar_chirality]);
  if (CalcAnomericRefTorsion(t_ar, sugar.AnomericRefAtom(), sugar.RingOxygenAtom(), sugar.RingEndAtom(),
                             sugar.RingAtoms(), topIn, frameIn))
  {
    return 1; 
  }
  bool t_ar_up = (t_ar > 0);

  // If config. C is not the anomeric reference, need the previous
  // carbon in the chain, next carbon in the chain, and config. C
  // substituent.
  double t_cc;
  if (sugar.AnomericRefAtom() != sugar.HighestStereocenter()) {
    if (CalcConfigCarbonTorsion(t_cc, sugar.HighestStereocenter(),
                                sugar.ChainAtoms(), topIn, frameIn))
      return 1;
  } else
    t_cc = t_ar;
  bool t_cc_up = (t_cc > 0);

  // Determine index of anomeric atom (typically index 0 but not always).
  int aa_idx =  AtomIdxInArray(sugar.ChainAtoms(), sugar.AnomericAtom());
  int aa_pos = (aa_idx % 2);
  // Determine index of the anomeric reference atom in the chain.
  int ar_idx = AtomIdxInArray(sugar.ChainAtoms(), sugar.AnomericRefAtom());
  int cc_idx = AtomIdxInArray(sugar.ChainAtoms(), sugar.HighestStereocenter());

  // Determine form and chirality.
  // May need to adjust definitions based on the positions of the anomeric
  // reference and config. atoms in the sequence, which alternates.
  if ((ar_idx % 2) != aa_pos)
    t_ar_up = !t_ar_up;
  if ((cc_idx % 2) != aa_pos)
    t_cc_up = !t_cc_up;

  if ( debug_ > 0) {
    mprintf("DEBUG: Index of the anomeric reference atom is %i\n", ar_idx);
    mprintf("DEBUG: Index of the config. carbon atom is %i\n", cc_idx);
    mprintf("DEBUG: t_an_up=%i  t_ar_up=%i  t_cc_up=%i\n",
            (int)t_an_up, (int)t_ar_up, (int)t_cc_up);
  }

  // Same side is beta, opposite is alpha.
  if (t_an_up == t_ar_up) {
    stoken.SetForm(BETA); //form = IS_BETA;
    //mprintf("DEBUG: Form is Beta\n");
  } else {
    stoken.SetForm(ALPHA); //form = IS_ALPHA;
    //mprintf("DEBUG: Form is Alpha\n");
  }

  // By the atom ordering used by CalcAnomericRefTorsion and
  // CalcConfigCarbonTorsion, D is a negative (down) torsion.
  if (!t_cc_up)
    stoken.SetChirality(IS_D);
  else
    stoken.SetChirality(IS_L);

  return 0;
}

/// \return Glycam linkage code for given linked atoms
std::string Exec_PrepareForLeap::GlycamLinkageCode(std::set<Link> const& linkages,
                                                   Topology const& topIn)
const
{
  std::string linkcode;

  // Try to create a link string based on link atom element and position.
  // Check for any unknown positions.
  std::string linkstr;
  for (std::set<Link>::const_iterator it = linkages.begin(); it != linkages.end(); ++it) {
    if (it->Position() < 1) {
      mprinterr("Error: Linkage for atom '%s' has undetermined position in sugar.\n",
                topIn.AtomMaskName(it->Idx()).c_str());
      return linkcode;
    }
    // Carbon is terminal
    if (topIn[it->Idx()].Element() == Atom::CARBON)
      linkstr.append("T");
    else
      linkstr.append( std::string(topIn[it->Idx()].ElementName()) +
                      integerToString(it->Position()) );
  }

  if (debug_ > 0)
    mprintf("DEBUG:\t  linkstr= '%s'\n", linkstr.c_str());
  if      (linkstr == "T") linkcode = "0";
  else if (linkstr == "O1") linkcode = "1";
  else if (linkstr == "TO2") linkcode = "2";
  else if (linkstr == "O2")  linkcode = "2"; // Furanose C2-O2-X
  else if (linkstr == "TO3") linkcode = "3";
  else if (linkstr == "TO4") linkcode = "4";
  else if (linkstr == "TO5") linkcode = "5";
  else if (linkstr == "TO6") linkcode = "6";
  else if (linkstr == "TO2O3") linkcode = "Z";
  else if (linkstr == "TO2O4") linkcode = "Y";
  else if (linkstr == "TO2O6") linkcode = "X";
  else if (linkstr == "TO3O4") linkcode = "W";
  else if (linkstr == "TO3O6") linkcode = "V";
  else if (linkstr == "TO4O6") linkcode = "U";
  else if (linkstr == "TO2O3O4") linkcode = "T";
  else if (linkstr == "TO2O3O6") linkcode = "S";
  else if (linkstr == "TO2O4O6") linkcode = "R";
  else if (linkstr == "TO3O4O6") linkcode = "Q";
  else if (linkstr == "TO2O3O4O6") linkcode = "P";
  if (linkcode.empty())
    mprintf("Warning: Could not determine link code for link atoms '%s'.\n", linkstr.c_str());
  return linkcode;
}

/** Determine linkages for the sugar. Non-sugar residues linked to sugars
  * with a recognized linkage will be marked as valid.
  */
std::string Exec_PrepareForLeap::DetermineSugarLinkages(Sugar const& sugar, CharMask const& cmask,
                                                        Topology& topIn, ResStatArray& resStatIn,
                                                        CpptrajFile* outfile,
                                                        std::set<BondType>& sugarBondsToRemove)
const
{
  int rnum = sugar.ResNum(topIn);
  Residue const& res = topIn.SetRes(rnum);

  // If we haven't identified carbon chain atoms we can't determine their
  // position and therefore the link code.
  if (sugar.ChainAtoms().empty()) {
    mprinterr("Error: No chain atoms determined for '%s', cannot determine link positions.\n",
              topIn.TruncResNameOnumId(rnum).c_str());
    resStatIn[rnum] = SUGAR_NO_CHAIN_FOR_LINK;
    return std::string("");
  }

  // Use set to store link atoms so it is ordered by position
  std::set<Link> linkages;

  // Bonds to non sugars to be removed since these will confuse tleap
  BondArray bondsToRemove;

  // Loop over sugar atoms
  for (int at = res.FirstAtom(); at != res.LastAtom(); at++)
  {
    // Ignore the sugar oxygen and any hydrogen atoms
    if (at == sugar.RingOxygenAtom() ||
        topIn[at].Element() == Atom::HYDROGEN)
      continue;
    // This will be the index of the carbon that is atom is bonded to (or of this atom itself).
    int atomChainPosition = -1;
    // Find position in carbon chain
    if (topIn[at].Element() == Atom::CARBON)
      atomChainPosition = AtomIdxInArray(sugar.ChainAtoms(), at);
    if (atomChainPosition == -1) {
      // Check if any bonded atoms are in carbon chain
      for (Atom::bond_iterator bat = topIn[at].bondbegin();
                               bat != topIn[at].bondend(); ++bat)
      {
        if (topIn[*bat].Element() == Atom::CARBON) {
          atomChainPosition = AtomIdxInArray(sugar.ChainAtoms(), *bat);
          if (atomChainPosition != -1) break;
        }
      }
    }

    // Check for bonds to other residues
    for (Atom::bond_iterator bat = topIn[at].bondbegin();
                             bat != topIn[at].bondend(); ++bat)
    {
      if (topIn[*bat].ResNum() != rnum) {
        // This atom is bonded to another residue.
        linkages.insert(Link(at, atomChainPosition+1));
        // Check if the other residue is a sugar or not
        if (!cmask.AtomInCharMask(*bat)) {
          // Atom is bonded to non-sugar residue.
          mprintf("\t  Sugar %s %s (%i) bonded to non-sugar %s %s (%i) at position %i\n",
                  topIn.TruncResNameOnumId(rnum).c_str(), *(topIn[at].Name()), rnum+1,
                  topIn.TruncResNameOnumId(topIn[*bat].ResNum()).c_str(),
                  *(topIn[*bat].Name()), topIn[*bat].ResNum()+1, atomChainPosition+1);
          bondsToRemove.push_back( BondType(at, *bat, -1) );
          // Check if this is a recognized linkage to non-sugar
          Residue& pres = topIn.SetRes( topIn[*bat].ResNum() );
          NameMapType::const_iterator lname = pdb_glycam_linkageRes_map_.find( pres.Name() );
          if (lname != pdb_glycam_linkageRes_map_.end()) {
            if (debug_ > 0)
              mprintf("DEBUG: Link residue name for %s found: %s\n", *(lname->first), *(lname->second));
            ChangeResName( pres, lname->second );
            resStatIn[topIn[*bat].ResNum()] = VALIDATED;
          } else if (pres.Name() == "ROH") {
            if (debug_ > 0)
              mprintf("DEBUG: '%s' is terminal hydroxyl.\n", *(pres.Name()));
            resStatIn[topIn[*bat].ResNum()] = VALIDATED;
          } else if (pres.Name() == "SO3") {
            if (debug_ > 0)
              mprintf("DEBUG: '%s' is a sulfate group.\n", *(pres.Name()));
            resStatIn[topIn[*bat].ResNum()] = VALIDATED;
          } else if (pres.Name() == "MEX") {
            if (debug_ > 0)
              mprintf("DEBUG: '%s' is a methyl group.\n", *(pres.Name()));
            resStatIn[topIn[*bat].ResNum()] = VALIDATED;
          } else if (pres.Name() == "OME") {
            if (debug_ > 0)
              mprintf("DEBUG: '%s' is an O-methyl group.\n", *(pres.Name()));
            resStatIn[topIn[*bat].ResNum()] = VALIDATED;
          } else {
            mprintf("Warning: Unrecognized link residue %s, not modifying name.\n", *pres.Name());
            resStatIn[topIn[*bat].ResNum()] = SUGAR_UNRECOGNIZED_LINK_RES;
          }
        } else {
          // Atom is bonded to sugar residue
          mprintf("\t  Sugar %s %s bonded to sugar %s %s at position %i\n",
                  topIn.TruncResNameOnumId(rnum).c_str(), *(topIn[at].Name()),
                  topIn.TruncResNameOnumId(topIn[*bat].ResNum()).c_str(),
                  *(topIn[*bat].Name()), atomChainPosition+1);
          // Also remove inter-sugar bonds since leap cant handle branching
          if (at < *bat)
            sugarBondsToRemove.insert( BondType(at, *bat, -1) );
          else
            sugarBondsToRemove.insert( BondType(*bat, at, -1) );
        }
      }
    } // END loop over bonded atoms

  } // END loop over sugar atoms

  // Determine linkage
  //if (debug_ > 0) {
    mprintf("\t  Link atoms:");
    for (std::set<Link>::const_iterator it = linkages.begin();
                                        it != linkages.end(); ++it)
      mprintf(" %s(%i)", *(topIn[it->Idx()].Name()), it->Position());
    mprintf("\n");
  //}
  std::string linkcode;
  if (linkages.empty()) {
    mprintf("\t  No linkages (may be missing atoms).\n");
    resStatIn[rnum] = SUGAR_NO_LINKAGE;
    return linkcode;
  } else {
    linkcode = GlycamLinkageCode(linkages, topIn);
    if (debug_ > 0) mprintf("\t  Linkage code: %s\n", linkcode.c_str());
    if (linkcode.empty()) {
      mprinterr("Error: Unrecognized sugar linkage.\n");
      resStatIn[rnum] = SUGAR_UNRECOGNIZED_LINKAGE;
      return linkcode;
    }
  }
  // Remove bonds to other residues 
  for (BondArray::const_iterator bnd = bondsToRemove.begin();
                                 bnd != bondsToRemove.end(); ++bnd)
  {
    LeapBond(bnd->A1(), bnd->A2(), topIn, outfile);
    topIn.RemoveBond(bnd->A1(), bnd->A2());
  }

  return linkcode;
}

/** Attempt to identify sugar residue, form, and linkages. */
int Exec_PrepareForLeap::IdentifySugar(Sugar& sugarIn, Topology& topIn,
                                       Frame const& frameIn, CharMask const& cmask,
                                       CpptrajFile* outfile, std::set<BondType>& sugarBondsToRemove)
{
  Sugar const& sugar = sugarIn;
  const std::string sugarName = topIn.TruncResNameOnumId(sugar.ResNum(topIn));
  const std::string resName = topIn.Res(sugar.ResNum(topIn)).Name().Truncated();
  //if (sugar.NotSet()) {
  //  mprintf("Warning: Sugar %s is not set up. Skipping sugar identification.\n", sugarName.c_str());
  //  return 0; // TODO return 1?
  //}

  int rnum = sugar.ResNum(topIn);
  Residue& res = topIn.SetRes(rnum);
  // Try to ID the base sugar type from the input name.
  //char resChar = ' ';

  MapType::const_iterator pdb_glycam;
  if (!hasGlycam_) {
    pdb_glycam = pdb_to_glycam_.find( res.Name() );
    //resChar = pdb_glycam->second;
  } else {
    mprintf("\tGlycam res: '%s'\n", resName.c_str());
    // Assume residue name is already glycam. Find it. Get sugar from 2nd char.
    std::string rChar(1, resName[1]);
    // Determine chirality from upper/lowercase
    ChirTypeEnum cType = UNKNOWN_CHIR;
    if (islower(rChar.front()))
      cType = IS_L;
    else
      cType = IS_D;
    // Now force uppercase
    rChar.assign( 1, toupper(rChar.front()) );
    std::string fChar(1, resName[2]);
    mprintf("DEBUG: Searching for glycam char '%s' '%s'\n", rChar.c_str(), fChar.c_str());
    // Get form from 3rd char
    FormTypeEnum fType = UNKNOWN_FORM;
    if (fChar == "U" || fChar == "B")
      fType = BETA;
    else if (fChar == "D" || fChar == "A")
      fType = ALPHA;
    for (pdb_glycam = pdb_to_glycam_.begin(); pdb_glycam != pdb_to_glycam_.end(); ++pdb_glycam)
      if (pdb_glycam->second.GlycamCode() == rChar &&
          pdb_glycam->second.Form() == fType &&
          pdb_glycam->second.Chirality() == cType)
        break;
  }

  if ( pdb_glycam == pdb_to_glycam_.end() ) { 
    mprinterr("Error: Could not identify sugar from residue name '%s'\n", *res.Name());
    return 1;
  }

  mprintf("\tSugar %s %s\n", sugarName.c_str(), pdb_glycam->second.InfoStr().c_str());

  SugarToken sugarInfo;
  int detect_err = 1;

  if (!sugar.NotSet()) {
    sugarInfo = SugarToken(sugar.RingType());

    // Determine alpha or beta and D or L
    if (sugarInfo.RingType() == FURANOSE)
      detect_err = DetermineUpOrDown(sugarInfo, sugar, topIn, frameIn);
    else if (sugarInfo.RingType() == PYRANOSE)
      detect_err = DetermineAnomericForm(sugarInfo, sugarIn, topIn, frameIn);
    else
      detect_err = 1;
  }

  if (detect_err != 0) {
    // This means there was an issue determining ring type, form,
    // chirality, or both.
    // Try to fill in info based on the residue name (pdb_glycam->second).
    if (sugarInfo.Form() == UNKNOWN_FORM) {
      mprintf("Warning: Could not determine anomer type from coordinates.\n");
      if (pdb_glycam->second.Form() == UNKNOWN_FORM)
        return 0;
      mprintf("Warning: Setting anomer type based on residue name.\n");
      sugarInfo.SetForm( pdb_glycam->second.Form() );
    }
    if (sugarInfo.Chirality() == UNKNOWN_CHIR) {
      mprintf("Warning: Could not determine configuration from coordinates.\n");
      if (pdb_glycam->second.Chirality() == UNKNOWN_CHIR)
        return 0;
      mprintf("Warning: Setting configuration based on residue name.\n");
      sugarInfo.SetChirality( pdb_glycam->second.Chirality() );
    }
    if (sugarInfo.RingType() == UNKNOWN_RING) {
      mprintf("Warning: Could not determine ring type from coordinates.\n");
      if (pdb_glycam->second.RingType() == UNKNOWN_RING)
        return 0;
      mprintf("Warning: Setting ring type based on residue name.\n");
      sugarInfo.SetRingType( pdb_glycam->second.RingType() );
    }
  }
  // Warn about form/chirality mismatches. If a mismatch occurs and the sugar
  // did not have all atoms present, defer to the name.
  if (pdb_glycam->second.Form() != UNKNOWN_FORM &&
      pdb_glycam->second.Form() != sugarInfo.Form())
  {
    mprintf("Warning: '%s' detected anomer type is %s but anomer type based on name is %s.\n",
            sugarName.c_str(), formstr_[sugarInfo.Form()],
            formstr_[pdb_glycam->second.Form()]);
    if (sugarInfo.Form() != UNKNOWN_FORM)
      resStat_[rnum] = SUGAR_NAME_MISMATCH;
    if (useSugarName_) {
      mprintf("\tSetting anomer type based on residue name.\n");
      sugarInfo.SetForm( pdb_glycam->second.Form() );
    } else if (sugar.IsMissingAtoms()) {
      mprintf("Warning: Residue was missing atoms; setting anomer type based on residue name.\n");
      sugarInfo.SetForm( pdb_glycam->second.Form() );
    }
  }
  if (pdb_glycam->second.Chirality() != UNKNOWN_CHIR &&
      pdb_glycam->second.Chirality() != sugarInfo.Chirality())
  {
    mprintf("Warning: '%s' detected configuration is %s but configuration based on name is %s.\n",
            sugarName.c_str(), chirstr_[sugarInfo.Chirality()],
            chirstr_[pdb_glycam->second.Chirality()]);
    if (sugarInfo.Chirality() != UNKNOWN_CHIR)
      resStat_[rnum] = SUGAR_NAME_MISMATCH;
    if (useSugarName_) {
      mprintf("\tSetting configuration based on residue name.\n");
      sugarInfo.SetChirality( pdb_glycam->second.Chirality() );
    } else if (sugar.IsMissingAtoms()) {
      mprintf("Warning: Residue was missing atoms; setting configuration based on residue name.\n");
      sugarInfo.SetChirality( pdb_glycam->second.Chirality() );
    }
  }
  if (pdb_glycam->second.RingType() != UNKNOWN_RING &&
      pdb_glycam->second.RingType() != sugarInfo.RingType())
  {
    mprintf("Warning: '%s' detected ring type is %s but ring type based on name is %s.\n",
            sugarName.c_str(), ringstr_[sugarInfo.RingType()],
            ringstr_[pdb_glycam->second.RingType()]);
    if (sugarInfo.RingType() != UNKNOWN_RING)
      resStat_[rnum] = SUGAR_NAME_MISMATCH;
    if (useSugarName_) {
      mprintf("Setting ring type bsaed on residue name.\n");
      sugarInfo.SetRingType( pdb_glycam->second.RingType() );
    } else if (sugar.IsMissingAtoms()) {
      mprintf("Warning: Residue was missing atoms; setting ring type based on residue name.\n");
      sugarInfo.SetRingType( pdb_glycam->second.RingType() );
    }
  }
    
  // Get glycam form string
  std::string formStr;
  if (sugarInfo.RingType() == PYRANOSE) {
    if (sugarInfo.Form() == ALPHA)
      formStr = "A";
    else
      formStr = "B";
  } else if (sugarInfo.RingType() == FURANOSE) {
    if (sugarInfo.Form() == ALPHA)
      formStr = "D";
    else
      formStr = "U";
  }
  // Sanity check
  if (formStr.empty()) {
    mprinterr("Internal Error: Could not set anomer type string.\n");
    return 1;
  }

  mprintf("\t  %s detected anomer type is %s(%s)-%s-%s\n", sugarName.c_str(),
         formstr_[sugarInfo.Form()], formStr.c_str(), chirstr_[sugarInfo.Chirality()],
         ringstr_[sugarInfo.RingType()]);

  // Identify linkages to other residues.
  std::string linkcode = DetermineSugarLinkages(sugar, cmask, topIn, resStat_,
                                                outfile, sugarBondsToRemove);
  if (linkcode.empty()) {
    //resStat_[rnum] = SUGAR_UNRECOGNIZED_LINKAGE;
    mprintf("Warning: Determination of sugar linkages failed.\n");
    return 0;
  }

  // Change PDB names to Glycam ones
  std::string resCode = pdb_glycam->second.GlycamCode();
  if (!hasGlycam_) {
    if (ChangePdbAtomNamesToGlycam(resCode, res, topIn, sugarInfo.Form()))
    {
      mprinterr("Error: Changing PDB atom names to Glycam failed.\n");
      return 1;
    }
  }

  // Modify residue char to indicate L form if necessary.
  // We do this here and not above so as not to mess with the
  // linkage determination.
  if (sugarInfo.Chirality() == IS_L) {
    for (std::string::iterator it = resCode.begin(); it != resCode.end(); ++it)
      *it = tolower( *it );
  }
  // Set new residue name
  NameType newResName( linkcode + resCode + formStr );
  if (!hasGlycam_) {
    mprintf("\t  Changing %s to Glycam resname: %s\n",
      topIn.TruncResNameOnumId(rnum).c_str(), *newResName);
    ChangeResName(res, newResName);
  } else {
    if (newResName.Truncated() != resName)
      mprintf("Warning: Detected glycam name '%s' differs from original name '%s'\n",
              *newResName, resName.c_str());
  }
  if (resStat_[rnum] == UNKNOWN)
    resStat_[rnum] = VALIDATED;
  return 0;
}

/** Attempt to find any missing linkages to the anomeric carbon in sugar. */
int Exec_PrepareForLeap::FindSugarC1Linkages(int rnum1, int c_beg,
                                             Topology& topIn, Frame const& frameIn)
const
{
  Residue const& res1 = topIn.SetRes(rnum1);
  // If the anomeric atom is already bonded to another residue, skip this.
  for (Atom::bond_iterator bat = topIn[c_beg].bondbegin();
                           bat != topIn[c_beg].bondend(); ++bat)
  {
    if (topIn[*bat].ResNum() != rnum1) {
      if (debug_ > 0)
        mprintf("\tSugar %s anomeric carbon is already bonded to another residue, skipping.\n",
                topIn.TruncResNameOnumId(rnum1).c_str());
      return 0;
    }
  }

  // residue first atom to residue first atom cutoff^2
  const double rescut2 = 64.0;
  // bond cutoff offset
  const double offset = 0.2;
  // index of atom to be bonded to c_beg
  int closest_at = -1;
  // distance^2 of atom to be bonded to c_beg
  double closest_d2 = -1.0;

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
            mprintf("DEBUG: Residue %s to %s = %f\n",
                    topIn.TruncResNameOnumId(rnum1).c_str(), topIn.TruncResNameOnumId(rnum2).c_str(),
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
                if (debug_ > 1)
                  mprintf("DEBUG: Atom %s to %s = %f\n",
                          topIn.AtomMaskName(c_beg).c_str(), topIn.AtomMaskName(at2).c_str(), sqrt(D2));
                if (closest_at == -1) {
                  closest_at = at2;
                  closest_d2 = D2;
                } else if (D2 < closest_d2) {
                  mprintf("\t  Atom %s (%f Ang.) is closer than %s (%f Ang.).\n",
                          topIn.ResNameNumAtomNameNum(at2).c_str(), sqrt(D2),
                          topIn.ResNameNumAtomNameNum(closest_at).c_str(), sqrt(closest_d2));
                  closest_at = at2;
                  closest_d2 = D2;
                }
                //mprintf("\t  Adding bond between %s and %s\n",
                //        topIn.ResNameNumAtomNameNum(c_beg).c_str(),
                //        topIn.ResNameNumAtomNameNum(at2).c_str());
                //topIn.AddBond(c_beg, at2);
              }
            }
          } // END loop over res2 atoms
        } // END res1-res2 distance cutoff
      } // END res2 is not solvent
    } // END res1 != res2
  } // END res2 loop over other residues
  if (closest_at != -1) {
    mprintf("\t  Adding bond between %s and %s\n",
            topIn.ResNameNumAtomNameNum(c_beg).c_str(),
            topIn.ResNameNumAtomNameNum(closest_at).c_str());
    topIn.AddBond(c_beg, closest_at);
  }
  return 0;
}

/** \return Type of group represented by the atom atIdx. */
Exec_PrepareForLeap::FunctionalGroupType
  Exec_PrepareForLeap::IdFunctionalGroup_Silent(Iarray& selected, int rnum, int atIdx, int linkAtIdx, Topology const& topIn)
const
{
  selected.clear();
  // ----- Sulfur ----------------------
  if (topIn[atIdx].Element() == Atom::SULFUR &&
      topIn[atIdx].ResNum() == rnum &&
      topIn[atIdx].Nbonds() == 4) {
    //so3_idx = atIdx;
    selected.push_back(atIdx);
    // All 4 bonds must be to oxygen
    int bonds_to_o = 0;
    for (Atom::bond_iterator bat = topIn[atIdx].bondbegin();
                             bat != topIn[atIdx].bondend(); ++bat)
    {
      if (topIn[*bat].Element() == Atom::OXYGEN &&
          topIn[*bat].ResNum() == rnum) {
        bonds_to_o++;
        if (*bat != linkAtIdx)
          selected.push_back( *bat );
      }
    }
    if (bonds_to_o == 4) {
      return G_SO3;
    }
  // ----- Carbon ----------------------
  } else if (topIn[atIdx].Element() == Atom::CARBON &&
             topIn[atIdx].ResNum() == rnum) {
    //so3_idx = atIdx;
    selected.push_back(atIdx);
    if (topIn[atIdx].Nbonds() == 4) {
      // If 4 bonds, 3 must be to hydrogen for CH3
      int bonds_to_h = 0;
      for (Atom::bond_iterator bat = topIn[atIdx].bondbegin();
                               bat != topIn[atIdx].bondend(); ++bat)
      {
        if (topIn[*bat].Element() == Atom::HYDROGEN) {
          bonds_to_h++;
          selected.push_back(*bat);
        }
      }
      if (bonds_to_h == 3) {
        return G_CH3;
      }
    } else if (topIn[atIdx].Nbonds() == 1)
      return G_CH3;
    } else if (topIn[atIdx].Nbonds() == 3) {
      //             O
      //             ||
      // Check for - C - CH3
      bool has_o = false;
      bool has_ch3 = false;
      Iarray ch3_atoms;
      for (Atom::bond_iterator bat = topIn[atIdx].bondbegin();
                               bat != topIn[atIdx].bondend(); ++bat)
      {
        if (*bat != linkAtIdx) {
          if (topIn[*bat].Element() == Atom::OXYGEN) {
            has_o = true;
            selected.push_back(*bat);
          } else if (topIn[*bat].Element() == Atom::CARBON) {
            FunctionalGroupType fgt = IdFunctionalGroup_Silent(ch3_atoms, rnum, *bat, atIdx, topIn);
            has_ch3 = (fgt == G_CH3);
          } else {
            // If bonded to anything else, not Acetyl
            has_o = false;
            has_ch3 = false;
            break;
          }
        }
      }
      if (has_o && has_ch3) {
        for (Iarray::const_iterator cit = ch3_atoms.begin(); cit != ch3_atoms.end(); ++cit)
          selected.push_back(*cit);
        return G_ACX;
      }
  // ----- Oxygen ----------------------
  } else if (topIn[atIdx].Element() == Atom::OXYGEN &&
             topIn[atIdx].ResNum() == rnum) {
    selected.push_back( atIdx );
    if (topIn[atIdx].Nbonds() == 1)
      // If only 1 bond, -OH
      return G_OH;
    else if (topIn[atIdx].Nbonds() == 2) {
      int bonded_atom;
      if (topIn[atIdx].Bond(0) == linkAtIdx)
        bonded_atom = topIn[atIdx].Bond(1);
      else
        bonded_atom = topIn[atIdx].Bond(0);
      if (topIn[bonded_atom].Element() == Atom::HYDROGEN) {
        // -OH
        selected.push_back(bonded_atom);
        return G_OH;
      } else if (topIn[bonded_atom].Element() == Atom::CARBON) {
        // Might be -OCH3
        int c_idx = bonded_atom;
        selected.push_back(c_idx);
        // Check for only 1 bond or 3 bonds to hydrogen
        int bonds_to_h = 0;
        for (Atom::bond_iterator bat = topIn[c_idx].bondbegin();
                                 bat != topIn[c_idx].bondend(); ++bat)
        {
          if (topIn[*bat].Element() == Atom::HYDROGEN) {
            bonds_to_h++;
            selected.push_back(*bat);
          }
        }
        if (topIn[c_idx].Nbonds() == 1 || bonds_to_h == 3)
          return G_OME;
      }
    }

  }
  selected.clear();      
  return UNRECOGNIZED_GROUP;
}

/** Corresponds to FunctionalGroupType */
const char* Exec_PrepareForLeap::FunctionalGroupStr_[] = {
  "SO3", "CH3", "Acetyl", "OH", "OCH3", "Unrecognized"
};

/** Identify functional group and print to stdout. */
Exec_PrepareForLeap::FunctionalGroupType
  Exec_PrepareForLeap::IdFunctionalGroup(Iarray& selected, int rnum, int atIdx, int linkAtIdx, Topology const& topIn)
const
{
  FunctionalGroupType groupType = IdFunctionalGroup_Silent(selected, rnum, atIdx, linkAtIdx, topIn);
  if (groupType != UNRECOGNIZED_GROUP) {
    mprintf("\tFound %s group centered on atom '%s %s' bonded to '%s'\n",
            FunctionalGroupStr_[groupType],
            topIn.TruncResNameOnumId(rnum).c_str(),
            *(topIn[atIdx].Name()), *(topIn[linkAtIdx].Name()));
  }

  return groupType;
}

/** Check for functional groups that need to be separate residues. */
int Exec_PrepareForLeap::CheckForFunctionalGroups(Sugar& sugar,
                                               Topology& topIn, Frame& frameIn)
const
{
  int rnum = sugar.ResNum(topIn);
  int original_at0 = topIn.Res(rnum).FirstAtom();
  int original_at1 = topIn.Res(rnum).LastAtom();
  std::string sugarName = topIn.TruncResNameOnumId(rnum);
  if (debug_ > 0) mprintf("DEBUG: Functional group check: %s\n", sugarName.c_str());

  bool atomsRemain = true;
  while (atomsRemain) {
    // Even if the residue is split, rnum will always refer to the original
    // sugar since the split SO3 will come AFTER this residue.
    // Find an oxygen that is both bound to a chain carbon and an SO3 group
    // in this residue.
    int o_idx = -1;
    int so3_idx = -1;
    Iarray selected;
    FunctionalGroupType groupType = UNRECOGNIZED_GROUP;
    Iarray::const_iterator cat = sugar.ChainAtoms().begin();
    for (; cat != sugar.ChainAtoms().end(); ++cat)
    {
      //mprintf("\t%s\n", *(topIn[*cat].Name())); // DEBUG
      for (Atom::bond_iterator oat = topIn[*cat].bondbegin();
                               oat != topIn[*cat].bondend(); ++oat)
      {
        o_idx = -1;
        if (topIn[*oat].Element() == Atom::OXYGEN &&
            topIn[*oat].ResNum() == rnum &&
            topIn[*oat].Nbonds() > 1) {
          o_idx = *oat;
          // Is this oxygen bound to a recognized group?
          for (Atom::bond_iterator sat = topIn[*oat].bondbegin();
                                   sat != topIn[*oat].bondend(); ++sat)
          {
            groupType = IdFunctionalGroup(selected, rnum, *sat, o_idx, topIn);
            if (groupType != UNRECOGNIZED_GROUP) {
              so3_idx = *sat;
              break;
            }
          } // END loop over bonds to oxygen
          if (so3_idx != -1) break;
        } // END atom is oxygen
      } // END loop over bonds to carbon
      if (so3_idx != -1) break;
    } // END loop over chain atoms
    if (cat == sugar.ChainAtoms().end())
      atomsRemain = false;
    else {
      // sanity check
      if (so3_idx == -1 || o_idx == -1 || selected.empty()) {
        mprinterr("Internal Error: Functional group index is negative.\n");
        return 1;
      }
      
      // Change the atom names
      std::string newResName;
      if (groupType == G_SO3) {
        ChangeAtomName(topIn.SetAtom(selected[0]), "S1");
        ChangeAtomName(topIn.SetAtom(selected[1]), "O1");
        ChangeAtomName(topIn.SetAtom(selected[2]), "O2");
        ChangeAtomName(topIn.SetAtom(selected[3]), "O3");
        newResName = "SO3";
      } else if (groupType == G_CH3) {
        ChangeAtomName(topIn.SetAtom(selected[0]), "CH3");
        if (selected.size() > 1) {
          ChangeAtomName(topIn.SetAtom(selected[1]), "H1");
          ChangeAtomName(topIn.SetAtom(selected[2]), "H2");
          ChangeAtomName(topIn.SetAtom(selected[3]), "H3");
        }
        newResName = "MEX";
      } else if (groupType == G_ACX) {
        ChangeAtomName(topIn.SetAtom(selected[0]), "C1A");
        ChangeAtomName(topIn.SetAtom(selected[1]), "O1A");
        ChangeAtomName(topIn.SetAtom(selected[2]), "C2A");
        if (selected.size() > 3) {
          ChangeAtomName(topIn.SetAtom(selected[3]), "H1A");
          ChangeAtomName(topIn.SetAtom(selected[4]), "H2A");
          ChangeAtomName(topIn.SetAtom(selected[5]), "H3A");
        }
      } else {
        mprinterr("Internal Error: Unhandled group in CheckForFunctionalGroups()\n");
        return 1;
      }
      // Create array with group selected
      AtomMask FxnGrpMask(selected, topIn.Natom());
      // Split the sulfate into a new residue named newResName for Glycam.
      // This may involve reordering atoms within the residue, but not
      // any other atoms, so we should not have to update other sugars.
      //mprintf("DEBUG: Before split: %s\n", topIn.AtomMaskName(sugar.RingOxygenAtom()).c_str());
      Iarray atomMap;
      if (topIn.SplitResidue(FxnGrpMask, newResName, atomMap)) {
        mprinterr("Error: Could not split functional group from residue '%s'.\n",
                  sugarName.c_str());
        return 1;
      }
      // Set the split residue as terminal
      topIn.SetRes(rnum+1).SetTerminal(true);
      // Reorder the frame to match
      Frame oldFrame = frameIn;
      frameIn.SetCoordinatesByMap( oldFrame, atomMap );
      // Remap the sugar indices
      //mprintf("DEBUG: Before remap: %s\n", topIn.AtomMaskName(sugar.RingOxygenAtom()).c_str());
      sugar.RemapIndices( atomMap, original_at0, original_at1 );
      //mprintf("DEBUG: After remap: %s\n", topIn.AtomMaskName(sugar.RingOxygenAtom()).c_str());
    }
    //atomsRemain = false; // DEBUG
  } // END while atoms remain

  return 0;
}

/** See if the sugar anomeric carbon is actually terminal and needs
  * to be a separate ROH residue.
  */
int Exec_PrepareForLeap::CheckIfSugarIsTerminal(Sugar& sugar, Topology& topIn, Frame& frameIn)
const
{
  int rnum = sugar.ResNum(topIn);
  int anomericAtom = sugar.AnomericAtom();
  int ringOxygen = sugar.RingOxygenAtom();
  std::string sugarName = topIn.TruncResNameOnumId(rnum);

  // Is the anomeric carbon bonded to an oxygen that is part of this residue.
  int o1_atom = -1;
  for (Atom::bond_iterator bat = topIn[anomericAtom].bondbegin();
                           bat != topIn[anomericAtom].bondend(); ++bat)
  {
    if (topIn[*bat].ResNum() == rnum &&
        *bat != ringOxygen &&
        topIn[*bat].Element() == Atom::OXYGEN) {
      if (o1_atom != -1) {
        mprintf("Warning: Anomeric atom '%s %s' bonded to more than 1 oxygen.\n",
                sugarName.c_str(), *(topIn[anomericAtom].Name()));
        o1_atom = -1;
        break;
      } else
        o1_atom = *bat;
    }
  }
  if (o1_atom == -1) return 0;
  //mprintf("DEBUG: Terminal check: %s O1 atom: '%s'\n", sugarName.c_str(), topIn.AtomMaskName(o1_atom).c_str());

  Iarray selected;
  FunctionalGroupType groupType = IdFunctionalGroup(selected, rnum, o1_atom, anomericAtom, topIn);

  std::string newResName;
  if (groupType == G_OH) {
    newResName = "ROH";//terminalHydroxylName_;
    // Change atom names
    ChangeAtomName(topIn.SetAtom(selected[0]), "O1");
    if (selected.size() > 1)
      ChangeAtomName(topIn.SetAtom(selected[1]), "HO1");
  } else if (groupType == G_OME) {
    newResName = "OME";
    // Change atom names
    ChangeAtomName(topIn.SetAtom(selected[0]), "O");
    ChangeAtomName(topIn.SetAtom(selected[1]), "CH3");
    if (selected.size() > 2) {
      ChangeAtomName(topIn.SetAtom(selected[2]), "H1");
      ChangeAtomName(topIn.SetAtom(selected[3]), "H2");
      ChangeAtomName(topIn.SetAtom(selected[4]), "H3");
    }
  } else if (groupType != UNRECOGNIZED_GROUP) {
    mprinterr("Internal Error: Unhandled group in CheckIfSugarIsTerminal\n");
    return 1;
  }
  //mprintf("\t  Will split into %s group.\n", newResName.c_str());

  // If terminal group is unrecognized, this could just be a regular O1 linkage
  if (selected.empty()) {
    if (debug_ > 0)
      mprintf("DEBUG: Sugar '%s' has unrecognized terminal group.\n", sugarName.c_str());
    return 0;
  }
  
  // Create atom mask
  AtomMask ROH(selected, topIn.Natom());

  // Split the hydroxyl into a new residue named ROH for Glycam.
  // This may involve reordering atoms within the residue, but not
  // any other atoms, so we should not have to update SugarIndices.
  int original_at0 = topIn.Res(rnum).FirstAtom();
  int original_at1 = topIn.Res(rnum).LastAtom();
  Iarray atomMap;
  if (topIn.SplitResidue(ROH, newResName, atomMap)) {
    mprinterr("Error: Could not split the residue '%s'.\n", sugarName.c_str());
    return 1;
  }
  // Set the split residue as terminal
  topIn.SetRes(rnum+1).SetTerminal(true);
  // DEBUG
  //for (int at = original_at0; at != original_at1; at++)
  //  mprintf("DEBUG:\t\tAtomMap[%i] = %i\n", at, atomMap[at]);
  // Reorder the frame to match
  Frame oldFrame = frameIn;
  frameIn.SetCoordinatesByMap( oldFrame, atomMap );
  // Remap the sugar indices
  sugar.RemapIndices( atomMap, original_at0, original_at1 );
  return 0;
}


/** Try to fix issues with sugar structure before trying to identify. */
int Exec_PrepareForLeap::FixSugarsStructure(std::vector<Sugar>& sugarResidues,
                                            std::string const& sugarMaskStr,
                                            Topology& topIn, Frame& frameIn,
                                            bool c1bondsearch, bool splitres) 
const
{
  sugarResidues.clear();
  AtomMask sugarMask(sugarMaskStr);
  mprintf("\tLooking for sugars selected by '%s'\n", sugarMask.MaskString());
  if (topIn.SetupIntegerMask( sugarMask )) return 1;
  //sugarMask.MaskInfo();
  mprintf("\tSelected %i sugar atoms.\n", sugarMask.Nselected());
  if (sugarMask.None()) {
    mprintf("Warning: No sugar atoms selected by %s\n", sugarMask.MaskString());
    return 0;
  }
  //CharMask cmask( sugarMask.ConvertToCharMask(), sugarMask.Nselected() );
  // Get sugar residue numbers.
  Iarray sugarResNums = topIn.ResnumsSelectedBy( sugarMask );
  // Try to identify sugar rings.
  for (Iarray::const_iterator rnum = sugarResNums.begin();
                              rnum != sugarResNums.end(); ++rnum)
  {
    Sugar sugar = IdSugarRing(*rnum, topIn);
    if (sugar.Status() != Sugar::SETUP_OK) {
      mprintf("Warning: Problem identifying atoms for sugar '%s'\n",
              topIn.TruncResNameOnumId(*rnum).c_str());
    }
/*
    if (stat == ID_ERR) {
      if (errorsAreFatal_) {
        mprinterr("Error: Problem identifying sugar ring for %s\n",
                  topIn.TruncResNameOnumId(*rnum).c_str());
        return 1;
      } else
        mprintf("Warning: Problem identifying sugar ring for %s\n",
                topIn.TruncResNameOnumId(*rnum).c_str());
    }*/
    //if (!sugar.NotSet()) {
      sugarResidues.push_back( sugar );
      if (debug_ > 0)
        sugarResidues.back().PrintInfo(topIn);
    //}
  }

  if (c1bondsearch) {
    // Loop over sugar indices to see if anomeric C is missing bonds
    for (std::vector<Sugar>::const_iterator sugar = sugarResidues.begin();
                                            sugar != sugarResidues.end(); ++sugar)
    {
      if (sugar->NotSet()) continue;
      int anomericAtom = sugar->AnomericAtom();
      int rnum = sugar->ResNum(topIn);
      if (FindSugarC1Linkages(rnum, anomericAtom, topIn, frameIn)) {
        mprinterr("Error: Search for bonds to anomeric carbon '%s' failed.\n",
                  topIn.AtomMaskName(anomericAtom).c_str());
        return 1;
      }
    }
  }
  //DEBUG
  //for (std::vector<Sugar>::const_iterator sugar = sugarResidues.begin();
  //                                    sugar != sugarResidues.end(); ++sugar)
  //  sugar->PrintInfo(topIn);

  if (splitres) {
    // Loop over sugar indices to see if residues have ROH that must be split off
    for (std::vector<Sugar>::iterator sugar = sugarResidues.begin();
                                      sugar != sugarResidues.end(); ++sugar)
    {
      if (sugar->NotSet()) continue;
      if (CheckIfSugarIsTerminal(*sugar, topIn, frameIn)) {
        mprinterr("Error: Checking if sugar %s has terminal functional groups failed.\n",
                  topIn.TruncResNameOnumId(sugar->ResNum(topIn)).c_str());
        return 1;
      }
    } // End loop over sugar indices

    // Loop over chain indices to see if residues need to be split
    for (std::vector<Sugar>::iterator sugar = sugarResidues.begin();
                                      sugar != sugarResidues.end(); ++sugar)
    {
      if (sugar->NotSet()) continue;
      if (CheckForFunctionalGroups(*sugar, topIn, frameIn)) {
        mprinterr("Error: Checking if sugar %s has functional groups failed.\n",
                 topIn.TruncResNameOnumId( sugar->ResNum(topIn) ).c_str());
        return 1;
      }
    }
    //DEBUG
    //for (std::vector<Sugar>::const_iterator sugar = sugarResidues.begin();
    //                                    sugar != sugarResidues.end(); ++sugar)
    //  sugar->PrintInfo(topIn);
  }


  return 0;
}

/** Prepare sugars for leap. */
int Exec_PrepareForLeap::PrepareSugars(std::string const& sugarmaskstr,
                                       std::vector<Sugar>& Sugars,
                                       Topology& topIn,
                                       Frame const& frameIn, CpptrajFile* outfile)
                                
{
  // Need to set up the mask again since topology may have been modified.
  AtomMask sugarMask;
  if (sugarMask.SetMaskString( sugarmaskstr )) return 1;
  //mprintf("\tPreparing sugars selected by '%s'\n", sugarMask.MaskString());
  if (topIn.SetupIntegerMask( sugarMask )) return 1;
  //sugarMask.MaskInfo();
  mprintf("\t%i sugar atoms selected in %zu residues.\n", sugarMask.Nselected(), Sugars.size());
  if (sugarMask.None())
    mprintf("Warning: No sugar atoms selected by %s\n", sugarMask.MaskString());
  else {
    CharMask cmask( sugarMask.ConvertToCharMask(), sugarMask.Nselected() );
    if (debug_ > 0) {
      for (std::vector<Sugar>::const_iterator sugar = Sugars.begin();
                                              sugar != Sugars.end(); ++sugar)
        sugar->PrintInfo(topIn);
    }
    std::set<BondType> sugarBondsToRemove;
    // For each sugar residue, see if it is bonded to a non-sugar residue.
    // If it is, remove that bond but record it.
    for (unsigned int sidx = 0; sidx != Sugars.size(); sidx++)
    {
      Sugar const& sugar = Sugars[sidx];
      Sugar& sugarIn = Sugars[sidx];
      //if (sugar.NotSet()) {
      //  resStat_[sugar.ResNum(topIn)] = SUGAR_SETUP_FAILED;
      //  continue;
      //}
      // See if we recognize this sugar.
      if (IdentifySugar(sugarIn, topIn, frameIn, cmask, outfile, sugarBondsToRemove))
      {
        if (errorsAreFatal_)
          return 1;
        else
          mprintf("Warning: Preparation of sugar %s failed, skipping.\n",
                  topIn.TruncResNameOnumId( sugar.ResNum(topIn) ).c_str());
      }
    } // END loop over sugar residues
    // Remove bonds between sugars
    for (std::set<BondType>::const_iterator bnd = sugarBondsToRemove.begin();
                                            bnd != sugarBondsToRemove.end(); ++bnd)
    {
      LeapBond(bnd->A1(), bnd->A2(), topIn, outfile);
      topIn.RemoveBond(bnd->A1(), bnd->A2());
    }
    // Bonds to sugars have been removed, so regenerate molecule info
    topIn.DetermineMolecules();
    // Set each sugar as terminal
    for (std::vector<Sugar>::const_iterator sugar = Sugars.begin(); sugar != Sugars.end(); ++sugar)
    {
      int rnum = sugar->ResNum(topIn);
      topIn.SetRes(rnum).SetTerminal(true);
      if (rnum - 1 > -1)
        topIn.SetRes(rnum-1).SetTerminal(true);
    }
  }
  return 0;
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
        resStat_[topIn[at1].ResNum()] = VALIDATED;
      }
    }
    mprintf("\tDetected %i disulfide bonds.\n", nDisulfides);
  }
  return 0;
}

/** \return True if residue name is in pdb_to_glycam_ or pdb_res_names_,
  *               or is solvent.
  */
bool Exec_PrepareForLeap::IsRecognizedPdbRes(NameType const& rname) const {
  MapType::const_iterator glycamIt = pdb_to_glycam_.find( rname );
  if (glycamIt != pdb_to_glycam_.end())
    return true;
  SetType::const_iterator amberIt = pdb_res_names_.find( rname );
  if (amberIt != pdb_res_names_.end())
    return true;
  if (rname == solventResName_)
    return true;
  return false;
}

/** \return Array of residue numbers with unrecognized PDB res names. */
Exec_PrepareForLeap::Iarray Exec_PrepareForLeap::GetUnrecognizedPdbResidues(Topology const& topIn)
const
{
  Iarray rnums;
  for (int ires = 0; ires != topIn.Nres(); ires++)
  {
    if (!IsRecognizedPdbRes( topIn.Res(ires).Name() ))
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

void Exec_PrepareForLeap::PrintAtomNameMap(const char* title,
                                           std::vector<NameMapType> const& namemap)
{
  mprintf("\t%s:\n", title);
  for (std::vector<NameMapType>::const_iterator it = namemap.begin();
                                                it != namemap.end(); ++it)
  {
    mprintf("\t  %li)", it - namemap.begin());
    for (NameMapType::const_iterator mit = it->begin(); mit != it->end(); ++mit)
      mprintf(" %s:%s", *(mit->first), *(mit->second));
    mprintf("\n");
  }
}

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
        mprintf("\t  %4s -> %s\n", *(mit->first), mit->second.GlycamCode().c_str());
      // DEBUG - print atom name maps
      mprintf("\tRes char to atom map index map:\n");
      for (ResIdxMapType::const_iterator mit = glycam_res_idx_map_.begin(); mit != glycam_res_idx_map_.end(); ++mit)
        mprintf("\t  %s -> %i\n", mit->first.c_str(), mit->second);
      PrintAtomNameMap("Atom name maps", pdb_glycam_name_maps_);
      PrintAtomNameMap("Atom name maps (alpha)", pdb_glycam_name_maps_A_);
      PrintAtomNameMap("Atom name maps (beta)", pdb_glycam_name_maps_B_);
      // DEBUG - print linkage res map
      mprintf("\tLinkage res name map:\n");
      for (NameMapType::const_iterator mit = pdb_glycam_linkageRes_map_.begin(); mit != pdb_glycam_linkageRes_map_.end(); ++mit)
        mprintf("\t  %s -> %s\n", *(mit->first), *(mit->second));
    }
  }

  Iarray pdbResToRemove;
  std::string removeArg = argIn.GetStringKey("remove");
  if (!removeArg.empty()) {
    if (removeArg == "unrecognized") {
      mprintf("\tRemoving unrecognized PDB residues.\n");
      pdbResToRemove = GetUnrecognizedPdbResidues( topIn );
    } else if (removeArg == "isolated") {
      mprintf("\tRemoving unrecognized and isolated PDB residues.\n");
      Iarray unrecognizedPdbRes = GetUnrecognizedPdbResidues( topIn );
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

  // Get sugar mask or default sugar mask
  std::string sugarmaskstr = argIn.GetStringKey("sugarmask");
  if (!sugarmaskstr.empty()) {
    if (!prepare_sugars) {
      mprinterr("Error: Cannot specify 'nosugars' and 'sugarmask'\n");
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
  }
  hasGlycam_ = argIn.hasKey("hasglycam");
  if (hasGlycam_)
    mprintf("\tAssuming sugars already have glycam residue names.\n");
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
  resStat_.assign( topIn.Nres(), UNKNOWN );

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
        resStat_[res-topIn.ResStart()] = VALIDATED;
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
  for (ResStatArray::iterator it = resStat_.begin(); it != resStat_.end(); ++it)
  {
    LeapFxnGroupWarning(topIn, it-resStat_.begin());
    //if ( *it == VALIDATED )
    //  mprintf("\t\t%s VALIDATED\n", topIn.TruncResNameOnumId(it-resStat_.begin()).c_str());
    //else
    //  mprintf("\t\t%s UNKNOWN\n", topIn.TruncResNameOnumId(it-resStat_.begin()).c_str());
    // ----- Warnings --------
    if ( *it == UNKNOWN ) {
      SetType::const_iterator pname = pdb_res_names_.find( topIn.Res(it-resStat_.begin()).Name() );
      if (pname == pdb_res_names_.end())
        mprintf("\t%s%s is an unrecognized name and may not have parameters.\n",
                msg1, topIn.TruncResNameOnumId(it-resStat_.begin()).c_str());
      else
        *it = VALIDATED;
    } else if ( *it == SUGAR_NAME_MISMATCH ) {
        mprintf("\t%s%s sugar anomer type and/or configuration is not consistent with name.\n",
                msg1, topIn.TruncResNameOnumId(it-resStat_.begin()).c_str());
    // ----- Fatal Errors ----
    } else if ( *it == SUGAR_UNRECOGNIZED_LINK_RES ) {
        mprintf("\t%s%s is linked to a sugar but has no sugar-linkage form.\n",
                msg2, topIn.TruncResNameOnumId(it-resStat_.begin()).c_str());
        fatal_errors++;
    } else if ( *it == SUGAR_UNRECOGNIZED_LINKAGE ) {
        mprintf("\t%s%s is a sugar with an unrecognized linkage.\n",
                msg2, topIn.TruncResNameOnumId(it-resStat_.begin()).c_str());
        fatal_errors++;
    } else if ( *it == SUGAR_NO_LINKAGE ) {
        mprintf("\t%s%s is an incomplete sugar with no linkages.\n",
                msg2, topIn.TruncResNameOnumId(it-resStat_.begin()).c_str());
        fatal_errors++;
    } else if ( *it == SUGAR_NO_CHAIN_FOR_LINK ) {
        mprintf("\t%s%s could not identify chain atoms for determining linkages.\n",
                msg2, topIn.TruncResNameOnumId(it-resStat_.begin()).c_str());
        fatal_errors++;
/*    } else if ( *it == SUGAR_MISSING_C1X ) { // TODO should this be a warning
        mprintf("\t%s%s Sugar is missing anomeric carbon substituent.\n",
                msg2, topIn.TruncResNameOnumId(it-resStat_.begin()).c_str());
        fatal_errors++;*/
    } else if ( *it == SUGAR_SETUP_FAILED ) {
        mprintf("\t%s%s Sugar setup failed and could not be identified.\n",
                msg2, topIn.TruncResNameOnumId(it-resStat_.begin()).c_str());
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
