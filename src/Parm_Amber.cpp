#include <cstring> // strncmp
#include <cstdio> // sscanf
#include <cctype> // toupper, isdigit
#include <cstdlib>  // atoi
#include <cmath> // sqrt
#include <algorithm> // std::min, std::max
#include "Parm_Amber.h"
#include "CpptrajStdio.h"
#include "Constants.h" // ELECTOAMBER, AMBERTOELEC
#include "StringRoutines.h" // NoTrailingWhitespace
#include "ExclusionArray.h"

// ---------- Constants and Enumerated types -----------------------------------
const int Parm_Amber::AMBERPOINTERS_ = 31;
/// Enumerated type for FLAG_POINTERS section
/** These variables are part of the POINTERS section of the topology.
  NATOM;    total number of atoms in the system
  NTYPES;   number of AMBER atom types used, max is 60
  NBONH;    number of bonds containing hydrogen
  NBONA;    number of bonds without hydrogen
  NTHETH;   number of angles containing hydrogen
  NTHETA;   number of angles not containing hydrogen
  NPHIH;    number of dihedrals containing hydrogen
  NPHIA;    number of dihedrals not containing hydrogen
  NHPARM;   NOT USED
  NPARM;    1 for LES parm
  NNB;      total number of excluded atoms
  NRES;     total number of residues
  MBONA;    NBONA + number of constraint bonds
  MTHETA;   NTHETA + number of constraint angles
  MPHIA;    NPHIA + number of constraint dihedral angles
  NUMBND;   total number of unique bond types
  NUMANG;   total number of unique angle types
  NPTRA;    total number of unique dihedral types
  NATYP;    number of atom types defined in parameter file
  NPHB;     number of types of hydrogen bonded pair interactions
  IFPERT;   =1 if perturbation info is to be read =0 otherwise
  NBPER;    number of bonds to be perturbed
  NGPER;    number of angles to be perturbed
  NDPER;    number of dihedrals to be perturbed
  MBPER;    num of pert bonds across boundary to non-pert groups
  MGPER;    num of pert angles across boundary to non-pert groups
  MDPER;    num of pert dihedrals across bndry to non-pert groups
  IFBOX;    >0 if periodic box info to be read =0 otherwise
  NMXRS;    number of atoms in the largest residue
  IFCAP;    =1 if CAP option was used in edit, =0 otherwise
  NUMEXTRA; number of extra points (aka lone pairs)
*/
enum topValues {
//0         1       2      3       4       5       6       7      8       9
  NATOM=0,  NTYPES, NBONH, NBONA,  NTHETH, NTHETA, NPHIH,  NPHIA, NHPARM, NPARM,
  NNB,      NRES,   MBONA, MTHETA, MPHIA,  NUMBND, NUMANG, NPTRA, NATYP,  NPHB,
  IFPERT,   NBPER,  NGPER, NDPER,  MBPER,  MGPER,  MDPER,  IFBOX, NMXRS,  IFCAP,
  NEXTRA
};
// FORTRAN format strings
static const char* F10I8 = "%FORMAT(10I8)";
static const char* F20a4 = "%FORMAT(20a4)";
static const char* F5E16 = "%FORMAT(5E16.8)";
static const char* F3I8  = "%FORMAT(3I8)";
static const char* F1a80 = "%FORMAT(1a80)";
static const char* F1I8  = "%FORMAT(1I8)";
/// Constant strings for Amber parm flags and fortran formats. Enumerated by FlagType
const Parm_Amber::ParmFlag Parm_Amber::FLAGS_[] = {
  { "POINTERS",                   F10I8 }, ///< Described above in topValues
  { "ATOM_NAME",                  F20a4 }, ///< Atom names
  { "CHARGE",                     F5E16 }, ///< Atom charges
  { "MASS",                       F5E16 }, ///< Atom masses
  { "RESIDUE_LABEL",              F20a4 }, ///< Residue names
  { "RESIDUE_POINTER",            F10I8 }, ///< Residue boundaries (atoms)
  { "AMBER_ATOM_TYPE",            F20a4 }, ///< Atom types
  { "BONDS_INC_HYDROGEN",         F10I8 }, ///< Bonds to hydrogen
  { "BONDS_WITHOUT_HYDROGEN",     F10I8 }, ///< Bonds not including hydrogen
  { "SOLVENT_POINTERS",           F3I8  },
  { "ATOMS_PER_MOLECULE",         F10I8 },
  { "BOX_DIMENSIONS",             F5E16 },
  { "ATOM_TYPE_INDEX",            F10I8 },
  { "NUMBER_EXCLUDED_ATOMS",      F10I8 },
  { "NONBONDED_PARM_INDEX",       F10I8 },
  { "LENNARD_JONES_ACOEF",        F5E16 },
  { "LENNARD_JONES_BCOEF",        F5E16 },
  { "LENNARD_JONES_CCOEF",        F5E16 },
  { "EXCLUDED_ATOMS_LIST",        F10I8 },
  { "RADII",                      F5E16 },
  { "SCREEN",                     F5E16 },
  { "BOND_FORCE_CONSTANT",        F5E16 },
  { "BOND_EQUIL_VALUE",           F5E16 },
  { "ANGLE_FORCE_CONSTANT",       F5E16 },
  { "ANGLE_EQUIL_VALUE",          F5E16 },
  { "DIHEDRAL_FORCE_CONSTANT",    F5E16 },
  { "DIHEDRAL_PERIODICITY",       F5E16 },
  { "DIHEDRAL_PHASE",             F5E16 },
  { "SCEE_SCALE_FACTOR",          F5E16 },
  { "SCNB_SCALE_FACTOR",          F5E16 },
  { "SOLTY",                      F5E16 },
  { "ANGLES_INC_HYDROGEN",        F10I8 },
  { "ANGLES_WITHOUT_HYDROGEN",    F10I8 },
  { "DIHEDRALS_INC_HYDROGEN",     F10I8 },
  { "DIHEDRALS_WITHOUT_HYDROGEN", F10I8 },
  { "HBOND_ACOEF",                F5E16 },
  { "HBOND_BCOEF",                F5E16 },
  { "HBCUT",                      F5E16 },
  { "TREE_CHAIN_CLASSIFICATION",  F20a4 },
  { "JOIN_ARRAY",                 F10I8 },
  { "IROTAT",                     F10I8 },
  { "ATOMIC_NUMBER",              F10I8 },
  { "TITLE",                      F20a4 },
  { "RADIUS_SET",                 F1a80 },
  { "LES_NTYP",                   F10I8 }, // Number of LES region types
  { "LES_TYPE",                   F10I8 }, // LES type for each atom
  { "LES_FAC",                    F5E16 }, // Scaling factor for typeA * typeB  
  { "LES_CNUM",                   F10I8 }, // Copy number for each atom; 0==in all
  { "LES_ID",                     F10I8 }, // LES region ID
  { "CAP_INFO",                   F10I8 },
  { "CAP_INFO2",                  F5E16 },
  { "IPOL",                       F1I8  }, // 0 for fixed charge, 1 for polarizable
  { "POLARIZABILITY",             F5E16 }, // Hold atom polarazabilities in Ang^3
  // CHAMBER parameters
  { "CTITLE",                             "%FORMAT(a80)" },
  { "CHARMM_UREY_BRADLEY_COUNT",          "%FORMAT(2I8)" }, // # UB terms and types
  { "CHARMM_UREY_BRADLEY",                F10I8          }, // UB: 2 atoms and param index
  { "CHARMM_UREY_BRADLEY_FORCE_CONSTANT", F5E16          },
  { "CHARMM_UREY_BRADLEY_EQUIL_VALUE",    F5E16          },
  { "CHARMM_NUM_IMPROPERS",               F10I8          }, // # improper terms
  { "CHARMM_IMPROPERS",                   F10I8          }, // Imp: 4 atoms and param index
  { "CHARMM_NUM_IMPR_TYPES",              "%FORMAT(I8)"  }, // # improper types
  { "CHARMM_IMPROPER_FORCE_CONSTANT",     F5E16          },
  { "CHARMM_IMPROPER_PHASE",              F5E16          },
  { "LENNARD_JONES_14_ACOEF",             "%FORMAT(3E24.16)" },
  { "LENNARD_JONES_14_BCOEF",             "%FORMAT(3E24.16)" },
  { "CHARMM_CMAP_COUNT",                  "%FORMAT(2I8)" }, // # CMAP terms, # unique CMAP params
  { "CHARMM_CMAP_RESOLUTION",             "%FORMAT(20I4)"}, // # steps along each Phi/Psi CMAP axis
  { "CHARMM_CMAP_PARAMETER_",             "%FORMAT(8(F9.5))"}, // CMAP grid
  { "CHARMM_CMAP_INDEX",                  "%FORMAT(6I8)" }, // Atom i,j,k,l,m of cross term and idx
  { "FORCE_FIELD_TYPE",                   "%FORMAT(i2,a78)"},// NOTE: Cannot use with SetFortranType
  // PDB extra info
  { "RESIDUE_NUMBER", "%FORMAT(20I4)" },                // PDB residue number
  { "RESIDUE_CHAINID", F20a4 },                         // PDB chain ID
  { "RESIDUE_ICODE", F20a4 },                           // PDB residue insertion code
  { "ATOM_ALTLOC", F20a4 },                             // PDB atom alt location indicator FIXME: format is guess
  { "ATOM_BFACTOR", "%FORMAT(10F8.2)"},                 // PDB atom B-factors
  { "ATOM_OCCUPANCY", "%FORMAT(10F8.2)"},               // PDB atom occupancies
  { "ATOM_NUMBER", F10I8},                              // PDB original atom serial #s
  // CHARMM CMAP
  { "CMAP_COUNT",                  "%FORMAT(2I8)" },    // # CMAP terms, # unique CMAP params
  { "CMAP_RESOLUTION",             "%FORMAT(20I4)"},    // # steps along each Phi/Psi CMAP axis
  { "CMAP_PARAMETER_",             "%FORMAT(8(F9.5))"}, // CMAP grid
  { "CMAP_INDEX",                  "%FORMAT(6I8)" },    // Atom i,j,k,l,m of cross term and idx
  { 0, 0 }
};

const double Parm_Amber::ELECTOCHAMBER_ = sqrt(332.0716);

const double Parm_Amber::CHAMBERTOELEC_ = 1.0 / ELECTOCHAMBER_;

// -----------------------------------------------------------------------------
// CONSTRUCTOR
Parm_Amber::Parm_Amber() :
  ptype_(OLDPARM),
  elec_to_parm_(Constants::ELECTOAMBER),
  parm_to_elec_(Constants::AMBERTOELEC),
  numLJparm_(0),
  SCEE_set_(false),
  SCNB_set_(false),
  atProblemFlag_(false),
  N_impropers_(0),
  N_impTerms_(0),
  writeChamber_(true),
  writeEmptyArrays_(false),
  writePdbInfo_(true)
{
  UB_count_[0] = 0;
  UB_count_[1] = 0;
}

/** It seems that %Xi in sscanf() can get confused when integers take up
  * the entire X characters. Not sure if this is a system bug or not, but
  * this function will first ensure only 6 characters are read before
  * converting to integer.
  * \param ptr Input string to parse. Must not be NULL.
  * \param IVALS Output integer array
  * \return Number of values actually read.
  */
static inline int Get12I6(const char* ptr, int* IVALS) {
  char SVALS[12][7];
  int nvals = sscanf(ptr, "%6c%6c%6c%6c%6c%6c%6c%6c%6c%6c%6c%6c",
                     SVALS[0], SVALS[1], SVALS[2], SVALS[3], SVALS[ 4], SVALS[ 5],
                     SVALS[6], SVALS[7], SVALS[8], SVALS[9], SVALS[10], SVALS[11]);
  // If less than 12 values read the line was short, so the final value will be
  // a newline. Ignore that.
  if (nvals > 0 && nvals < 12) nvals--;
  for (int i = 0; i < nvals; i++) {
    // Only allow integers. Right-aligned, so check final character.
    if (!isdigit(SVALS[i][5])) return 0;
    SVALS[i][6] = '\0';
    //mprintf("DEBUG: %2i = %6s\n", i, SVALS[i]);
    IVALS[i] = atoi( SVALS[i] );
  }
  return nvals;
}

// Parm_Amber::ID_ParmFormat()
bool Parm_Amber::ID_ParmFormat(CpptrajFile& fileIn) {
  int iamber[12];
  char lineBuf[BUF_SIZE];
  // Assumes already set up for READ
  if (fileIn.OpenFile()) return false;
  fileIn.Gets(lineBuf, BUF_SIZE);
  // Check for %VERSION
  if (strncmp(lineBuf,"%VERSION",8)==0) {
    fileIn.Gets(lineBuf, BUF_SIZE);
    // Check for %FLAG
    if (strncmp(lineBuf,"%FLAG",5)==0) {
      if (debug_>0) mprintf("  AMBER TOPOLOGY file\n");
      ptype_ = NEWPARM;
      fileIn.CloseFile();
      return true;
    }
  } else {
    // Since %VERSION not present, If first line is 81 bytes and the second 
    // line has 12 numbers in 12I6 format, assume old-style Amber topology
    // NOTE: Could also be less than 81? Only look for 12 numbers?
    int line1size = (int)strlen(lineBuf);
    if (line1size == (81 + fileIn.IsDos())) {
      fileIn.Gets(lineBuf, BUF_SIZE);
      if ( Get12I6(lineBuf, iamber) == 12)
      {
        if (debug_>0) mprintf("  AMBER TOPOLOGY, OLD FORMAT\n");
        ptype_ = OLDPARM;
        fileIn.CloseFile();
        return true;
      }
    }
  }
  fileIn.CloseFile();
  return false;
}

// -----------------------------------------------------------------------------
// Parm_Amber::ReadParm()
int Parm_Amber::ReadParm(FileName const& fname, Topology& TopIn ) {
  if (file_.OpenRead( fname )) return 1;
  elec_to_parm_ = Constants::ELECTOAMBER;
  parm_to_elec_ = Constants::AMBERTOELEC;
  int err = 0;
  if (ptype_ == OLDPARM)
    err = ReadOldParm( TopIn );
  else
    err = ReadNewParm( TopIn );
  if (err != 0) return 1;
  // Determine Atom elements
  if (atomicNums_.empty()) {
    mprintf("\tThis Amber topology does not include atomic numbers.\n"
            "\tAssigning elements from atom masses/names.\n");
    atomicNums_.assign( values_[NATOM], 0 );
  }
  for (int idx = 0; idx != values_[NATOM]; idx++)
    TopIn.SetAtom(idx).DetermineElement( atomicNums_[idx] );
  // Set Atom residue numbers
  for (Topology::res_iterator res = TopIn.ResStart(); res != TopIn.ResEnd(); ++res)
    for (int at = res->FirstAtom(); at != res->LastAtom(); ++at)
      TopIn.SetAtom(at).SetResNum( res - TopIn.ResStart() );
  // If SCEE/SCNB not read, set defaults
  if (!SCEE_set_) {
    mprintf("\tNo SCEE section: setting Amber default (1.2)\n");
    for (unsigned int idx = 0; idx != TopIn.DihedralParm().size(); idx++)
      TopIn.SetDihedralParm(idx).SetSCEE( 1.2 );
  }
  if (!SCNB_set_) {
    mprintf("\tNo SCNB section: setting Amber default (2.0)\n");
    for (unsigned int idx = 0; idx != TopIn.DihedralParm().size(); idx++)
      TopIn.SetDihedralParm(idx).SetSCNB( 2.0 );
  }
  // Check box info
  if (values_[IFBOX] > 0) {
    if (!parmbox_.HasBox()) {
      if (ptype_ != CHAMBER) mprintf("Warning: Prmtop missing Box information.\n");
      // Check for IFBOX/BoxType mismatch
      if (values_[IFBOX]==2 && parmbox_.CellShape() != Box::OCTAHEDRAL) {
        mprintf("Warning: Amber Parm Box should be Truncated Octahedron (ifbox==2)\n"
                "         but BOX_DIMENSIONS indicate %s - may cause imaging problems.\n",
                parmbox_.CellShapeName());
      }
    }
  }

  TopIn.SetParmBox( parmbox_ );

  // DEBUG
  //for (Topology::atom_iterator atom = TopIn.begin(); atom != TopIn.end(); ++atom)
  //  mprintf("%u: %s Res %i\n", atom-TopIn.begin(), atom->c_str(), atom->ResNum());
  //for (Topology::res_iterator res = TopIn.ResStart(); res != TopIn.ResEnd(); ++res)
  //  mprintf("%u: %s FirstAt=%i EndAt=%i Original=%i\n", res-TopIn.ResStart(),
  //          res->c_str(), res->FirstAtom(), res->LastAtom(), res->OriginalResNum());
  return 0;
}

// Parm_Amber::ReadOldParm()
int Parm_Amber::ReadOldParm(Topology& TopIn) {
  mprintf("\tReading old (<v7) Amber Topology file.\n");
  std::string title = NoTrailingWhitespace( file_.GetLine() );
  TopIn.SetParmName( title, file_.Filename() );
  int Npointers = 0; // Number of pointers not known, but definitely no NEXTRA etc
  const FortranData DBL(FDOUBLE, 5, 16, 0);
  const FortranData INT(FINT, 12, 6, 0);
  const FortranData CHAR(FCHAR, 20, 4, 0);
  if ( ReadPointers( Npointers, TopIn, INT ) ) return 1;
  if ( ReadAtomNames( TopIn, CHAR ) ) return 1;
  if ( ReadAtomCharges( TopIn, DBL ) ) return 1;
  if ( ReadAtomicMass( TopIn, DBL ) ) return 1;
  if ( ReadAtomTypeIndex( TopIn, INT ) ) return 1;
  // Skip past NUMEX
  if (SetupBuffer(F_NUMEX, values_[NATOM], INT)) return 1;
  if ( ReadNonbondIndices(TopIn, INT) ) return 1;
  if ( ReadResidueNames(TopIn, CHAR) ) return 1;
  if ( ReadResidueAtomNums(TopIn, INT) ) return 1;
  if ( ReadBondRK(TopIn, DBL) ) return 1;
  if ( ReadBondREQ(TopIn, DBL) ) return 1;
  if ( ReadAngleTK(TopIn, DBL) ) return 1;
  if ( ReadAngleTEQ(TopIn, DBL) ) return 1;
  if ( ReadDihedralPK(TopIn, DBL) ) return 1;
  if ( ReadDihedralPN(TopIn, DBL) ) return 1;
  if ( ReadDihedralPHASE(TopIn, DBL) ) return 1;
  // Skip past SOLTY
  if (SetupBuffer(F_SOLTY, values_[NATYP], DBL)) return 1;
  if ( ReadLJA(TopIn, DBL) ) return 1;
  if ( ReadLJB(TopIn, DBL) ) return 1;
  if ( ReadBondsH(TopIn, INT) ) return 1;
  if ( ReadBonds(TopIn, INT) ) return 1;
  if ( ReadAnglesH(TopIn, INT) ) return 1;
  if ( ReadAngles(TopIn, INT) ) return 1;
  if ( ReadDihedralsH(TopIn, INT) ) return 1;
  if ( ReadDihedrals(TopIn, INT) ) return 1;
  // Skip past EXCLUDE
  if (SetupBuffer(F_EXCLUDE, values_[NNB], INT)) return 1;
  if ( ReadAsol(TopIn, DBL) ) return 1;
  if ( ReadBsol(TopIn, DBL) ) return 1;
  if ( ReadHBcut(TopIn, DBL) ) return 1;
  if ( ReadAtomTypes(TopIn, CHAR) ) return 1;
  if ( ReadItree(TopIn, CHAR) ) return 1;
  if ( ReadJoin(TopIn, INT) ) return 1;
  if ( ReadIrotat(TopIn, INT) ) return 1;
  // Solvent info
  if (values_[IFBOX] > 0) {
    // Read SOLVENT_POINTERS, only need number of molecules (2nd value).
    if (SetupBuffer(F_SOLVENT_POINTER, 3, INT)) return 1;
    file_.NextElement(); // Final solute residue
    int nmolecules = atoi(file_.NextElement());
    // Skip past ATOMS_PER_MOL
    if (SetupBuffer(F_ATOMSPERMOL, nmolecules, INT)) return 1;
    if ( ReadBox(DBL) ) return 1;
  }
  // Water cap info
  if (values_[IFCAP]) {
    if ( ReadCapInfo(TopIn, INT) ) return 1;
    if ( ReadCapInfo2(TopIn, DBL) ) return 1;
  }
  // LES info
  if (values_[NPARM] == 1) {
    if ( ReadLESntyp(TopIn, INT) ) return 1;
    if ( ReadLEStypes(TopIn, INT) ) return 1;
    if ( ReadLESfac(TopIn, DBL) ) return 1;
    if ( ReadLEScnum(TopIn, INT) ) return 1;
    if ( ReadLESid(TopIn, INT) ) return 1;
  }

  return 0;
}

static inline bool IsFLAG(const char* ptr) {
  return (ptr[1] == 'F' && ptr[2] == 'L' && ptr[3] == 'A' && ptr[4] == 'G');
}

// Parm_Amber::ReadNewParm()
int Parm_Amber::ReadNewParm(Topology& TopIn) {
  FortranData FMT; // Used to hold fortran format from format string
  const char* ptr = file_.NextLine();
  if (ptr == 0) {
    mprinterr("Error: Unexpected EOF encountered.\n");
    return 1;
  }
  // Main loop for reading file.
  while (ptr != 0) {
    //mprintf("DEBUG: LINE: %s", ptr);
    if ( ptr[0] == '%' ) {
      if (ptr[1] == 'V' && ptr[2] == 'E' && ptr[3] == 'R') {
        // %VERSION line. Skip it.
        ptr = file_.NextLine();
      } else if (IsFLAG(ptr)) {
        // %FLAG <type> line. Determine the flag type.
        std::string flagType = NoTrailingWhitespace(ptr+6);
        //mprintf("DEBUG: Flag type: %s\n", flagType.c_str());
        int flagIdx = -1;
        if ( flagType.compare(0, 22, "CHARMM_CMAP_PARAMETER_") == 0 ) {
          // Special case. This flag has a 2 digit extension.
          flagIdx = (int)F_CHM_CMAPP;
        } else if ( flagType.compare(0, strlen("CMAP_PARAMETER_"), "CMAP_PARAMETER_") == 0 ) {
          // Special case. This flag has a 2 digit extension.
          flagIdx = (int)F_CMAPP;
        } else {
          for (ParmPtr P = FLAGS_; P->Flag != 0; ++P) {
            if (flagType.compare(P->Flag) == 0) {
              flagIdx = (int)(P - FLAGS_);
              break;
            }
          }
        }
        // Read the format line. Do this even if FLAG is not recognized.
        if (ReadFormatLine(FMT)) return 1;
        // Process the FLAG
        if (flagIdx == -1) {
          mprintf("Warning: Amber topology flag '%s' is unrecognized and will be skipped.\n",
                  flagType.c_str());
          ptr = SkipToNextFlag();
        } else {
          int err = 0;
          FlagType ftype = (FlagType)flagIdx;
          switch (ftype) {
            case F_CTITLE: ptype_ = CHAMBER; // Fall through to F_TITLE
                           elec_to_parm_ = ELECTOCHAMBER_;
                           parm_to_elec_ = CHAMBERTOELEC_;
            case F_TITLE:     err = ReadTitle(TopIn); break;
            case F_POINTERS:  err = ReadPointers(AMBERPOINTERS_, TopIn, FMT); break;
            case F_NAMES:     err = ReadAtomNames(TopIn, FMT); break;
            case F_CHARGE:    err = ReadAtomCharges(TopIn, FMT); break;
            case F_ATOMICNUM: err = ReadAtomicNum(FMT); break;
            case F_MASS:      err = ReadAtomicMass(TopIn, FMT); break;
            case F_ATYPEIDX:  err = ReadAtomTypeIndex(TopIn, FMT); break;
            // NOTE: CPPTRAJ sets up its own exclusion list so reading this is skipped.
            case F_NUMEX: ptr = SkipToNextFlag(); break;
            case F_NB_INDEX:  err = ReadNonbondIndices(TopIn, FMT); break;
            case F_RESNAMES:  err = ReadResidueNames(TopIn, FMT); break;
            case F_RESNUMS:   err = ReadResidueAtomNums(TopIn, FMT); break;
            case F_BONDRK:    err = ReadBondRK(TopIn, FMT); break;
            case F_BONDREQ:   err = ReadBondREQ(TopIn, FMT); break;
            case F_ANGLETK:   err = ReadAngleTK(TopIn, FMT); break;
            case F_ANGLETEQ:  err = ReadAngleTEQ(TopIn, FMT); break;
            case F_DIHPK:     err = ReadDihedralPK(TopIn, FMT); break;
            case F_DIHPN:     err = ReadDihedralPN(TopIn, FMT); break;
            case F_DIHPHASE:  err = ReadDihedralPHASE(TopIn, FMT); break;
            case F_SCEE:      err = ReadDihedralSCEE(TopIn, FMT); break;
            case F_SCNB:      err = ReadDihedralSCNB(TopIn, FMT); break;
            case F_SOLTY: ptr = SkipToNextFlag(); break;
            case F_LJ_A:      err = ReadLJA(TopIn, FMT); break;
            case F_LJ_B:      err = ReadLJB(TopIn, FMT); break;
            case F_LJ_C:      err = ReadLJC(TopIn, FMT); break;
            case F_BONDSH:    err = ReadBondsH(TopIn, FMT); break;
            case F_BONDS:     err = ReadBonds(TopIn, FMT); break;
            case F_ANGLESH:   err = ReadAnglesH(TopIn, FMT); break;
            case F_ANGLES:    err = ReadAngles(TopIn, FMT); break;
            case F_DIHH:      err = ReadDihedralsH(TopIn, FMT); break;
            case F_DIH:       err = ReadDihedrals(TopIn, FMT); break;
            case F_EXCLUDE: ptr = SkipToNextFlag(); break;
            case F_ASOL:      err = ReadAsol(TopIn, FMT); break;
            case F_BSOL:      err = ReadBsol(TopIn, FMT); break;
            case F_HBCUT:     err = ReadHBcut(TopIn, FMT); break;
            case F_TYPES:     err = ReadAtomTypes(TopIn, FMT); break;
            case F_ITREE:     err = ReadItree(TopIn, FMT); break;
            case F_JOIN:      err = ReadJoin(TopIn, FMT); break;
            case F_IROTAT:    err = ReadIrotat(TopIn, FMT); break;
            case F_SOLVENT_POINTER: ptr = SkipToNextFlag(); break;
            case F_ATOMSPERMOL: ptr = SkipToNextFlag(); break;
            case F_PARMBOX:   err = ReadBox(FMT); break;
            case F_CAP_INFO:  err = ReadCapInfo(TopIn, FMT); break;
            case F_CAP_INFO2: err = ReadCapInfo2(TopIn, FMT); break;
            case F_RADSET:    err = ReadGBradiiSet(TopIn); break;
            case F_RADII:     err = ReadGBradii(TopIn, FMT); break;
            case F_SCREEN:    err = ReadGBscreen(TopIn, FMT); break;
            case F_IPOL:      err = ReadIpol(TopIn, FMT); break;
            case F_POLAR:     err = ReadPolar(TopIn, FMT); break;
            // Extra PDB Info
            case F_PDB_RES:   err = ReadPdbRes(TopIn, FMT); break;
            case F_PDB_CHAIN: err = ReadPdbChainID(TopIn, FMT); break;
            case F_PDB_ICODE: err = ReadPdbIcode(TopIn, FMT); break;
            case F_PDB_ALT:   err = ReadPdbAlt(TopIn, FMT); break;
            case F_PDB_BFAC:  err = ReadPdbBfactor(TopIn, FMT); break;
            case F_PDB_OCC:   err = ReadPdbOccupancy(TopIn, FMT); break;
            case F_PDB_NUM:   err = ReadPdbNumbers(TopIn, FMT); break;
            // CHAMBER
            case F_FF_TYPE:   err = ReadChamberFFtype(TopIn, FMT); break;
            case F_CHM_UBC:   err = ReadChamberUBCount(TopIn, FMT); break;
            case F_CHM_UB:    err = ReadChamberUBTerms(TopIn, FMT); break;
            case F_CHM_UBFC:  err = ReadChamberUBFC(TopIn, FMT); break;
            case F_CHM_UBEQ:  err = ReadChamberUBEQ(TopIn, FMT); break;
            case F_CHM_NIMP:  err = ReadChamberNumImpropers(TopIn, FMT); break;
            case F_CHM_NIMPT: err = ReadChamberNumImpTerms(TopIn, FMT); break;
            case F_CHM_IMP:   err = ReadChamberImpropers(TopIn, FMT); break;
            case F_CHM_IMPFC: err = ReadChamberImpFC(TopIn, FMT); break;
            case F_CHM_IMPP:  err = ReadChamberImpPHASE(TopIn, FMT); break;
            case F_LJ14A:     err = ReadChamberLJ14A(TopIn, FMT); break;
            case F_LJ14B:     err = ReadChamberLJ14B(TopIn, FMT); break;
            // LES
            case F_LES_NTYP:  err = ReadLESntyp(TopIn, FMT); break;
            case F_LES_TYPE:  err = ReadLEStypes(TopIn, FMT); break;
            case F_LES_FAC:   err = ReadLESfac(TopIn, FMT); break;
            case F_LES_CNUM:  err = ReadLEScnum(TopIn, FMT); break;
            case F_LES_ID:    err = ReadLESid(TopIn, FMT); break;
            // CMAP
            case F_CHM_CMAPC: // fallthrough
            case F_CMAPC:     err = ReadCmapCounts(ftype, FMT); break;
            case F_CHM_CMAPR: // fallthrough
            case F_CMAPR:     err = ReadCmapRes(ftype, TopIn, FMT); break;
            case F_CHM_CMAPP: // fallthrough
            case F_CMAPP:     err = ReadCmapGrid(ftype, flagType.c_str(), TopIn, FMT); break;
            case F_CHM_CMAPI: // fallthrough
            case F_CMAPI:     err = ReadCmapTerms(ftype, TopIn, FMT); break;
            // Sanity check
            default:
              mprinterr("Internal Error: Unhandled FLAG '%s'.\n", flagType.c_str());
              return 1;
          }
          if (err != 0) {
            mprinterr("Error: Reading format FLAG '%s'\n", flagType.c_str());
            return 1;
          }
          if (atProblemFlag_) {
            // We hit a problematic flag and need to read past it.
            ptr = SkipToNextFlag();
            atProblemFlag_ = false;
          }
        }
      } else {
        // Unknown '%' tag. Read past it.
        ptr = file_.NextLine();
      }  
    } else {
      // Unknown line. Read past it.
      ptr = file_.NextLine();
    }
  }

  file_.CloseFile();
  return 0;
}

// Parm_Amber::SkipToNextFlag()
const char* Parm_Amber::SkipToNextFlag() {
  const char* ptr = file_.NextLine();
  while (ptr != 0 && !IsFLAG(ptr)) ptr = file_.NextLine();
  return ptr;
}

// Parm_Amber::ReadFormatLine()
int Parm_Amber::ReadFormatLine(FortranData& FMT) {
  // Read next line; can be either a COMMENT or FORMAT. If COMMENT, 
  // read past until you get to the FORMAT line.
  const char* ptr = file_.NextLine();
  if (ptr == 0) {
    mprinterr("Error: Unexpected EOF in Amber Topology when looking for FORMAT.\n");
    return 1;
  }
  while ( ptr != 0 && strncmp(ptr, "%FORMAT", 7) !=0) {
    ptr = file_.NextLine();
    // Sanity check
    if (IsFLAG(ptr)) {
      mprinterr("Error: Missing FORMAT line.\n");
      return 1;
    }
  }
  if (debug_>1) mprintf("DEBUG: Format line [%s]\n", ptr);
  // Parse format string
  if (FMT.ParseFortranFormat( ptr )) return 1;
  
  return 0;
}

/** If an error occurred at FLAG, this can be used to reset the file
  * and scan past the problem flag.
  */
void Parm_Amber::ResetFileToFlag(FlagType currFlag) {
  mprintf("Info: Scanning past problematic flag %s\n", FLAGS_[currFlag].Flag);
  file_.Rewind();
  const char* ptr = file_.NextLine();
  // TODO trap null?
  atProblemFlag_ = false;
  while (ptr != 0) {
    if ( ptr[0] == '%' && IsFLAG(ptr) ) {
      // %FLAG <type> line. Determine the flag type.
      std::string flagType = NoTrailingWhitespace(ptr+6);
      if (flagType.compare( FLAGS_[currFlag].Flag ) == 0) {
        // Problem flag found. Set the problemFlag variable so the parser knows to skip ahead.
        atProblemFlag_ = true;
        break;
      }
    }
    ptr = file_.NextLine();
  }
}

/** Print a warning that a problem was encountered reading an element
  * of the specified FLAG. Sets the atProblemFlag_ variable to true
  * and resets the file to the bad FLAG so it can be skipped.
  */
void Parm_Amber::ProblemFlagWarning(FlagType currFlag, unsigned int idx, unsigned int nExpected) {
   mprintf("Warning: Bad conversion detected: %s\n", FLAGS_[currFlag].Flag);
   mprintf("Warning: Issue reading element %u of %u\n", idx+1, nExpected);
   atProblemFlag_ = true;
   ResetFileToFlag( currFlag );
}

/** Attempt to convert next element in file buffer to double, with
  * some error checking.
  */
double Parm_Amber::FileBufferToDouble(FlagType currFlag, unsigned int idx, unsigned int nExpected) {
  char* endptr;
  const char* elt = file_.NextElement();
  double dval = strtod( elt, &endptr );
  if (elt == endptr) {
    ProblemFlagWarning(currFlag, idx, nExpected);
    return 0.0;
  }
  return dval;
}

// Parm_Amber::ReadTitle()
int Parm_Amber::ReadTitle(Topology& TopIn) {
  std::string title = NoTrailingWhitespace( file_.GetLine() );
  if (debug_ > 0)
    mprintf("\tAmberParm Title: \"%s\"\n",title.c_str());
  TopIn.SetParmName( title, file_.Filename() );
  if (file_.NextLine() == 0) return 1; // Advance to next line
  return 0;
}

// Parm_Amber::ReadPointers()
int Parm_Amber::ReadPointers(int Npointers, Topology& TopIn, FortranData const& FMT) {
  if (Npointers > 0) {
    // New >= version 7 topology: number of pointers is consistent.
    file_.SetupFrameBuffer( Npointers, FMT.Width(), FMT.Ncols() );
    if (file_.ReadFrame()) return 1;
    values_.reserve(Npointers);
    for (int idx = 0; idx != Npointers; idx++)
      values_.push_back( atoi( file_.NextElement() ) );
  } else {
    // Old version. 3 lines, but number of pointers is not consistent.
    int IVALS[12];
    int nPointers = 0;
    for (int line = 0; line != 3; line++) {
      const char* ptr = file_.NextLine();
      if (ptr == 0) return 1;
      // Old pointers format is 12I6
      int nvals = Get12I6(ptr, IVALS);
      nPointers += nvals;
      // First two lines should always have 12 values.
      if (line < 2 && nvals < 12) {
        mprinterr("Error: In old topology file, not enough POINTERS (%i).\n", nPointers);
        return 1;
      }
      for (int ip = 0; ip != nvals; ip++)
        values_.push_back( IVALS[ip] );
    }
    // Make sure there are at least AMBERPOINTERS_ pointers.
    if (values_.size() < AMBERPOINTERS_)
      values_.resize( AMBERPOINTERS_, 0 );
  }
  if (debug_ > 4) {
    mprintf("DEBUG: POINTERS\n");
    for (Iarray::const_iterator it = values_.begin(); it != values_.end(); ++it)
      mprintf("%li\t%i\n", it-values_.begin(), *it);
  }
  TopIn.Resize( Topology::Pointers(values_[NATOM], values_[NRES],
                                   values_[NUMBND], values_[NUMANG], values_[NPTRA]) );

  if (values_[IFPERT] > 0)
    mprintf("Warning: '%s' contains perturbation information.\n"
            "Warning:  Cpptraj currently does not read or write perturbation information.\n",
            file_.Filename().base());

  TopIn.SetNonbond().SetupLJforNtypes( values_[NTYPES] );
  numLJparm_ = TopIn.Nonbond().NBarray().size();
  TopIn.SetNonbond().SetNHBterms( values_[NPHB] );
  TopIn.SetNatyp( values_[NATYP] );
  return 0;
}

// Parm_Amber::SetupBuffer()
int Parm_Amber::SetupBuffer(FlagType ftype, int nvals, FortranData const& FMT) {
  if (values_.empty()) {
    mprinterr("Error: Flag '%s' encountered before POINTERS.\n", FLAGS_[ftype].Flag);
    return 1;
  }
  if (nvals > 0) {
    if (debug_>0) mprintf("DEBUG: Set up buffer for '%s', %i vals.\n", FLAGS_[ftype].Flag, nvals);
    file_.SetupFrameBuffer( nvals, FMT.Width(), FMT.Ncols() );
    if (file_.ReadFrame()) return 1;
    if (debug_ > 5) {
      mprintf("DEBUG: '%s':\n", FLAGS_[ftype].Flag);
      if (debug_ > 6) mprintf("FileBuffer=[%s]", file_.Buffer());
    }
  } else {
    if (debug_>5) mprintf("DEBUG: No values for flag '%s'\n", FLAGS_[ftype].Flag);
    // Read blank line
    file_.NextLine();
  }
  return 0;
}

// Parm_Amber::ReadAtomNames()
int Parm_Amber::ReadAtomNames(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_NAMES, values_[NATOM], FMT)) return 1;
  for (int idx = 0; idx != values_[NATOM]; idx++)
    TopIn.SetAtom(idx).SetName( NameType(file_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadAtomCharges()
/** Read atomic charges. Convert units to elec. */
int Parm_Amber::ReadAtomCharges(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_CHARGE, values_[NATOM], FMT)) return 1;
  for (int idx = 0; idx != values_[NATOM]; idx++)
    TopIn.SetAtom(idx).SetCharge( atof(file_.NextElement()) * parm_to_elec_ );
  return 0;
}

// Parm_Amber::ReadAtomicNum()
/** Read atomic numbers to a temporary array. This is done because some
  * topology files do not have atomic number information.
  */
int Parm_Amber::ReadAtomicNum(FortranData const& FMT) {
  if (SetupBuffer(F_ATOMICNUM, values_[NATOM], FMT)) return 1;
  for (int idx = 0; idx != values_[NATOM]; idx++)
    atomicNums_.push_back( atoi(file_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadAtomicMass()
int Parm_Amber::ReadAtomicMass(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_MASS, values_[NATOM], FMT)) return 1;
  for (int idx = 0; idx != values_[NATOM]; idx++)
    TopIn.SetAtom(idx).SetMass( atof(file_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadAtomTypeIndex()
/** Read atom type indices. Shift by -1 to match CPPTRAJ internal indexing.
  * NOTE: If ever used, shift atom #s in excludedAtoms by -1 so they start from 0
  */
int Parm_Amber::ReadAtomTypeIndex(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_ATYPEIDX, values_[NATOM], FMT)) return 1;
  for (int idx = 0; idx != values_[NATOM]; idx++)
    TopIn.SetAtom(idx).SetTypeIndex( atoi(file_.NextElement()) - 1 );
  return 0;
}

// Parm_Amber::ReadNonbondIndices()
int Parm_Amber::ReadNonbondIndices(Topology& TopIn, FortranData const& FMT) {
  int nvals = values_[NTYPES]*values_[NTYPES];
  if (SetupBuffer(F_NB_INDEX, nvals, FMT)) return 1;
  for (int idx = 0; idx != nvals; idx++)
  {
    // Shift positive indices in NONBONDED index array by -1.
    int nbidx = atoi(file_.NextElement());
    if (nbidx > 0)
      nbidx -= 1;
    TopIn.SetNonbond().SetNbIdx(idx, nbidx);
  }
  return 0;
}

// Parm_Amber::ReadResidueNames()
int Parm_Amber::ReadResidueNames(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_RESNAMES, values_[NRES], FMT)) return 1;
  for (int idx = 0; idx != values_[NRES]; idx++)
    TopIn.SetRes(idx).SetName( NameType(file_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadResidueAtomNums()
int Parm_Amber::ReadResidueAtomNums(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_RESNUMS, values_[NRES], FMT)) return 1;
  // Each entry in the array is the start atom (indexed from 1) for
  // the corresponding residue. Do all but the final residue.
  int numRes = values_[NRES];
  int firstAtom = atoi(file_.NextElement()) - 1;
  for (int idx = 1; idx < numRes; idx++) {
    Residue& currentRes = TopIn.SetRes(idx-1);
    int lastAtom = atoi(file_.NextElement()) - 1;
    currentRes.SetFirstAtom( firstAtom );
    currentRes.SetLastAtom( lastAtom );
    currentRes.SetOriginalNum( idx );
    firstAtom = lastAtom;
  }
  // Do the final residue
  TopIn.SetRes(numRes-1).SetFirstAtom( firstAtom );
  TopIn.SetRes(numRes-1).SetLastAtom( values_[NATOM] );
  TopIn.SetRes(numRes-1).SetOriginalNum( numRes );
  return 0;
}

// Parm_Amber::ReadBondRK()
int Parm_Amber::ReadBondRK(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_BONDRK, values_[NUMBND], FMT)) return 1;
  for (int idx = 0; idx != values_[NUMBND]; idx++)
    TopIn.SetBondParm(idx).SetRk( atof(file_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadBondREQ()
int Parm_Amber::ReadBondREQ(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_BONDREQ, values_[NUMBND], FMT)) return 1;
  for (int idx = 0; idx != values_[NUMBND]; idx++)
    TopIn.SetBondParm(idx).SetReq( atof(file_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadAngleTK()
int Parm_Amber::ReadAngleTK(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_ANGLETK, values_[NUMANG], FMT)) return 1;
  for (int idx = 0; idx != values_[NUMANG]; idx++)
    TopIn.SetAngleParm(idx).SetTk( atof(file_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadAngleTEQ()
int Parm_Amber::ReadAngleTEQ(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_ANGLETEQ, values_[NUMANG], FMT)) return 1;
  for (int idx = 0; idx != values_[NUMANG]; idx++)
    TopIn.SetAngleParm(idx).SetTeq( atof(file_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadDihedralPK()
int Parm_Amber::ReadDihedralPK(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_DIHPK, values_[NPTRA], FMT)) return 1;
  for (int idx = 0; idx != values_[NPTRA]; idx++)
    TopIn.SetDihedralParm(idx).SetPk( atof(file_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadDihedralPN()
int Parm_Amber::ReadDihedralPN(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_DIHPN, values_[NPTRA], FMT)) return 1;
  for (int idx = 0; idx != values_[NPTRA]; idx++)
    TopIn.SetDihedralParm(idx).SetPn( atof(file_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadDihedralPHASE()
int Parm_Amber::ReadDihedralPHASE(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_DIHPHASE, values_[NPTRA], FMT)) return 1;
  for (int idx = 0; idx != values_[NPTRA]; idx++)
    TopIn.SetDihedralParm(idx).SetPhase( atof(file_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadDihedralSCEE()
int Parm_Amber::ReadDihedralSCEE(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_SCEE, values_[NPTRA], FMT)) return 1;
  for (int idx = 0; idx != values_[NPTRA]; idx++)
    TopIn.SetDihedralParm(idx).SetSCEE( atof(file_.NextElement()) );
  SCEE_set_ = true;
  return 0;
}

// Parm_Amber::ReadDihedralSCNB()
int Parm_Amber::ReadDihedralSCNB(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_SCNB, values_[NPTRA], FMT)) return 1;
  for (int idx = 0; idx != values_[NPTRA]; idx++)
    TopIn.SetDihedralParm(idx).SetSCNB( atof(file_.NextElement()) );
  SCNB_set_ = true;
  return 0;
}

// Parm_Amber::ReadLJA()
int Parm_Amber::ReadLJA(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_LJ_A, numLJparm_, FMT)) return 1;
  for (int idx = 0; idx != numLJparm_; idx++)
  {
    TopIn.SetNonbond().SetLJ(idx).SetA( FileBufferToDouble(F_LJ_A, idx, numLJparm_) );
    if (atProblemFlag_) break;
    //TopIn.SetNonbond().SetLJ(idx).SetA( atof(file_.NextElement()) );
  }

  return 0;
}

// Parm_Amber::ReadLJB()
int Parm_Amber::ReadLJB(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_LJ_B, numLJparm_, FMT)) return 1;
  for (int idx = 0; idx != numLJparm_; idx++)
  {
    TopIn.SetNonbond().SetLJ(idx).SetB( FileBufferToDouble(F_LJ_B, idx, numLJparm_) );
    if (atProblemFlag_) break;
  }
  return 0;
}

// Parm_Amber::ReadLJC()
int Parm_Amber::ReadLJC(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_LJ_C, numLJparm_, FMT)) return 1;
  for (int idx = 0; idx != numLJparm_; idx++)
  {
    TopIn.SetNonbond().AddLJC( FileBufferToDouble(F_LJ_C, idx, numLJparm_) );
    if (atProblemFlag_) break;
  }
  return 0;
}

// Parm_Amber::GetBond()
/** Amber bond indices are * 3, bond parm indices are +1.
  * NOTE: Since NextElement() only returns a pointer to the file buffer,
  *       have to convert a number as soon as it is available.
  */
BondType Parm_Amber::GetBond() {
  int a1 = atoi(file_.NextElement());
  int a2 = atoi(file_.NextElement());
  int bidx = atoi(file_.NextElement());
  return BondType( a1 / 3, a2 / 3, bidx - 1 );
}

// Parm_Amber::ReadBondsH()
int Parm_Amber::ReadBondsH(Topology& TopIn, FortranData const& FMT) {
  int nvals = values_[NBONH]*3;
  if (SetupBuffer(F_BONDSH, nvals, FMT)) return 1;
  for (int idx = 0; idx != nvals; idx += 3)
    TopIn.AddBond( GetBond(), true );
  return 0;
}

// Parm_Amber::ReadBonds()
int Parm_Amber::ReadBonds(Topology& TopIn, FortranData const& FMT) {
  int nvals = values_[MBONA]*3;
  if (SetupBuffer(F_BONDS, nvals, FMT)) return 1;
  for (int idx = 0; idx != nvals; idx += 3)
    TopIn.AddBond( GetBond(), false );
  return 0;
}

// Parm_Amber::GetAngle()
AngleType Parm_Amber::GetAngle() {
  int a1 = atoi(file_.NextElement());
  int a2 = atoi(file_.NextElement());
  int a3 = atoi(file_.NextElement());
  int aidx = atoi(file_.NextElement());
  return AngleType( a1 / 3, a2 / 3, a3 / 3, aidx - 1 );
}

// Parm_Amber::ReadAnglesH()
int Parm_Amber::ReadAnglesH(Topology& TopIn, FortranData const& FMT) {
  int nvals = values_[NTHETH]*4;
  if (SetupBuffer(F_ANGLESH, nvals, FMT)) return 1;
  for (int idx = 0; idx != nvals; idx += 4)
    TopIn.AddAngle( GetAngle(), true );
  return 0;
}

// Parm_Amber::ReadAngles()
int Parm_Amber::ReadAngles(Topology& TopIn, FortranData const& FMT) {
  int nvals = values_[MTHETA]*4;
  if (SetupBuffer(F_ANGLES, nvals, FMT)) return 1;
  for (int idx = 0; idx != nvals; idx += 4)
    TopIn.AddAngle( GetAngle(), false );
  return 0;
}

// Parm_Amber::GetDihedral()
DihedralType Parm_Amber::GetDihedral() {
  int a1 = atoi(file_.NextElement());
  int a2 = atoi(file_.NextElement());
  int a3 = atoi(file_.NextElement());
  int a4 = atoi(file_.NextElement());
  int didx = atoi(file_.NextElement());
  return DihedralType( a1 / 3, a2 / 3, a3 / 3, a4 / 3, didx - 1 );
}

// Parm_Amber::ReadDihedralsH()
int Parm_Amber::ReadDihedralsH(Topology& TopIn, FortranData const& FMT) {
  int nvals = values_[NPHIH]*5;
  if (SetupBuffer(F_DIHH, nvals, FMT)) return 1;
  for (int idx = 0; idx != nvals; idx += 5)
    TopIn.AddDihedral( GetDihedral(), true );
  return 0;
}

// Parm_Amber::ReadDihedrals()
int Parm_Amber::ReadDihedrals(Topology& TopIn, FortranData const& FMT) {
  int nvals = values_[MPHIA]*5;
  if (SetupBuffer(F_DIH, nvals, FMT)) return 1;
  for (int idx = 0; idx != nvals; idx += 5)
    TopIn.AddDihedral( GetDihedral(), false );
  return 0;
}

// Parm_Amber::ReadAsol()
int Parm_Amber::ReadAsol(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_ASOL, values_[NPHB], FMT)) return 1;
  for (int idx = 0; idx != values_[NPHB]; idx++)
    TopIn.SetNonbond().SetHB(idx).SetAsol( atof(file_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadBsol()
int Parm_Amber::ReadBsol(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_BSOL, values_[NPHB], FMT)) return 1;
  for (int idx = 0; idx != values_[NPHB]; idx++)
    TopIn.SetNonbond().SetHB(idx).SetBsol( atof(file_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadHBcut()
int Parm_Amber::ReadHBcut(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_HBCUT, values_[NPHB], FMT)) return 1;
  for (int idx = 0; idx != values_[NPHB]; idx++)
    TopIn.SetNonbond().SetHB(idx).SetHBcut( atof(file_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadAtomTypes()
int Parm_Amber::ReadAtomTypes(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_TYPES, values_[NATOM], FMT)) return 1;
  for (int idx = 0; idx != values_[NATOM]; idx++)
    TopIn.SetAtom(idx).SetTypeName( NameType(file_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadItree()
int Parm_Amber::ReadItree(Topology& TopIn, FortranData const& FMT) {
  TopIn.AllocTreeChainClassification();
  if (SetupBuffer(F_ITREE, values_[NATOM], FMT)) return 1;
  for (int idx = 0; idx != values_[NATOM]; idx++)
    TopIn.SetTreeChainClassification(idx, file_.NextElement());
  return 0;
}

// Parm_Amber::ReadJoin()
int Parm_Amber::ReadJoin(Topology& TopIn, FortranData const& FMT) {
  TopIn.AllocJoinArray();
  if (SetupBuffer(F_JOIN, values_[NATOM], FMT)) return 1;
  for (int idx = 0; idx != values_[NATOM]; idx++)
    TopIn.SetJoinArray(idx, atoi(file_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadIrotat()
int Parm_Amber::ReadIrotat(Topology& TopIn, FortranData const& FMT) {
  TopIn.AllocRotateArray();
  if (SetupBuffer(F_IROTAT, values_[NATOM], FMT)) return 1;
  for (int idx = 0; idx != values_[NATOM]; idx++)
    TopIn.SetRotateArray(idx, atoi(file_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadBox()
int Parm_Amber::ReadBox(FortranData const& FMT) {
  if (SetupBuffer(F_PARMBOX, 4, FMT)) return 1;
  double xyzabg[6];
  xyzabg[Box::BETA] = atof(file_.NextElement());
  xyzabg[Box::X   ] = atof(file_.NextElement());
  xyzabg[Box::Y   ] = atof(file_.NextElement());
  xyzabg[Box::Z   ] = atof(file_.NextElement());
  //parmbox_.SetBetaLengths( beta, bx, by, bz );
  // Only beta angle is set (e.g. from Amber topology).
  if (xyzabg[Box::BETA] == 90.0) {
    xyzabg[Box::ALPHA] = 90.0;
    xyzabg[Box::GAMMA] = 90.0;
    if (debug_>0) mprintf("\tAmber topology box is orthogonal.\n");
  } else if ( Box::IsTruncOct( xyzabg[Box::BETA] ) ) {
    // Use trunc oct angle from Box; higher precision
    //xyzabg[Box::BETA ] = Box::TruncatedOctAngle();
    xyzabg[Box::ALPHA] = xyzabg[Box::BETA];
    xyzabg[Box::GAMMA] = xyzabg[Box::BETA];
    if (debug_>0) mprintf("\tAmber topology box is truncated octahedron.\n");
  } else if (xyzabg[Box::BETA] == 60.0) {
    xyzabg[Box::ALPHA] = 60.0; 
    xyzabg[Box::BETA]  = 90.0; 
    xyzabg[Box::GAMMA] = 60.0;
    if (debug_>0) mprintf("\tAmber topology box is rhombic dodecahedron, alpha=gamma=60.0, beta=90.0.\n");
  } else {
    mprintf("Warning: AmberParm: Unrecognized beta (%g); setting all angles to beta.\n", xyzabg[Box::BETA]);
    xyzabg[Box::ALPHA] = xyzabg[Box::BETA];
    xyzabg[Box::GAMMA] = xyzabg[Box::BETA];
  }
  parmbox_.SetupFromXyzAbg( xyzabg );
 
  return 0;
}

// Parm_Amber::ReadCapInfo()
int Parm_Amber::ReadCapInfo(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_CAP_INFO, 1, FMT)) return 1;
  TopIn.SetCap().SetNatcap( atoi(file_.NextElement())-1 );
  return 0;
}

// Parm_Amber::ReadCapInfo2()
int Parm_Amber::ReadCapInfo2(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_CAP_INFO2, 4, FMT)) return 1;
  TopIn.SetCap().SetCutCap( atof(file_.NextElement()) );
  TopIn.SetCap().SetXcap( atof(file_.NextElement()) );
  TopIn.SetCap().SetYcap( atof(file_.NextElement()) );
  TopIn.SetCap().SetZcap( atof(file_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadGBradiiSet()
int Parm_Amber::ReadGBradiiSet(Topology& TopIn) {
  std::string radius_set = NoTrailingWhitespace(file_.GetLine());
  mprintf("\tRadius Set: %s\n",radius_set.c_str());
  TopIn.SetGBradiiSet( radius_set );
  return 0;
}

// Parm_Amber::ReadGBradii()
int Parm_Amber::ReadGBradii(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_RADII, values_[NATOM], FMT)) return 1;
  for (int idx = 0; idx != values_[NATOM]; idx++)
    TopIn.SetAtom(idx).SetGBradius( atof(file_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadGBscreen()
int Parm_Amber::ReadGBscreen(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_SCREEN, values_[NATOM], FMT)) return 1;
  for (int idx = 0; idx != values_[NATOM]; idx++)
    TopIn.SetAtom(idx).SetGBscreen( atof(file_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadIpol()
int Parm_Amber::ReadIpol(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_IPOL, 1, FMT)) return 1;
  TopIn.SetIpol( atoi(file_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadPolar()
int Parm_Amber::ReadPolar(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_POLAR, values_[NATOM], FMT)) return 1;
  for (int idx = 0; idx != values_[NATOM]; idx++)
    TopIn.SetAtom(idx).SetPolar( atof(file_.NextElement()) );
  return 0;
}

// ----- Extra PDB Info --------------------------
int Parm_Amber::ReadPdbRes(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_PDB_RES, values_[NRES], FMT)) return 1;
  for (int idx = 0; idx != values_[NRES]; idx++)
    TopIn.SetRes(idx).SetOriginalNum( atoi(file_.NextElement()) );
  return 0;
}

int Parm_Amber::ReadPdbChainID(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_PDB_CHAIN, values_[NRES], FMT)) return 1;
  for (int idx = 0; idx != values_[NRES]; idx++)
    TopIn.SetRes(idx).SetChainID( *(file_.NextElement()) );
  return 0;
}

int Parm_Amber::ReadPdbIcode(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_PDB_ICODE, values_[NRES], FMT)) return 1;
  for (int idx = 0; idx != values_[NRES]; idx++)
    TopIn.SetRes(idx).SetIcode( *(file_.NextElement()) );
  return 0;
}

int Parm_Amber::ReadPdbAlt(Topology& TopIn, FortranData const& FMT) {
  TopIn.AllocAtomAltLoc();
  if (SetupBuffer(F_PDB_ALT, values_[NATOM], FMT)) return 1;
  for (int idx = 0; idx != values_[NATOM]; idx++)
    TopIn.SetAtomAltLoc(idx, *(file_.NextElement()) );
  return 0;
}

int Parm_Amber::ReadPdbBfactor(Topology& TopIn, FortranData const& FMT) {
  TopIn.AllocBfactor();
  if (SetupBuffer(F_PDB_BFAC, values_[NATOM], FMT)) return 1;
  for (int idx = 0; idx != values_[NATOM]; idx++)
    TopIn.SetBfactor(idx, atof(file_.NextElement()) );
  return 0;
}

int Parm_Amber::ReadPdbOccupancy(Topology& TopIn, FortranData const& FMT) {
  TopIn.AllocOccupancy();
  if (SetupBuffer(F_PDB_OCC, values_[NATOM], FMT)) return 1;
  for (int idx = 0; idx != values_[NATOM]; idx++)
    TopIn.SetOccupancy(idx, atof(file_.NextElement()) );
  return 0;
}

int Parm_Amber::ReadPdbNumbers(Topology& TopIn, FortranData const& FMT) {
  TopIn.AllocPdbSerialNum();
  if (SetupBuffer(F_PDB_NUM, values_[NATOM], FMT)) return 1;
  for (int idx = 0; idx != values_[NATOM]; idx++)
    TopIn.SetPdbSerialNum(idx, atoi(file_.NextElement()) );
  return 0;
}

// ----- CHAMBER ---------------------------------
static inline int fftype_err(const char* ptr, const char* Flag) {
  if (ptr == 0) {
    mprinterr("Error: Unexpected EOF when reading '%s'\n", Flag);
    return 1;
  }
  if (ptr[0] == '%') {
    mprintf("Warning: Section '%s' appears to have incorrect # lines.\n", Flag);
    return 2;
  }
  return 0;
}

// ReadChamberFFtype()
/** Format should be (nlines, text). Only the first nlines value is used
  * to determine how many lines of text to read.
  * Also set the number of LJ 1-4 terms.
  * FIXME Should the LJ setup be somewhere else?
  * NOTE: This fortran format is different than the others since it has
  *       two mixed types. Assume the format is IX,A80-X (currently X=2).
  *       If this fails print a warning and try to continue.
  */
int Parm_Amber::ReadChamberFFtype(Topology& TopIn, FortranData const& FMT) {
  mprintf("\tCHAMBER topology:\n");
  if (FMT.Ftype() != FINT)
    mprintf("Warning: In '%s' expected format to begin with integer. Skipping.\n",
            FLAGS_[F_FF_TYPE].Flag);
  else {
    const char* Flag = FLAGS_[F_FF_TYPE].Flag;
    const char* ptr = file_.NextLine(); // Read first line
    int err = fftype_err(ptr, Flag);
    if (err == 1) return 1;
    if (err == 0) {
      char* tmpbuf = new char[ FMT.Width()+1 ];
      tmpbuf[FMT.Width()] = '\0';
      std::copy(ptr, ptr+FMT.Width(), tmpbuf);
      int nlines = atoi(tmpbuf);
      delete[] tmpbuf;
      if (nlines > 0) {
        std::string ff_desc = NoTrailingWhitespace( ptr + FMT.Width() );
        mprintf("  %s\n", ff_desc.c_str());
        TopIn.SetChamber().AddDescription( ff_desc );
        for (int line = 1; line < nlines; line++) {
          ptr = file_.NextLine();
          err = fftype_err(ptr, Flag);
          if (err == 1)
            return 1;
          else if (err == 2) break;
          ff_desc = NoTrailingWhitespace( ptr + FMT.Width() );
          mprintf("  %s\n", ff_desc.c_str());
          TopIn.SetChamber().AddDescription( ff_desc );
        }
      }
    }
  }
  TopIn.SetChamber().SetNLJ14terms( numLJparm_ );
  return 0;
}

// Parm_Amber::ReadChamberUBCount()
int Parm_Amber::ReadChamberUBCount(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_CHM_UBC, 2, FMT)) return 1;
  UB_count_[0] = atoi(file_.NextElement()); // Number of bonds
  UB_count_[1] = atoi(file_.NextElement()); // Number of parameters
  TopIn.SetChamber().ReserveUBterms( UB_count_[0] );
  TopIn.SetChamber().ResizeUBparm( UB_count_[1] );
  UB_count_[0] *= 3; // Number of bond terms (bonds x 3)
  return 0;
}

// Parm_Amber::ReadChamberUBTerms()
int Parm_Amber::ReadChamberUBTerms(Topology& TopIn, FortranData const& FMT) {
  // NOTE: UB_count_[0] already multiplied by 3
  if (SetupBuffer(F_CHM_UB, UB_count_[0], FMT)) return 1;
  for (int idx = 0; idx != UB_count_[0]; idx += 3) {
    int a1 = atoi(file_.NextElement()) - 1;
    int a2 = atoi(file_.NextElement()) - 1;
    int bidx = atoi(file_.NextElement()) - 1;
    TopIn.SetChamber().AddUBterm( BondType(a1, a2, bidx) );
  }
  return 0;
}

// Parm_Amber::ReadChamberUBFC()
int Parm_Amber::ReadChamberUBFC(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_CHM_UBFC, UB_count_[1], FMT)) return 1;
  for (int idx = 0; idx != UB_count_[1]; idx++)
    TopIn.SetChamber().SetUBparm(idx).SetRk( atof(file_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadChamberUBEQ()
int Parm_Amber::ReadChamberUBEQ(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_CHM_UBEQ, UB_count_[1], FMT)) return 1;
  for (int idx = 0; idx != UB_count_[1]; idx++)
    TopIn.SetChamber().SetUBparm(idx).SetReq( atof(file_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadChamberNumImpropers()
int Parm_Amber::ReadChamberNumImpropers(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_CHM_NIMP, 1, FMT)) return 1;
  N_impropers_ = atoi(file_.NextElement());
  TopIn.SetChamber().ReserveImproperTerms( N_impropers_ );
  N_impropers_ *= 5; // Number of improper terms (impropers * 5)
  return 0;
}

// Parm_Amber::ReadChamberNumImpTerms()
int Parm_Amber::ReadChamberNumImpTerms(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_CHM_NIMPT, 1, FMT)) return 1;
  N_impTerms_ = atoi(file_.NextElement());
  TopIn.SetChamber().ResizeImproperParm( N_impTerms_ );
  return 0;
}

// Parm_Amber::ReadChamberImpropers()
int Parm_Amber::ReadChamberImpropers(Topology& TopIn, FortranData const& FMT) {
  // NOTE: N_impropers_ already multiplied by 5
  if (SetupBuffer(F_CHM_IMP, N_impropers_, FMT)) return 1;
  for (int idx = 0; idx != N_impropers_; idx +=5) {
    int a1 = atoi(file_.NextElement()) - 1;
    int a2 = atoi(file_.NextElement()) - 1;
    int a3 = atoi(file_.NextElement()) - 1;
    int a4 = atoi(file_.NextElement()) - 1;
    int didx = atoi(file_.NextElement()) - 1;
    TopIn.SetChamber().AddImproperTerm( DihedralType(a1, a2, a3, a4, didx) );
  }
  return 0;
}

// Parm_Amber::ReadChamberImpFC()
int Parm_Amber::ReadChamberImpFC(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_CHM_IMPFC, N_impTerms_, FMT)) return 1;
  for (int idx = 0; idx != N_impTerms_; idx++)
    TopIn.SetChamber().SetImproperParm(idx).SetPk( atof(file_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadChamberImpPHASE()
int Parm_Amber::ReadChamberImpPHASE(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_CHM_IMPP, N_impTerms_, FMT)) return 1;
  for (int idx = 0; idx != N_impTerms_; idx++)
    TopIn.SetChamber().SetImproperParm(idx).SetPhase( atof(file_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadChamberLJ14A()
int Parm_Amber::ReadChamberLJ14A(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_LJ14A, numLJparm_, FMT)) return 1;
  for (int idx = 0; idx != numLJparm_; idx++) {
    TopIn.SetChamber().SetLJ14(idx).SetA( FileBufferToDouble(F_LJ14A, idx, numLJparm_) );
    if (atProblemFlag_) break;
  }
  return 0;
}

// Parm_Amber::ReadChamberLJ14B()
int Parm_Amber::ReadChamberLJ14B(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_LJ14B, numLJparm_, FMT)) return 1;
  for (int idx = 0; idx != numLJparm_; idx++) {
    TopIn.SetChamber().SetLJ14(idx).SetB( FileBufferToDouble(F_LJ14B, idx, numLJparm_) );
    if (atProblemFlag_) break;
  }
  return 0;
}

// Parm_Amber::ReadChamberCmapCounts()
int Parm_Amber::ReadCmapCounts(FlagType ftype, FortranData const& FMT) {
  if (SetupBuffer(ftype, 2, FMT)) return 1;
  n_cmap_terms_ = atoi( file_.NextElement() );
  n_cmap_grids_ = atoi( file_.NextElement() );
  return 0;
}

// Parm_Amber::ReadChamberCmapRes()
/** Get CMAP resolutions for each grid and allocate grids. */
int Parm_Amber::ReadCmapRes(FlagType ftype, Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(ftype, n_cmap_grids_, FMT)) return 1;
  for (int i = 0; i != n_cmap_grids_; i++)
    TopIn.AddCmapGrid( CmapGridType( atoi(file_.NextElement()) ) );
  return 0;
}

// Parm_Amber::ReadChamberCmapGrid()
/** Read CMAP grid. It is expected that the grid number is between 1 and CHARMM_CMAP_COUNT */
int Parm_Amber::ReadCmapGrid(FlagType ftype, const char* CmapFlag, Topology& TopIn, FortranData const& FMT)
{
  if (CmapFlag == 0) {
    mprinterr("Internal Error: ReadCmapGrid: CmapFlag is null.\n");
    return 1;
  }
  // Figure out which grid this is.
  //           11111111112222
  // 012345678901234567890123
  // CHARMM_CMAP_PARAMETER_XX or CMAP_PARAMETER_XX
  // The former is what is specified by the Amber parm/top file format spec.
  // The latter appears to be what is output from Charmm GUI
  // Advance to the final underscore.
  const char* ptr = CmapFlag;
  unsigned int underscore_idx = 0;
  while (*ptr != '\0') {
    if (*ptr == '_') underscore_idx = (ptr - CmapFlag);
    ptr++;
  }
  // Get string
  std::string cmap_index_str( CmapFlag+underscore_idx+1 );
  // Sanity check
  if (cmap_index_str.size() != 2 || !validInteger(cmap_index_str)) {
    mprinterr("Error: CMAP index flag %s does not appear to contain an integer: %s\n",
             CmapFlag, cmap_index_str.c_str());
    return 1;
  }
  //mprintf("DEBUG: cmap index string: %s\n", cmap_index_str.c_str());
  int gridnum = convertToInteger( cmap_index_str ) - 1;
  if (gridnum < 0 || gridnum >= (int)TopIn.CmapGrid().size()) {
    mprintf("Warning: CMAP grid '%s' out of range.\n", CmapFlag);
    if (TopIn.HasCmap())
      mprintf("Warning: Expected grid between 1 and %zu, got %i\n",
              TopIn.CmapGrid().size(), gridnum+1);
    else
      mprintf("Warning: Missing previous %s section.\n", FLAGS_[ftype].Flag);
    mprintf("Warning: Skipping read of CMAP grid.\n");
    return 0;
  }
  CmapGridType& GRID = TopIn.SetCmapGrid( gridnum );
  if (SetupBuffer(ftype, GRID.Size(), FMT)) return 1;
  for (int idx = 0; idx != GRID.Size(); idx++)
    GRID.SetGridPt( idx, atof(file_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadChamberCmapTerms()
int Parm_Amber::ReadCmapTerms(FlagType ftype, Topology& TopIn, FortranData const& FMT) {
  int nvals = n_cmap_terms_ * 6;
  if (SetupBuffer(ftype, nvals, FMT)) return 1;
  for (int idx = 0; idx != nvals; idx += 6) {
    int a1 = atoi(file_.NextElement()) - 1;
    int a2 = atoi(file_.NextElement()) - 1;
    int a3 = atoi(file_.NextElement()) - 1;
    int a4 = atoi(file_.NextElement()) - 1;
    int a5 = atoi(file_.NextElement()) - 1;
    int gidx = atoi(file_.NextElement()) - 1;
    TopIn.AddCmapTerm( CmapType(a1, a2, a3, a4, a5, gidx) );
  }
  return 0;
}

// ----- LES -------------------------------------
int Parm_Amber::ReadLESntyp(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_LES_NTYP, 1, FMT)) return 1;
  nlestyp_ = atoi(file_.NextElement());
  // This allocates the FAC array (nlestyp*nlestyp) and reserves LES array (natoms)
  TopIn.SetLES().Allocate( values_[NATOM], nlestyp_ );
  return 0;
}

int Parm_Amber::ReadLESfac(Topology& TopIn, FortranData const& FMT) {
  int nvals = nlestyp_ * nlestyp_;
  if (SetupBuffer(F_LES_FAC, nvals, FMT)) return 1;
  for (int idx = 0; idx != nvals; idx++)
    TopIn.SetLES().SetFAC( idx, atof(file_.NextElement()) );
  return 0;
}

int Parm_Amber::ReadLEStypes(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_LES_TYPE, values_[NATOM], FMT)) return 1;
  for (int idx = 0; idx != values_[NATOM]; idx++)
    TopIn.SetLES().SetType( idx, atoi(file_.NextElement()) );
  return 0;
}

int Parm_Amber::ReadLEScnum(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_LES_CNUM, values_[NATOM], FMT)) return 1;
  for (int idx = 0; idx != values_[NATOM]; idx++)
    TopIn.SetLES().SetCopy( idx, atoi(file_.NextElement()) );
  return 0;
}

int Parm_Amber::ReadLESid(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_LES_ID, values_[NATOM], FMT)) return 1;
  for (int idx = 0; idx != values_[NATOM]; idx++)
    TopIn.SetLES().SetID( idx, atoi(file_.NextElement()) );
  return 0;
}

// -----------------------------------------------------------------------------
void Parm_Amber::WriteHelp() {
  mprintf("\tnochamber  : Do not write CHAMBER information to topology (useful for e.g. using"
          "\t             topology for visualization with VMD).\n");
  mprintf("\twriteempty : Write Amber tree, join, and rotate info even if not present.\n");
  mprintf("\tnopdbinfo  : Do not write \"PDB\" info (e.g. chain IDs, original res #s, etc).\n");
}

int Parm_Amber::processWriteArgs(ArgList& argIn) {
  writeChamber_ = !argIn.hasKey("nochamber");
  writeEmptyArrays_ = argIn.hasKey("writeempty");
  writePdbInfo_ = !argIn.hasKey("nopdbinfo");
  return 0;
}

/** \return Format for given flag. */
Parm_Amber::FortranData Parm_Amber::WriteFormat(FlagType fflag) const {
  FortranData FMT;
  // For chamber, certain flags have different format (boo).
  if (ptype_ == CHAMBER) {
    if      (fflag == F_CHARGE  ) FMT.ParseFortranFormat("%FORMAT(3E24.16)");
    else if (fflag == F_ANGLETEQ) FMT.ParseFortranFormat("%FORMAT(3E25.17)");
    else if (fflag == F_LJ_A    ) FMT.ParseFortranFormat("%FORMAT(3E24.16)");
    else if (fflag == F_LJ_B    ) FMT.ParseFortranFormat("%FORMAT(3E24.16)");
  }
  if (FMT.Ftype() == UNKNOWN_FTYPE)
    FMT.ParseFortranFormat( FLAGS_[fflag].Fmt );
  //mprintf("DEBUG: Flag '%s' format '%s'\n", FLAGS_[fflag].Flag, FMT.Fstr());
  return FMT;
}

/** This version uses the default flag to allocate the format string. */
int Parm_Amber::BufferAlloc(FlagType ftype, int nvals, int idx) {
  FortranData FMT = WriteFormat( ftype );
  return BufferAlloc(ftype, FMT, nvals, idx);
}

// Parm_Amber::BufferAlloc()
/** Allocate internal buffer for writing given flag with given format.
  * \param ftype The FLAG to write.
  * \param FMT The Fortran format string.
  * \param nvals The number of values that will be written.
  * \param idx Index for flags that can be specified multiple times (e.g. CMAP grids).
  */
int Parm_Amber::BufferAlloc(FlagType ftype, FortranData const& FMT, int nvals, int idx) {
  if ( FMT.Ftype() == UNKNOWN_FTYPE) {
    mprinterr("Interal Error: Could not set up format string.\n");
    return 1;
  }
  // Write FLAG and FORMAT lines
  if (idx < 0)
    file_.Printf("%%FLAG %-74s\n%-80s\n", FLAGS_[ftype].Flag, FMT.Fstr());
  else {
    // NOTE: Currently only needed for CMAP grid flags
    std::string fflag( FLAGS_[ftype].Flag );
    fflag.append( integerToString( idx, 2 ) );
    file_.Printf("%%FLAG %-74s\n%-80s\n", fflag.c_str(), FMT.Fstr());
  }
  if (nvals > 0) {
    TextFormat WriteFmt;
    //mprintf("DEBUG: Set up write buffer for '%s', %i vals.\n", FLAGS_[ftype].Flag, nvals);
    if      (FMT.Ftype() == FINT)
      WriteFmt = TextFormat(TextFormat::INTEGER, FMT.Width());
    else if (FMT.Ftype() == FDOUBLE)
      WriteFmt = TextFormat(TextFormat::SCIENTIFIC, FMT.Width(), FMT.Precision());
    else if (FMT.Ftype() == FCHAR)
      WriteFmt = TextFormat(TextFormat::STRING, FMT.Width(), TextFormat::LEFT);
    else if (FMT.Ftype() == FFLOAT)
      WriteFmt = TextFormat(TextFormat::DOUBLE, FMT.Width(), FMT.Precision());
    //mprintf("DEBUG: Write format: \"%s\"\n", WriteFmt.fmt());
    file_.SetupFrameBuffer( nvals, WriteFmt, FMT.Ncols() );
  } else {
    //mprintf("DEBUG: No values for flag '%s'\n", FLAGS_[ftype].Flag);
    // Write blank line
    file_.Printf("\n");
  }
  return 0;
}

// Parm_Amber::WriteBondParm()
int Parm_Amber::WriteBondParm(FlagType RKflag, FlagType REQflag, BondParmArray const& BP) {
  // BOND RK
  if (BufferAlloc(RKflag, BP.size())) return 1;
  for (BondParmArray::const_iterator it = BP.begin(); it != BP.end(); ++it)
    file_.DblToBuffer( it->Rk() );
  file_.FlushBuffer();
  // BOND REQ
  if (BufferAlloc(REQflag, BP.size())) return 1;
  for (BondParmArray::const_iterator it = BP.begin(); it != BP.end(); ++it)
    file_.DblToBuffer( it->Req() );
  file_.FlushBuffer();
  return 0;
}

// Parm_Amber::WriteLJ()
int Parm_Amber::WriteLJ(FlagType Aflag, FlagType Bflag, NonbondArray const& NB) {
  // LJ A terms
  if (BufferAlloc(Aflag, NB.size())) return 1;
  for (NonbondArray::const_iterator it = NB.begin(); it != NB.end(); ++it)
    file_.DblToBuffer( it->A() );
  file_.FlushBuffer();
  // LJ B terms
  if (BufferAlloc(Bflag, NB.size())) return 1;
  for (NonbondArray::const_iterator it = NB.begin(); it != NB.end(); ++it)
    file_.DblToBuffer( it->B() );
  file_.FlushBuffer();
  return 0;
}

/** Amber bond array. Indices must be x3, parameter index +1 */
int Parm_Amber::WriteBonds(FlagType flag, BondArray const& BND) {
  if (BufferAlloc(flag, BND.size()*3)) return 1;
  for (BondArray::const_iterator it = BND.begin(); it != BND.end(); ++it) {
    file_.IntToBuffer( it->A1()*3 );
    file_.IntToBuffer( it->A2()*3 );
    file_.IntToBuffer( it->Idx()+1 );
  }
  file_.FlushBuffer();
  return 0;
}

/** Amber angle array. Indices must be x3, parameter index +1 */
int Parm_Amber::WriteAngles(FlagType flag, AngleArray const& ANG) {
  if (BufferAlloc(flag, ANG.size()*4)) return 1;
  for (AngleArray::const_iterator it = ANG.begin(); it != ANG.end(); ++it) {
    file_.IntToBuffer( it->A1()*3 );
    file_.IntToBuffer( it->A2()*3 );
    file_.IntToBuffer( it->A3()*3 );
    file_.IntToBuffer( it->Idx()+1 );
  }
  file_.FlushBuffer();
  return 0;
}

/** Amber dihedral array. Indices must be x3, parameters index +1. End 
  * dihedrals have the third atom index negative, impropers have fourth.
  */
int Parm_Amber::WriteDihedrals(FlagType flag, DihedralArray const& DIH) {
  if (BufferAlloc(flag, DIH.size()*5)) return 1;
  for (DihedralArray::const_iterator it = DIH.begin(); it != DIH.end(); ++it) {
    file_.IntToBuffer( it->A1()*3 );
    file_.IntToBuffer( it->A2()*3 );
    if ( it->Type() == DihedralType::BOTH || it->Type() == DihedralType::END)
      file_.IntToBuffer( -(it->A3()*3) );
    else
      file_.IntToBuffer( it->A3()*3 );
    if ( it->Type() == DihedralType::BOTH || it->Type() == DihedralType::IMPROPER)
      file_.IntToBuffer( -(it->A4()*3) );
    else
      file_.IntToBuffer( it->A4()*3 );
    file_.IntToBuffer( it->Idx()+1 );
  }
  file_.FlushBuffer();
  return 0;
}

/** Write a single line to the topology file. */
void Parm_Amber::WriteLine(FlagType flag, std::string const& lineIn) {
  std::string title = lineIn;
  // Resize to max 80 char
  if (title.size() > 80)
    title.resize(80);
  file_.Printf("%%FLAG %-74s\n%-80s\n%-80s\n", FLAGS_[flag].Flag, FLAGS_[flag].Fmt, title.c_str());
}

/** Write TREE_CHAIN_CLASSIFICATION array. */
int Parm_Amber::WriteTreeChainClassification(std::vector<NameType> const& tree) {
  if (BufferAlloc(F_ITREE, tree.size())) return 1;
  for (std::vector<NameType>::const_iterator it = tree.begin(); it != tree.end(); ++it)
    file_.CharToBuffer( *(*it) );
  file_.FlushBuffer();
  return 0;
}

/** Write IJOIN array */
int Parm_Amber::WriteIjoin(std::vector<int> const& ijoin) {
  if (BufferAlloc(F_JOIN, ijoin.size())) return 1;
  for (std::vector<int>::const_iterator it = ijoin.begin(); it != ijoin.end(); ++it)
    file_.IntToBuffer( *it );
  file_.FlushBuffer();
  return 0;
}

/** Write IROTAT array */
int Parm_Amber::WriteIrotat(std::vector<int> const& irotat) {
  if (BufferAlloc(F_IROTAT, irotat.size())) return 1;
  for (std::vector<int>::const_iterator it = irotat.begin(); it != irotat.end(); ++it)
    file_.IntToBuffer( *it );
  file_.FlushBuffer();
  return 0;
}

// Parm_Amber::WriteExtra()
/** Write "extra" Amber info that is deprecated: the tree, join, and rotate
  * arrays. If the info is not present but it is indicated the user
  * wants them written, generate blank arrays of size natom.
  */
int Parm_Amber::WriteExtra(Topology const& topOut, int natom) {
  if (topOut.TreeChainClassification().empty() ||
      topOut.JoinArray().empty() ||
      topOut.RotateArray().empty())
  {
    mprintf("Warning: Topology does not contain tree, join, and/or rotate arrays.\n");
    if (writeEmptyArrays_)
      mprintf("Warning: Empty entries will be created for missing arrays. Care should be\n"
              "Warning:   taken if this topology is used for anything besides analysis or\n"
              "Warning:   visualization.\n");
    else
      mprintf("Warning: Missing arrays will not be written. Care should be taken if this\n"
              "Warning:   topology is used for anything besides analysis or visualization.\n"
              "Warning: To create an empty array specify the 'writeempty' keyword.\n");
  }
  // TREE CHAIN CLASSIFICATION
  if (topOut.TreeChainClassification().empty()) {
    if (writeEmptyArrays_) {
      if (WriteTreeChainClassification( std::vector<NameType>(natom, "BLA") )) return 1;
    }
  } else {
    if (WriteTreeChainClassification( topOut.TreeChainClassification() )) return 1;
  }
  // JOIN
  if (topOut.JoinArray().empty()) {
    if (writeEmptyArrays_) {
      if (WriteIjoin( std::vector<int>(natom, 0) )) return 1;
    }
  } else {
    if (WriteIjoin( topOut.JoinArray() )) return 1;
  } 
  // IROTAT
  if (topOut.RotateArray().empty()) {
    if (writeEmptyArrays_) {
      if (WriteIrotat( std::vector<int>(natom, 0) )) return 1;
    }
  } else {
    if (WriteIrotat( topOut.RotateArray() )) return 1;
  }
  return 0;
}

// Parm_Amber::WriteParm()
int Parm_Amber::WriteParm(FileName const& fname, Topology const& TopOut) {
  if (file_.OpenWrite( fname )) return 1;
  elec_to_parm_ = Constants::ELECTOAMBER;
  parm_to_elec_ = Constants::AMBERTOELEC;
  // Determine if this is a CHAMBER topology
  ptype_ = NEWPARM;
  FlagType titleFlag = F_TITLE;
  if (TopOut.Chamber().HasChamber()) {
    if (!writeChamber_)
      mprintf("\tnochamber: Removing CHAMBER info from topology.\n");
    else {
      titleFlag = F_CTITLE;
      ptype_ = CHAMBER;
      elec_to_parm_ = ELECTOCHAMBER_;
      parm_to_elec_ = CHAMBERTOELEC_;
    }
  }
  // Warn about empty parameters
  if (TopOut.Nbonds()>0 && TopOut.BondParm().empty())
    mprintf("Warning: No bond parameters.\n");
  if (TopOut.Nangles()>0 && TopOut.AngleParm().empty())
    mprintf("Warning: No angle parameters.\n");
  if (TopOut.Ndihedrals()>0 && TopOut.DihedralParm().empty())
    mprintf("Warning: No dihedral parameters.\n");
  if (!TopOut.Nonbond().HasNonbond())
    mprintf("Warning: No non-bonded parameters.\n");
  // HEADER AND TITLE (4 lines, version, flag, format, title)
  file_.Printf("%-44s%s                  \n",
               "%VERSION  VERSION_STAMP = V0001.000  DATE = ",
               TimeString().c_str());
  WriteLine( titleFlag, TopOut.ParmName() );

  // Generate atom exclusion list. Do this here since POINTERS needs the size.
  ExclusionArray exclusionArray;
  if (exclusionArray.SetupExcluded(TopOut.Atoms(), 4,
                                   ExclusionArray::NO_EXCLUDE_SELF,
                                   ExclusionArray::ONLY_GREATER_IDX))
  {
    mprinterr("Error: Parm_Amber: Could not set up exclusion array for topology write.\n");
    return 1;
  }
  Iarray Excluded;
  for (ExclusionArray::const_iterator exList = exclusionArray.begin();
                                      exList != exclusionArray.end();
                                    ++exList)
  {
    if (exList->empty())
      Excluded.push_back( 0 );
    else {
      for (ExclusionArray::ExListType::const_iterator ex = exList->begin();
                                                      ex != exList->end(); ex++)
        // Amber atom #s start from 1
        Excluded.push_back( (*ex) + 1 );
    }
  }

  // Determine if atoms have bfactors and/or occupancy.
  bool hasBfac = !TopOut.Bfactor().empty();
  bool hasOcc  = !TopOut.Occupancy().empty();
  bool hasNum  = !TopOut.PdbSerialNum().empty();
  if (hasBfac)
    mprintf("\tTopology has atomic B-factors.\n");
  if (hasOcc)
    mprintf("\tTopology has atomic occupancies.\n");
  if (hasNum)
    mprintf("\tTopology has original PDB serial numbers.\n");
  // Determine max residue size. Also determine if extra info needs to be
  // written such as PDB residue number, chain ID, etc.
  // Since original res num gets set no matter what, only print out
  // if chain IDs are present or original res nums != res index + 1.
  int maxResSize = 0;
  bool hasOrigResNums = false;
  bool hasChainID = false;
  bool hasIcodes = false;
  for (Topology::res_iterator res = TopOut.ResStart(); res != TopOut.ResEnd(); ++res)
  {
    long int oresidx = res - TopOut.ResStart() + 1;
    if (oresidx != (long int)res->OriginalResNum())
      hasOrigResNums = true;
    if (res->HasChainID())
      hasChainID = true;
    if (res->Icode() != ' ')
      hasIcodes = true;
    maxResSize = std::max(maxResSize, res->NumAtoms());
  }
  if (hasOrigResNums)
    mprintf("\tTopology has alternative residue numbering (from e.g PDB, stripping, reordering, etc).\n");
  if (hasChainID)
    mprintf("\tTopology has chain IDs.\n");
  if (hasIcodes)
    mprintf("\tTopology has residue insertion codes.\n");
  // Turn PDB info off if specified
  if (!writePdbInfo_) {
    if (hasBfac)
      mprintf("\tnopdbinfo : Not writing atomic B-factors.\n");
    if (hasOcc)
      mprintf("\tnopdbinfo : Not writing atomic occupancies.\n");
    if (hasNum)
      mprintf("\tnopdbinfo : Not writing original PDB serial numbers.\n");
    if (hasOrigResNums)
      mprintf("\tnopdbinfo : Not writing alternative residue numbering.\n");
    if (hasChainID)
      mprintf("\tnopdbinfo : Not writing chain IDs.\n");
    if (hasIcodes)
      mprintf("\tnopdbinfo : Not writing residue insertion codes.\n");
    hasOrigResNums = false;
    hasChainID = false;
    hasIcodes = false;
    hasBfac = false;
    hasOcc = false;
    hasNum = false;
  }

  // Determine value of ifbox
  int ifbox = 0;
  Box::CellShapeType cellShape = TopOut.ParmBox().CellShape();
  if ( cellShape != Box::NO_SHAPE ) {
    if (cellShape == Box::CUBIC || cellShape == Box::TETRAGONAL || cellShape == Box::ORTHORHOMBIC) {
      // Orthorhombic
      if (!TopOut.ParmBox().Is_X_Aligned_Ortho())
        mprintf("Warning: Parm box shape is orthorhombic but cell is not X-aligned.\n");
      ifbox = 1;
    } else if (cellShape == Box::OCTAHEDRAL) {
      // Truncated octahedron
      ifbox = 2;
    } else {
      // General triclinic
      ifbox = 3;
    }
  }

  // POINTERS
  if (BufferAlloc(F_POINTERS, AMBERPOINTERS_)) return 1;
  file_.IntToBuffer( TopOut.Natom() ); // NATOM
  file_.IntToBuffer( TopOut.Nonbond().Ntypes() ); // NTYPES
  file_.IntToBuffer( TopOut.BondsH().size() ); // NBONH
  file_.IntToBuffer( TopOut.Bonds().size() ); // NBONA
  file_.IntToBuffer( TopOut.AnglesH().size() ); // NTHETH
  file_.IntToBuffer( TopOut.Angles().size() ); // NTHETA
  file_.IntToBuffer( TopOut.DihedralsH().size() ); // NPHIH
  file_.IntToBuffer( TopOut.Dihedrals().size() ); // NPHIA
  file_.IntToBuffer( 0 ); // NHPARM, not used
  // FIXME: Currently LES info not 100% correct, in particular the excluded list
  if (TopOut.LES().HasLES()) { // NPARM
    mprintf("Warning: Excluded atom list for LES info is not correct.\n");
    file_.IntToBuffer( 1 );
  } else
    file_.IntToBuffer( 0 );
  file_.IntToBuffer( Excluded.size() ); // NNB
  file_.IntToBuffer( TopOut.Nres() ); // NRES
  //   NOTE: Assuming MBONA == NBONA etc
  file_.IntToBuffer( TopOut.Bonds().size() ); // MBONA
  file_.IntToBuffer( TopOut.Angles().size() ); // MTHETA
  file_.IntToBuffer( TopOut.Dihedrals().size() ); // MPHIA
  file_.IntToBuffer( TopOut.BondParm().size() ); // NUMBND
  file_.IntToBuffer( TopOut.AngleParm().size() ); // NUMANG
  file_.IntToBuffer( TopOut.DihedralParm().size() ); // NPTRA
  file_.IntToBuffer( TopOut.NatomTypes() ); // NATYP, only for SOLTY
  file_.IntToBuffer( TopOut.Nonbond().HBarray().size() ); // NPHB
  file_.IntToBuffer( 0 ); // IFPERT
  file_.IntToBuffer( 0 ); // NBPER
  file_.IntToBuffer( 0 ); // NGPER
  file_.IntToBuffer( 0 ); // NDPER
  file_.IntToBuffer( 0 ); // MBPER
  file_.IntToBuffer( 0 ); // MGPER
  file_.IntToBuffer( 0 ); // MDPER
  file_.IntToBuffer( ifbox ); // IFBOX
  file_.IntToBuffer( maxResSize ); // NMXRS
  if (TopOut.Cap().NatCap() > 0) // IFCAP
    file_.IntToBuffer( 1 );
  else
    file_.IntToBuffer( 0 );
  file_.IntToBuffer( TopOut.NextraPts() ); // NEXTRA
  file_.FlushBuffer();
 
  // CHAMBER only - FF type description 
  if (ptype_ == CHAMBER) {
    file_.Printf("%%FLAG %-74s\n%-80s\n", FLAGS_[F_FF_TYPE].Flag, FLAGS_[F_FF_TYPE].Fmt);
    int nlines = (int)TopOut.Chamber().Description().size();
    if (nlines > 99) {
      mprintf("Warning: Number of CHAMBER description lines > 99. Only writing 99.\n", nlines);
      nlines = 99;
    }
    if (nlines > 0) {
      for (int line = 0; line != nlines; line++)
        file_.Printf("%2i%78s\n", nlines, TopOut.Chamber().Description()[line].c_str());
    } else {
      // No description. Write a placeholder.
      file_.Printf("%2i%78s\n", 1, "CHARMM:");
    }
  }

  // NAMES
  if (BufferAlloc(F_NAMES, TopOut.Natom())) return 1;
  for (Topology::atom_iterator atm = TopOut.begin(); atm != TopOut.end(); ++atm)
    file_.CharToBuffer( atm->c_str() );
  file_.FlushBuffer();

  // CHARGES
  if (BufferAlloc(F_CHARGE, TopOut.Natom())) return 1;
  for (Topology::atom_iterator atm = TopOut.begin(); atm != TopOut.end(); ++atm)
    file_.DblToBuffer( atm->Charge() * elec_to_parm_ );
  file_.FlushBuffer();

  // ATOMIC NUMBER
  if (BufferAlloc(F_ATOMICNUM, TopOut.Natom())) return 1;
  for (Topology::atom_iterator atm = TopOut.begin(); atm != TopOut.end(); ++atm)
    file_.IntToBuffer( atm->AtomicNumber() );
  file_.FlushBuffer();

  // MASS
  if (BufferAlloc(F_MASS, TopOut.Natom())) return 1;
  for (Topology::atom_iterator atm = TopOut.begin(); atm != TopOut.end(); ++atm)
    file_.DblToBuffer( atm->Mass() );
  file_.FlushBuffer();

  // TYPE INDEX
  if (BufferAlloc(F_ATYPEIDX, TopOut.Natom())) return 1;
  for (Topology::atom_iterator atm = TopOut.begin(); atm != TopOut.end(); ++atm)
    file_.IntToBuffer( atm->TypeIndex()+1 );
  file_.FlushBuffer();

  // NUMEX
  if (BufferAlloc(F_NUMEX, TopOut.Natom())) return 1;
  for (ExclusionArray::const_iterator exList = exclusionArray.begin();
                                      exList != exclusionArray.end(); ++exList)
  {
    if (exList->empty())
      file_.IntToBuffer( 1 );
    else
      file_.IntToBuffer( (int)exList->size() );
  }
  file_.FlushBuffer();

  // NONBONDED INDICES - positive needs to be shifted by +1 for fortran
  if (BufferAlloc(F_NB_INDEX, TopOut.Nonbond().NBindex().size())) return 1;
  for (Iarray::const_iterator it = TopOut.Nonbond().NBindex().begin();
                              it != TopOut.Nonbond().NBindex().end(); ++it)
    if (*it > -1)
      file_.IntToBuffer( *it + 1 );
    else
      file_.IntToBuffer( *it );
  file_.FlushBuffer();

  // RESIDUE NAME
  if (BufferAlloc(F_RESNAMES, TopOut.Nres())) return 1;
  for (Topology::res_iterator res = TopOut.ResStart(); res != TopOut.ResEnd(); ++res)
    file_.CharToBuffer( res->c_str() );
  file_.FlushBuffer();

  // RESIDUE POINTERS
  if (BufferAlloc(F_RESNUMS, TopOut.Nres())) return 1;
  for (Topology::res_iterator res = TopOut.ResStart(); res != TopOut.ResEnd(); ++res)
    file_.IntToBuffer( res->FirstAtom()+1 );
  file_.FlushBuffer();

  // BOND RK and REQ
  if (WriteBondParm(F_BONDRK, F_BONDREQ, TopOut.BondParm())) return 1;

  // ANGLE TK
  if (BufferAlloc(F_ANGLETK, TopOut.AngleParm().size())) return 1;
  for (AngleParmArray::const_iterator it = TopOut.AngleParm().begin();
                                      it != TopOut.AngleParm().end(); ++it)
    file_.DblToBuffer( it->Tk() );
  file_.FlushBuffer();

  // ANGLE TEQ
  if (BufferAlloc(F_ANGLETEQ, TopOut.AngleParm().size())) return 1;
  for (AngleParmArray::const_iterator it = TopOut.AngleParm().begin();
                                      it != TopOut.AngleParm().end(); ++it)
    file_.DblToBuffer( it->Teq() );
  file_.FlushBuffer();

  // CHARMM only - Urey-Bradley
  if (ptype_ == CHAMBER) {
    // UB COUNT
    if (BufferAlloc(F_CHM_UBC, 2)) return 1;
    file_.IntToBuffer( TopOut.Chamber().UB().size() );
    file_.IntToBuffer( TopOut.Chamber().UBparm().size() );
    file_.FlushBuffer();
    // UB terms
    if (BufferAlloc(F_CHM_UB, TopOut.Chamber().UB().size()*3)) return 1;
    for (BondArray::const_iterator it = TopOut.Chamber().UB().begin();
                                   it != TopOut.Chamber().UB().end(); ++it)
    {
      file_.IntToBuffer( it->A1()+1 );
      file_.IntToBuffer( it->A2()+1 );
      file_.IntToBuffer( it->Idx()+1 );
    }
    file_.FlushBuffer();
    // UB FORCE CONSTANTS and EQ
    if (WriteBondParm(F_CHM_UBFC, F_CHM_UBEQ, TopOut.Chamber().UBparm())) return 1;
  }

  // DIHEDRAL PK
  if (BufferAlloc(F_DIHPK, TopOut.DihedralParm().size())) return 1;
  for (DihedralParmArray::const_iterator it = TopOut.DihedralParm().begin();
                                         it != TopOut.DihedralParm().end(); ++it)
    file_.DblToBuffer( it->Pk() );
  file_.FlushBuffer();

  // DIHEDRAL PN 
  if (BufferAlloc(F_DIHPN, TopOut.DihedralParm().size())) return 1;
  for (DihedralParmArray::const_iterator it = TopOut.DihedralParm().begin();
                                         it != TopOut.DihedralParm().end(); ++it)
    file_.DblToBuffer( it->Pn() );
  file_.FlushBuffer();

  // DIHEDRAL PHASE
  if (BufferAlloc(F_DIHPHASE, TopOut.DihedralParm().size())) return 1;
  for (DihedralParmArray::const_iterator it = TopOut.DihedralParm().begin();
                                         it != TopOut.DihedralParm().end(); ++it)
    file_.DblToBuffer( it->Phase() );
  file_.FlushBuffer();

  // DIHEDRAL SCEE
  if (BufferAlloc(F_SCEE, TopOut.DihedralParm().size())) return 1;
  for (DihedralParmArray::const_iterator it = TopOut.DihedralParm().begin();
                                         it != TopOut.DihedralParm().end(); ++it)
    file_.DblToBuffer( it->SCEE() );
  file_.FlushBuffer();

  // DIHEDRAL SCNB
  if (BufferAlloc(F_SCNB, TopOut.DihedralParm().size())) return 1;
  for (DihedralParmArray::const_iterator it = TopOut.DihedralParm().begin();
                                         it != TopOut.DihedralParm().end(); ++it)
    file_.DblToBuffer( it->SCNB() );
  file_.FlushBuffer();

  // CHAMBER only - Impropers
  if (ptype_ == CHAMBER) {
    // NUM IMPROPERS
    if (BufferAlloc(F_CHM_NIMP, 1)) return 1;
    file_.IntToBuffer( TopOut.Chamber().Impropers().size() );
    file_.FlushBuffer();
    // IMPROPER TERMS
    if (BufferAlloc(F_CHM_IMP, TopOut.Chamber().Impropers().size()*5)) return 1;
    for (DihedralArray::const_iterator it = TopOut.Chamber().Impropers().begin();
                                       it != TopOut.Chamber().Impropers().end(); ++it)
    {
      file_.IntToBuffer( it->A1() + 1 );
      file_.IntToBuffer( it->A2() + 1 );
      if ( it->Type()  == DihedralType::BOTH || it->Type() == DihedralType::END )
        file_.IntToBuffer( -(it->A3() + 1) );
      else
        file_.IntToBuffer( it->A3() + 1 );
      if ( it->Type() == DihedralType::BOTH || it->Type() == DihedralType::IMPROPER )
        file_.IntToBuffer( -(it->A4() + 1) );
      else
        file_.IntToBuffer( it->A4() + 1 );
      file_.IntToBuffer( it->Idx() + 1 );
    }
    file_.FlushBuffer();
    // NUM IMPROPER PARAMS
    if (BufferAlloc(F_CHM_NIMPT, 1)) return 1;
    file_.IntToBuffer( TopOut.Chamber().ImproperParm().size() );
    file_.FlushBuffer();
    // IMPROPER FORCE CONSTANTS
    if (BufferAlloc(F_CHM_IMPFC, TopOut.Chamber().ImproperParm().size())) return 1;
    for (DihedralParmArray::const_iterator it = TopOut.Chamber().ImproperParm().begin();
                                           it != TopOut.Chamber().ImproperParm().end(); ++it)
      file_.DblToBuffer( it->Pk() );
    file_.FlushBuffer();
    // IMPROPER PHASES
    if (BufferAlloc(F_CHM_IMPP, TopOut.Chamber().ImproperParm().size())) return 1;
    for (DihedralParmArray::const_iterator it = TopOut.Chamber().ImproperParm().begin();
                                           it != TopOut.Chamber().ImproperParm().end(); ++it)
      file_.DblToBuffer( it->Phase() );
    file_.FlushBuffer();
  }

  // SOLTY - Currently unused but must be written.
  if (BufferAlloc(F_SOLTY, TopOut.NatomTypes())) return 1;
  for (int idx = 0; idx != TopOut.NatomTypes(); idx++)
    file_.DblToBuffer( 0.0 );
  file_.FlushBuffer();

  // LJ A and B terms
  if (WriteLJ(F_LJ_A, F_LJ_B, TopOut.Nonbond().NBarray())) return 1;
 
  // CHAMBER only - LJ 1-4 terms
  if (ptype_ == CHAMBER) {
    if (WriteLJ(F_LJ14A, F_LJ14B, TopOut.Chamber().LJ14())) return 1;
  }

  // BONDSH and BONDS
  if (WriteBonds(F_BONDSH, TopOut.BondsH())) return 1;
  if (WriteBonds(F_BONDS,  TopOut.Bonds()) ) return 1;
  // ANGLESH and ANGLES
  if (WriteAngles(F_ANGLESH, TopOut.AnglesH())) return 1;
  if (WriteAngles(F_ANGLES,  TopOut.Angles()) ) return 1;
  // DIHEDRALSH and DIHEDRALS
  if (WriteDihedrals(F_DIHH, TopOut.DihedralsH())) return 1;
  if (WriteDihedrals(F_DIH,  TopOut.Dihedrals()) ) return 1;

  // EXCLUDED ATOMS LIST
  if (BufferAlloc(F_EXCLUDE, Excluded.size())) return 1;
  for (Iarray::const_iterator it = Excluded.begin(); it != Excluded.end(); ++it)
    file_.IntToBuffer( *it );
  file_.FlushBuffer();
  Excluded.clear();

  // HBOND ASOL
  if (BufferAlloc(F_ASOL, TopOut.Nonbond().HBarray().size())) return 1;
  for (HB_ParmArray::const_iterator it = TopOut.Nonbond().HBarray().begin();
                                    it != TopOut.Nonbond().HBarray().end(); ++it)
    file_.DblToBuffer( it->Asol() );
  file_.FlushBuffer();

  // HBOND BSOL
  if (BufferAlloc(F_BSOL, TopOut.Nonbond().HBarray().size())) return 1;
  for (HB_ParmArray::const_iterator it = TopOut.Nonbond().HBarray().begin();
                                    it != TopOut.Nonbond().HBarray().end(); ++it)
    file_.DblToBuffer( it->Bsol() );
  file_.FlushBuffer();

  // HBOND HBCUT 
  if (BufferAlloc(F_HBCUT, TopOut.Nonbond().HBarray().size())) return 1;
  for (HB_ParmArray::const_iterator it = TopOut.Nonbond().HBarray().begin();
                                    it != TopOut.Nonbond().HBarray().end(); ++it)
    file_.DblToBuffer( it->HBcut() );
  file_.FlushBuffer();

  // AMBER ATOM TYPES
  if (BufferAlloc(F_TYPES, TopOut.Natom())) return 1;
  // Check that atoms actually have a type.
  int NemptyTypes = 0;
  for (Topology::atom_iterator atm = TopOut.begin(); atm != TopOut.end(); ++atm) {
    if (atm->Type().len() < 1) ++NemptyTypes;
    file_.CharToBuffer( *(atm->Type()) );
  }
  file_.FlushBuffer();
  if (NemptyTypes > 0) mprintf("Warning: %i empty atom type names.\n", NemptyTypes);

  // TODO: Generate automatically
  // NOTE: These are required by Amber. If they are not present it
  // may mean this topology was not originally read from an Amber topology
  // and could cause problems if used for simulations. Create "empty" values,
  // but warn.
  if (WriteExtra(TopOut, TopOut.Natom())) return 1;

  // CHAMBER only - write CMAP parameters
  if (TopOut.HasCmap()) {
    // CMAP COUNT
    const FlagType f_count = (ptype_ == CHAMBER) ? F_CHM_CMAPC : F_CMAPC;
    if (BufferAlloc(f_count, 2)) return 1;
    file_.IntToBuffer( TopOut.Cmap().size() );     // CMAP terms
    file_.IntToBuffer( TopOut.CmapGrid().size() ); // CMAP grids
    file_.FlushBuffer();
    // CMAP GRID RESOLUTIONS
    const FlagType f_grid_res = (ptype_ == CHAMBER) ? F_CHM_CMAPR : F_CMAPR;
    if (BufferAlloc(f_grid_res, TopOut.CmapGrid().size())) return 1;
    for (CmapGridArray::const_iterator grid = TopOut.CmapGrid().begin();
                                       grid != TopOut.CmapGrid().end(); ++grid)
      file_.IntToBuffer( grid->Resolution() );
    file_.FlushBuffer();
    // CMAP GRIDS
    int ngrid = 1;
    const FlagType f_grid = (ptype_ == CHAMBER) ? F_CHM_CMAPP : F_CMAPP;
    for (CmapGridArray::const_iterator grid = TopOut.CmapGrid().begin();
                                       grid != TopOut.CmapGrid().end();
                                       ++grid, ++ngrid)
    {
      if (BufferAlloc(f_grid, grid->Size(), ngrid)) return 1;
      for (std::vector<double>::const_iterator it = grid->Grid().begin();
                                               it != grid->Grid().end(); ++it)
        file_.DblToBuffer( *it );
      file_.FlushBuffer();
    }
    // CMAP parameters
    const FlagType f_param = (ptype_ == CHAMBER) ? F_CHM_CMAPI : F_CMAPI;
    if (BufferAlloc(f_param, TopOut.Cmap().size())) return 1;
    for (CmapArray::const_iterator it = TopOut.Cmap().begin();
                                   it != TopOut.Cmap().end(); ++it)
    {
      file_.IntToBuffer( it->A1() + 1 );
      file_.IntToBuffer( it->A2() + 1 );
      file_.IntToBuffer( it->A3() + 1 );
      file_.IntToBuffer( it->A4() + 1 );
      file_.IntToBuffer( it->A5() + 1 );
      file_.IntToBuffer( it->Idx() + 1 );
    }
    file_.FlushBuffer();
  }

  // Write solvent info if IFBOX > 0
  if (ifbox > 0) {
    // Need to check if molecules are contiguous or not via
    // how many segments they have.
    bool molsAreContiguous = true;
    for (Topology::mol_iterator mol = TopOut.MolStart(); mol != TopOut.MolEnd(); ++mol)
      if (mol->MolUnit().nSegments() > 1) {
        molsAreContiguous = false;
        break;
      }
    if (molsAreContiguous) {
      // ----- Contiguous molecules --------------
      // Determine first solvent molecule 
      int firstSolventMol = -1;
      for (Topology::mol_iterator mol = TopOut.MolStart(); mol != TopOut.MolEnd(); ++mol) {
        if ( mol->IsSolvent() ) {
          firstSolventMol = (int)(mol - TopOut.MolStart());
          break;
        }
      }
      // Determine final solute residue based on first solvent molecule.
      int finalSoluteRes = 0;
      if (firstSolventMol == -1)
        finalSoluteRes = TopOut.Nres(); // No solvent Molecules
      else if (firstSolventMol > 0) {
        int finalSoluteAtom = TopOut.Mol(firstSolventMol).MolUnit().Front() - 1;
        finalSoluteRes = TopOut[finalSoluteAtom].ResNum() + 1;
      }
      // If no solvent, just set to 1 beyond # of molecules
      if (firstSolventMol == -1)
        firstSolventMol = TopOut.Nmol();
      // SOLVENT POINTERS
      if (BufferAlloc(F_SOLVENT_POINTER, 3)) return 1;
      file_.IntToBuffer( finalSoluteRes ); // Already +1
      file_.IntToBuffer( TopOut.Nmol() );
      file_.IntToBuffer( firstSolventMol + 1 );
      file_.FlushBuffer();
      // ATOMS PER MOLECULE
      if (BufferAlloc(F_ATOMSPERMOL, TopOut.Nmol())) return 1;
      for (Topology::mol_iterator mol = TopOut.MolStart(); mol != TopOut.MolEnd(); mol++)
        file_.IntToBuffer( mol->NumAtoms() );
      file_.FlushBuffer();
    } else {
      // ----- Non-contiguous molecules ----------
      // Molecules are not contiguous. In order to ensure that #s in the
      // ATOMS_PER_MOLECULE section match up, do segments instead.
      mprintf("Warning: Topology has non-contiguous molecules.\n"
              "Warning: Setting up SOLVENT_POINTERS and ATOMS_PER_MOLECULE\n"
              "Warning: based on contiguous segments. The resulting Topology\n"
              "Warning: will not strictly correspond to Amber specifications\n"
              "Warning: and should be used with caution.\n"
              "Warning: The 'fixatomorder' command can be used to reorder\n"
              "Warning: molecules into contiguous segments.\n");
      // Map segment to solvent status
      typedef std::pair<Segment, bool> SegPair;
      typedef std::vector<SegPair> SegArray;
      SegArray segments;
      for (Topology::mol_iterator mol = TopOut.MolStart(); mol != TopOut.MolEnd(); ++mol) {
        for (Unit::const_iterator seg = mol->MolUnit().segBegin();
                                  seg != mol->MolUnit().segEnd(); ++seg)
          segments.push_back( SegPair(*seg, mol->IsSolvent()) );
      }
      std::sort(segments.begin(), segments.end());
      mprintf("\t%zu segments.\n", segments.size());
      // Determine first solvent molecule 
      SegArray::const_iterator firstSolventSeg = segments.end();
      int firstSolventMol = -1;
      for (SegArray::const_iterator it = segments.begin(); it != segments.end(); ++it) {
        // it->second is true if segment was part of solvent molecule
        if ( it->second ) {
          firstSolventSeg = it;
          firstSolventMol = it - segments.begin();
          break;
        }
      }
      // Determine final solute residue based on first solvent molecule.
      int finalSoluteRes = 0;
      if (firstSolventSeg == segments.end())
        finalSoluteRes = TopOut.Nres(); // No solvent Molecules
      else if (firstSolventSeg != segments.begin()) {
        int finalSoluteAtom = firstSolventSeg->first.Begin() - 1;
        finalSoluteRes = TopOut[finalSoluteAtom].ResNum() + 1;
      }
      // If no solvent, just set to 1 beyond # of segments 
      if (firstSolventMol == -1)
        firstSolventMol = segments.size();
      // SOLVENT POINTERS
      if (BufferAlloc(F_SOLVENT_POINTER, 3)) return 1;
      file_.IntToBuffer( finalSoluteRes ); // Already +1
      file_.IntToBuffer( segments.size() );
      file_.IntToBuffer( firstSolventMol + 1 );
      file_.FlushBuffer();
      // ATOMS PER MOLECULE
      if (BufferAlloc(F_ATOMSPERMOL, segments.size())) return 1;
      for (SegArray::const_iterator it = segments.begin(); it != segments.end(); ++it)
        file_.IntToBuffer( it->first.Size() );
      file_.FlushBuffer();
    }
    // BOX DIMENSIONS
    if (BufferAlloc(F_PARMBOX, 4)) return 1;
    double beta;
    // Special case: Rhombic dodecahedron stores alpha/gamma (60) instead of beta (90)
    if (cellShape == Box::RHOMBIC_DODECAHEDRON)
      beta = 60.0;
    else
      beta = TopOut.ParmBox().Param(Box::BETA);
    file_.DblToBuffer( beta );
    file_.DblToBuffer( TopOut.ParmBox().Param(Box::X) );
    file_.DblToBuffer( TopOut.ParmBox().Param(Box::Y) );
    file_.DblToBuffer( TopOut.ParmBox().Param(Box::Z) );
    file_.FlushBuffer();
  }

  // CAP info
  if (TopOut.Cap().NatCap() > 0) {
    if (BufferAlloc(F_CAP_INFO, 1)) return 1;
    file_.IntToBuffer( TopOut.Cap().NatCap()+1 );
    file_.FlushBuffer();
    if (BufferAlloc(F_CAP_INFO2, 4)) return 1;
    file_.DblToBuffer( TopOut.Cap().CutCap() );
    file_.DblToBuffer( TopOut.Cap().xCap() );
    file_.DblToBuffer( TopOut.Cap().yCap() );
    file_.DblToBuffer( TopOut.Cap().zCap() );
    file_.FlushBuffer();
  }

  // Only write GB params if actually present. At least one atom must have something.
  bool hasGB = false;
  for (Topology::atom_iterator atm = TopOut.begin(); atm != TopOut.end(); ++atm)
    if (atm->GBRadius() > 0.0 || atm->Screen() > 0.0) { // TODO check negative?
      hasGB = true;
      break;
    }
  if (hasGB) {
    // GB RADIUS SET
    if (!TopOut.GBradiiSet().empty())
      WriteLine(F_RADSET, TopOut.GBradiiSet());
    // GB RADII
    if (BufferAlloc(F_RADII, TopOut.Natom())) return 1;
    for (Topology::atom_iterator atm = TopOut.begin(); atm != TopOut.end(); ++atm)
      file_.DblToBuffer( atm->GBRadius() );
    file_.FlushBuffer();
    // GB SCREENING PARAMS
    if (BufferAlloc(F_SCREEN, TopOut.Natom())) return 1;
    for (Topology::atom_iterator atm = TopOut.begin(); atm != TopOut.end(); ++atm)
      file_.DblToBuffer( atm->Screen() );
    file_.FlushBuffer();
  }

  // Polarizability - only write if it needs to be there
  if (TopOut.Ipol() > 0) {
    if (BufferAlloc(F_IPOL, 1)) return 1;
    file_.IntToBuffer( TopOut.Ipol() );
    file_.FlushBuffer();
    if (BufferAlloc(F_POLAR, TopOut.Natom())) return 1;
    for (Topology::atom_iterator atm = TopOut.begin(); atm != TopOut.end(); ++atm)
      file_.DblToBuffer( atm->Polar() );
    file_.FlushBuffer();
  }

  // LES parameters
  if (TopOut.LES().HasLES()) {
    // LES NTYP
    if (BufferAlloc(F_LES_NTYP, 1)) return 1;
    file_.IntToBuffer( TopOut.LES().Ntypes() );
    file_.FlushBuffer();
    // Sanity check.
    if ( (int)TopOut.LES().Array().size() != TopOut.Natom() ) {
      mprinterr("Internal Error: # LES atoms (%zu) != # atoms in topology (%i).\n",
                TopOut.LES().Array().size(), TopOut.Natom());
      return 1;
    }
    // LES TYPE
    if (BufferAlloc(F_LES_TYPE, TopOut.LES().Array().size())) return 1;
    for (LES_Array::const_iterator les = TopOut.LES().Array().begin();
                                   les != TopOut.LES().Array().end(); ++les)
      file_.IntToBuffer( les->Type() );
    file_.FlushBuffer();
    // LES FAC
    if (BufferAlloc(F_LES_FAC, TopOut.LES().FAC().size())) return 1;
    for (std::vector<double>::const_iterator it = TopOut.LES().FAC().begin();
                                             it != TopOut.LES().FAC().end(); ++it)
      file_.DblToBuffer( *it );
    file_.FlushBuffer();
    // LES CNUM
    if (BufferAlloc(F_LES_CNUM, TopOut.LES().Array().size())) return 1;
    for (LES_Array::const_iterator les = TopOut.LES().Array().begin();
                                   les != TopOut.LES().Array().end(); ++les)
      file_.IntToBuffer( les->Copy() );
    file_.FlushBuffer();
    // LES ID
    if (BufferAlloc(F_LES_ID, TopOut.LES().Array().size())) return 1;
    for (LES_Array::const_iterator les = TopOut.LES().Array().begin();
                                   les != TopOut.LES().Array().end(); ++les)
      file_.IntToBuffer( les->ID() );
    file_.FlushBuffer();
  }

  if (hasOrigResNums || hasChainID || hasIcodes) {
    // PDB residue numbers. Need to adjust format based on res # digit width.
    int maxWidth = std::max( DigitWidth(TopOut.Res(0).OriginalResNum() ),
                             DigitWidth(TopOut.Res(TopOut.Nres()-1).OriginalResNum()) );
    //mprintf("DEBUG: Max original res # digit width is %i\n", maxWidth);
    int err;
    if ( maxWidth < 5 )
      err = BufferAlloc(F_PDB_RES, TopOut.Nres());
    else if ( maxWidth < 9 ) {
      FortranData FMT;
      FMT.ParseFortranFormat("%FORMAT(10I8)");
      err = BufferAlloc(F_PDB_RES, FMT, TopOut.Nres(), -1);
    } else {
      // Max 32 bit int is 10 digits
      FortranData FMT;
      FMT.ParseFortranFormat("%FORMAT(5I16)");
      err = BufferAlloc(F_PDB_RES, FMT, TopOut.Nres(), -1);
    }
    if (err != 0) return 1;
    for (Topology::res_iterator res = TopOut.ResStart(); res != TopOut.ResEnd(); ++res)
      file_.IntToBuffer( res->OriginalResNum() );
    file_.FlushBuffer();
  }
  if (hasChainID) {
    // PDB chain IDs
    if (BufferAlloc(F_PDB_CHAIN, TopOut.Nres())) return 1;
    char cid[2];
    cid[1] = '\0';
    for (Topology::res_iterator res = TopOut.ResStart(); res != TopOut.ResEnd(); ++res)
    {
      cid[0] = res->ChainId();
      file_.CharToBuffer( cid );
    }
    file_.FlushBuffer();
  }
  if (hasIcodes) {
    // PDB residue insertion codes
    if (BufferAlloc(F_PDB_ICODE, TopOut.Nres())) return 1;
    char icode[2];
    icode[1] = '\0';
    for (Topology::res_iterator res = TopOut.ResStart(); res != TopOut.ResEnd(); ++res)
    {
      icode[0] = res->Icode();
      file_.CharToBuffer( icode );
    }
    file_.FlushBuffer();
  }
  if (hasOcc) {
    // PDB atomic occupancy
    if (BufferAlloc(F_PDB_OCC, TopOut.Natom())) return 1;
    for (std::vector<float>::const_iterator it = TopOut.Occupancy().begin();
                                            it != TopOut.Occupancy().end(); ++it)
      file_.DblToBuffer( *it );
    file_.FlushBuffer();
  }
  if (hasBfac) {
    // PDB atomic B-factors
    if (BufferAlloc(F_PDB_BFAC, TopOut.Natom())) return 1;
    for (std::vector<float>::const_iterator it = TopOut.Bfactor().begin();
                                            it != TopOut.Bfactor().end(); ++it)
      file_.DblToBuffer( *it );
    file_.FlushBuffer();
  }
  if (hasNum) {
    // PDB original serial numbers
    if (BufferAlloc(F_PDB_NUM, TopOut.Natom())) return 1;
    for (std::vector<int>::const_iterator it = TopOut.PdbSerialNum().begin();
                                          it != TopOut.PdbSerialNum().end(); ++it)
      file_.IntToBuffer( *it );
    file_.FlushBuffer();
  }

  // LJ C terms.
  // NOTE: I put this at the very end to match behavior of parmed add_12_6_4
  if (TopOut.Nonbond().Has_C_Coeff()) {
    std::vector<double> const& ccoef = TopOut.Nonbond().LJC_Array();
    if (BufferAlloc(F_LJ_C, ccoef.size())) return 1;
    for (std::vector<double>::const_iterator it = ccoef.begin(); it != ccoef.end(); ++it)
      file_.DblToBuffer( *it );
    file_.FlushBuffer();
  }

  return 0;
}

// =============================================================================
Parm_Amber::FortranData::FortranData(const char* ptrIn) :
  ftype_(UNKNOWN_FTYPE), fncols_(0), fwidth_(0), fprecision_(0)
{
  ParseFortranFormat(ptrIn);
}

/** Given a fortran-type format string, set the corresponding fortran
  * type. Set fncols (if present), fwidth, and fprecision (if present).
  */
int Parm_Amber::FortranData::ParseFortranFormat(const char* ptrIn) {
  if (ptrIn == 0) {
    mprinterr("Error: Empty format string.\n");
    return 1;
  }
  fstr_ = ptrIn;
  std::string fformat( NoTrailingWhitespace( ptrIn ) );
  if ( fformat.empty() ) return 1;
  //mprintf("DEBUG: Fortran format: %s\n", fformat.c_str());
  // Make sure characters are upper case.
  for (std::string::iterator p = fformat.begin(); p != fformat.end(); p++)
    *p = toupper(*p);
  // Advance past left parentheses
  std::string::iterator ptr = fformat.begin() + 7;
  while (*ptr=='(') ++ptr;
  // If digit, have number of data columns. Min # is 1
  std::string arg;
  fncols_ = 1;
  if (isdigit(*ptr)) {
    while (ptr!=fformat.end() && isdigit(*ptr)) {
      arg += *ptr;
      ++ptr;
    }
    fncols_ = atoi( arg.c_str() );
  }
  // Advance past any more left parentheses
  while (ptr!=fformat.end() && *ptr=='(') ++ptr;
  // Type
  if (ptr==fformat.end()) {
    mprinterr("Error: Malformed fortran format string (%s)\n", fformat.c_str());
    return 1;
  }
  switch (*ptr) {
    case 'I' : ftype_ = FINT;    break;
    case 'E' : ftype_ = FDOUBLE; break;
    case 'A' : ftype_ = FCHAR;   break;
    case 'F' : ftype_ = FFLOAT;  break;
    default  : ftype_ = UNKNOWN_FTYPE;
  }
  ++ptr;
  // Width
  fwidth_ = 0;
  arg.clear();
  while (isdigit(*ptr)) {
    arg += *ptr;
    ++ptr;
  }
  fwidth_ = atoi( arg.c_str() );
  // Precision
  fprecision_ = 0;
  if (*ptr == '.') {
    ++ptr;
    arg.clear();
    while (isdigit(*ptr)) {
      arg += *ptr;
      ++ptr;
    }
    fprecision_ = atoi( arg.c_str() );
  }
  //if (debug_ > 2)
  //  mprintf("[%s]: cols=%i type=%i width=%i precision=%i\n",fformat.c_str(),
  //          fncols_,(int)ftype_,fwidth_,fprecision_);

  return 0;
}
