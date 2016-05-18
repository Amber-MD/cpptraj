#include <cstring> // strncmp
#include <cstdio> // sscanf
#include <cctype> // toupper, isdigit
#include <cstdlib>  // atoi
#include "Parm_Amber.h"
#include "CpptrajStdio.h"
#include "Constants.h" // ELECTOAMBER, AMBERTOELEC
#include "StringRoutines.h" // NoTrailingWhitespace

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
/// Constant strings for Amber parm flags and fortran formats.
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
  { "RESIDUE_NUMBER", "%FORMAT(20I4)" }, // PDB residue number
  { "RESIDUE_CHAINID", F20a4 }, // PDB chain ID
  { "RESIDUE_ICODE", F20a4 },   // PDB residue insertion code
  { "ATOM_ALTLOC", F20a4 },     // PDB atom alt location indicator FIXME: format is guess
  { 0, 0 }
};

// -----------------------------------------------------------------------------
// CONSTRUCTOR
Parm_Amber::Parm_Amber() :
  ptype_(OLDPARM),
  numLJparm_(0),
  SCEE_set_(false),
  SCNB_set_(false),
  N_impropers_(0),
  N_impTerms_(0),
  nochamber_(false)
{
  UB_count_[0] = 0;
  UB_count_[1] = 0;
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
      if ( sscanf(lineBuf,"%6i%6i%6i%6i%6i%6i%6i%6i%6i%6i%6i%6i",
                  iamber,   iamber+1, iamber+2,  iamber+3,
                  iamber+4, iamber+5, iamber+6,  iamber+7,
                  iamber+8, iamber+9, iamber+10, iamber+11) == 12 )
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
  if (infile_.OpenRead( fname )) return 1;
  int err = 0;
  if (ptype_ == OLDPARM)
    err = ReadOldParm( TopIn );
  else
    err = ReadNewParm( TopIn );
  if (err != 0) return 1;
  // Determine Atom elements
  if (atomicNums_.empty()) {
    mprintf("Warning: Amber topology does not include atomic numbers.\n"
            "Warning: Will attempt to assign elements from atom masses/names.\n");
    atomicNums_.assign( values_[NATOM], 0 );
  }
  for (int idx = 0; idx != values_[NATOM]; idx++)
    TopIn.SetAtom(idx).DetermineElement( atomicNums_[idx] );
  // Set Atom residue numbers
  for (Topology::res_iterator res = TopIn.ResStart(); res != TopIn.ResEnd(); ++res)
    for (int at = res->FirstAtom(); at != res->LastAtom(); ++at)
      TopIn.SetAtom(at).SetResNum( res - TopIn.ResStart() );
  // If SCEE/SCNB not read, set defaults
  if (!SCEE_set_)
    for (unsigned int idx = 0; idx != TopIn.DihedralParm().size(); idx++)
      TopIn.SetDihedralParm(idx).SetSCEE( 1.2 );
  if (!SCNB_set_)
    for (unsigned int idx = 0; idx != TopIn.DihedralParm().size(); idx++)
      TopIn.SetDihedralParm(idx).SetSCNB( 2.0 );
  // Check box info
  if (values_[IFBOX] > 0) {
    if (parmbox_.Type() == Box::NOBOX) {
      if (ptype_ != CHAMBER) mprintf("Warning: Prmtop missing Box information.\n");
      // ifbox 2: truncated octahedron for certain
      if (values_[IFBOX] == 2)
        parmbox_.SetTruncOct();
    }
  }
  // Check for IFBOX/BoxType mismatch
  if (values_[IFBOX]==2 && parmbox_.Type() != Box::TRUNCOCT) {
    mprintf("Warning: Amber Parm Box should be Truncated Octahedron (ifbox==2)\n");
    mprintf("         but BOX_DIMENSIONS indicate %s - may cause imaging problems.\n",
            parmbox_.TypeName());
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
  std::string title = NoTrailingWhitespace( infile_.GetLine() );
  mprintf("DEBUG: Title= '%s'\n", title.c_str());
  int Npointers = 30; // No NEXTRA etc
  FortranData DBL(FDOUBLE, 5, 16, 0);
  FortranData INT(FINT, 12, 6, 0);
  FortranData CHAR(FCHAR, 20, 4, 0);
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
    infile_.NextElement(); // Final solute residue
    int nmolecules = atoi(infile_.NextElement());
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
  const char* ptr = infile_.NextLine();
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
        //if (debug_ > 0)
        mprintf("DEBUG: Version: %s\n", NoTrailingWhitespace(ptr).c_str());
        ptr = infile_.NextLine();
      } else if (IsFLAG(ptr)) {
        // %FLAG <type> line. Determine the flag type.
        std::string flagType = NoTrailingWhitespace(ptr+6);
        //mprintf("DEBUG: Flag type: %s\n", flagType.c_str());
        int flagIdx = -1;
        if ( flagType.compare(0, 22, "CHARMM_CMAP_PARAMETER_") == 0 ) {
          // Special case. This flag has a 2 digit extension.
          flagIdx = (int)F_CHM_CMAPP;
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
          switch ((AmberParmFlagType)flagIdx) {
            case F_CTITLE: ptype_ = CHAMBER; // Fall through to F_TITLE
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
            // CHAMBER
            case F_FF_TYPE:   err = ReadChamberFFtype(TopIn); break;
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
            case F_CHM_CMAPC: err = ReadChamberCmapCounts(FMT); break;
            case F_CHM_CMAPR: err = ReadChamberCmapRes(TopIn, FMT); break;
            case F_CHM_CMAPP: err = ReadChamberCmapGrid(flagType.c_str(), TopIn, FMT); break;
            case F_CHM_CMAPI: err = ReadChamberCmapTerms(TopIn, FMT); break;
            // LES
            case F_LES_NTYP:  err = ReadLESntyp(TopIn, FMT); break;
            case F_LES_TYPE:  err = ReadLEStypes(TopIn, FMT); break;
            case F_LES_FAC:   err = ReadLESfac(TopIn, FMT); break;
            case F_LES_CNUM:  err = ReadLEScnum(TopIn, FMT); break;
            case F_LES_ID:    err = ReadLESid(TopIn, FMT); break;
            // Sanity check
            default:
              mprinterr("Internal Error: Unhandled FLAG '%s'.\n",flagType.c_str());
              return 1;
          }
          if (err != 0) return 1;
        }
      } else {
        // Unknown '%' tag. Read past it.
        ptr = infile_.NextLine();
      }  
    } else {
      // Unknown line. Read past it.
      ptr = infile_.NextLine();
    }
  }

  infile_.CloseFile();
  return 0;
}

// Parm_Amber::SkipToNextFlag()
const char* Parm_Amber::SkipToNextFlag() {
  const char* ptr = infile_.NextLine();
  while (ptr != 0 && !IsFLAG(ptr)) ptr = infile_.NextLine();
  return ptr;
}

// Parm_Amber::ReadFormatLine()
int Parm_Amber::ReadFormatLine(FortranData& FMT) {
  // Read next line; can be either a COMMENT or FORMAT. If COMMENT, 
  // read past until you get to the FORMAT line.
  const char* ptr = infile_.NextLine();
  if (ptr == 0) {
    mprinterr("Error: Unexpected EOF in Amber Topology when looking for FORMAT.\n");
    return 1;
  }
  while ( ptr != 0 && strncmp(ptr, "%FORMAT", 7) !=0) {
    ptr = infile_.NextLine();
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

// Parm_Amber::ReadTitle()
int Parm_Amber::ReadTitle(Topology& TopIn) {
  std::string title = NoTrailingWhitespace( infile_.GetLine() );
  //if (debug_>0)
    mprintf("\tAmberParm Title: \"%s\"\n",title.c_str());
  TopIn.SetParmName( title, infile_.Filename() );
  if (infile_.NextLine() == 0) return 1; // Advance to next line
  return 0;
}

// Parm_Amber::ReadPointers()
int Parm_Amber::ReadPointers(int Npointers, Topology& TopIn, FortranData const& FMT) {
  infile_.SetupFrameBuffer( Npointers, FMT.Width(), FMT.Ncols() );
  if (infile_.ReadFrame()) return 1;
  values_.reserve(Npointers);
  for (int idx = 0; idx != Npointers; idx++)
    values_.push_back( atoi( infile_.NextElement() ) );

  mprintf("DEBUG: POINTERS\n");
  for (Iarray::const_iterator it = values_.begin(); it != values_.end(); ++it)
    mprintf("%u\t%i\n", it-values_.begin(), *it);

  TopIn.Resize( Topology::Pointers(values_[NATOM], values_[NRES], values_[NATOM],
                                   values_[NUMBND], values_[NUMANG], values_[NPTRA]) );

  if (values_[IFPERT] > 0)
    mprintf("Warning: '%s' contains perturbation information.\n"
            "Warning:  Cpptraj currently does not read or write perturbation information.\n",
            infile_.Filename().base());

  numLJparm_ = values_[NTYPES] * (values_[NTYPES]+1) / 2;
  TopIn.SetNonbond().SetNLJterms( numLJparm_ );
  TopIn.SetNonbond().SetNHBterms( values_[NPHB] );
  return 0;
}

// Parm_Amber::SetupBuffer()
int Parm_Amber::SetupBuffer(AmberParmFlagType ftype, int nvals, FortranData const& FMT) {
  if (values_.empty()) {
    mprinterr("Error: Flag '%s' encountered before POINTERS.\n", FLAGS_[ftype].Flag);
    return 1;
  }
  if (nvals > 0) {
    mprintf("DEBUG: Set up buffer for '%s', %i vals.\n", FLAGS_[ftype].Flag, nvals);
    infile_.SetupFrameBuffer( nvals, FMT.Width(), FMT.Ncols() );
    if (infile_.ReadFrame()) return 1;
    //mprintf("DEBUG: '%s':\n%s", FLAGS_[ftype].Flag, infile_.Buffer());
  } else {
    mprintf("DEBUG: No values for flag '%s'\n", FLAGS_[ftype].Flag);
    // Read blank line
    infile_.NextLine();
  }
  return 0;
}

// Parm_Amber::ReadAtomNames()
int Parm_Amber::ReadAtomNames(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_NAMES, values_[NATOM], FMT)) return 1;
  for (int idx = 0; idx != values_[NATOM]; idx++)
    TopIn.SetAtom(idx).SetName( NameType(infile_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadAtomCharges()
/** Read atomic charges. Convert units to elec. */
int Parm_Amber::ReadAtomCharges(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_CHARGE, values_[NATOM], FMT)) return 1;
  for (int idx = 0; idx != values_[NATOM]; idx++)
    TopIn.SetAtom(idx).SetCharge( atof(infile_.NextElement()) * Constants::AMBERTOELEC );
  return 0;
}

// Parm_Amber::ReadAtomicNum()
/** Read atomic numbers to a temporary array. This is done because some
  * topology files do not have atomic number information.
  */
int Parm_Amber::ReadAtomicNum(FortranData const& FMT) {
  if (SetupBuffer(F_ATOMICNUM, values_[NATOM], FMT)) return 1;
  for (int idx = 0; idx != values_[NATOM]; idx++)
    atomicNums_.push_back( atoi(infile_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadAtomicMass()
int Parm_Amber::ReadAtomicMass(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_MASS, values_[NATOM], FMT)) return 1;
  for (int idx = 0; idx != values_[NATOM]; idx++)
    TopIn.SetAtom(idx).SetMass( atof(infile_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadAtomTypeIndex()
/** Read atom type indices. Shift by -1 to match CPPTRAJ internal indexing.
  * NOTE: If ever used, shift atom #s in excludedAtoms by -1 so they start from 0
  */
int Parm_Amber::ReadAtomTypeIndex(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_ATYPEIDX, values_[NATOM], FMT)) return 1;
  for (int idx = 0; idx != values_[NATOM]; idx++)
    TopIn.SetAtom(idx).SetTypeIndex( atoi(infile_.NextElement()) - 1 );
  return 0;
}

// Parm_Amber::ReadNonbondIndices()
int Parm_Amber::ReadNonbondIndices(Topology& TopIn, FortranData const& FMT) {
  int nvals = values_[NTYPES]*values_[NTYPES];
  if (SetupBuffer(F_NB_INDEX, nvals, FMT)) return 1;
  TopIn.SetNonbond().SetNtypes( values_[NTYPES] );
  for (int idx = 0; idx != nvals; idx++)
  {
    // Shift positive indices in NONBONDED index array by -1.
    int nbidx = atoi(infile_.NextElement());
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
    TopIn.SetRes(idx).SetName( NameType(infile_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadResidueAtomNums()
int Parm_Amber::ReadResidueAtomNums(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_RESNUMS, values_[NRES], FMT)) return 1;
  for (int idx = 0; idx != values_[NRES]; idx++) {
    int atnum = atoi(infile_.NextElement()) - 1;
    TopIn.SetRes(idx).SetFirstAtom( atnum );
    if (idx > 0) TopIn.SetRes(idx-1).SetLastAtom( atnum );
    TopIn.SetRes(idx).SetOriginalNum( idx+1 );
  }
  TopIn.SetRes( values_[NRES]-1 ).SetLastAtom( values_[NATOM] );
  return 0;
}

// Parm_Amber::ReadBondRK()
int Parm_Amber::ReadBondRK(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_BONDRK, values_[NUMBND], FMT)) return 1;
  for (int idx = 0; idx != values_[NUMBND]; idx++)
    TopIn.SetBondParm(idx).SetRk( atof(infile_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadBondREQ()
int Parm_Amber::ReadBondREQ(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_BONDREQ, values_[NUMBND], FMT)) return 1;
  for (int idx = 0; idx != values_[NUMBND]; idx++)
    TopIn.SetBondParm(idx).SetReq( atof(infile_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadAngleTK()
int Parm_Amber::ReadAngleTK(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_ANGLETK, values_[NUMANG], FMT)) return 1;
  for (int idx = 0; idx != values_[NUMANG]; idx++)
    TopIn.SetAngleParm(idx).SetTk( atof(infile_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadAngleTEQ()
int Parm_Amber::ReadAngleTEQ(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_ANGLETEQ, values_[NUMANG], FMT)) return 1;
  for (int idx = 0; idx != values_[NUMANG]; idx++)
    TopIn.SetAngleParm(idx).SetTeq( atof(infile_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadDihedralPK()
int Parm_Amber::ReadDihedralPK(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_DIHPK, values_[NPTRA], FMT)) return 1;
  for (int idx = 0; idx != values_[NPTRA]; idx++)
    TopIn.SetDihedralParm(idx).SetPk( atof(infile_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadDihedralPN()
int Parm_Amber::ReadDihedralPN(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_DIHPN, values_[NPTRA], FMT)) return 1;
  for (int idx = 0; idx != values_[NPTRA]; idx++)
    TopIn.SetDihedralParm(idx).SetPn( atof(infile_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadDihedralPHASE()
int Parm_Amber::ReadDihedralPHASE(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_DIHPHASE, values_[NPTRA], FMT)) return 1;
  for (int idx = 0; idx != values_[NPTRA]; idx++)
    TopIn.SetDihedralParm(idx).SetPhase( atof(infile_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadDihedralSCEE()
int Parm_Amber::ReadDihedralSCEE(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_SCEE, values_[NPTRA], FMT)) return 1;
  for (int idx = 0; idx != values_[NPTRA]; idx++)
    TopIn.SetDihedralParm(idx).SetSCEE( atof(infile_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadDihedralSCNB()
int Parm_Amber::ReadDihedralSCNB(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_SCNB, values_[NPTRA], FMT)) return 1;
  for (int idx = 0; idx != values_[NPTRA]; idx++)
    TopIn.SetDihedralParm(idx).SetSCNB( atof(infile_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadLJA()
int Parm_Amber::ReadLJA(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_LJ_A, numLJparm_, FMT)) return 1;
  for (int idx = 0; idx != numLJparm_; idx++)
    TopIn.SetNonbond().SetLJ(idx).SetA( atof(infile_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadLJB()
int Parm_Amber::ReadLJB(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_LJ_B, numLJparm_, FMT)) return 1;
  for (int idx = 0; idx != numLJparm_; idx++)
    TopIn.SetNonbond().SetLJ(idx).SetB( atof(infile_.NextElement()) );
  return 0;
}

// Parm_Amber::GetBond()
/** Amber bond indices are * 3, bond parm indices are +1.
  * NOTE: Since NextElement() only returns a pointer to the file buffer,
  *       have to convert a number as soon as it is available.
  */
BondType Parm_Amber::GetBond() {
  int a1 = atoi(infile_.NextElement());
  int a2 = atoi(infile_.NextElement());
  int bidx = atoi(infile_.NextElement());
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
  int a1 = atoi(infile_.NextElement());
  int a2 = atoi(infile_.NextElement());
  int a3 = atoi(infile_.NextElement());
  int aidx = atoi(infile_.NextElement());
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
  int a1 = atoi(infile_.NextElement());
  int a2 = atoi(infile_.NextElement());
  int a3 = atoi(infile_.NextElement());
  int a4 = atoi(infile_.NextElement());
  int didx = atoi(infile_.NextElement());
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
    TopIn.SetNonbond().SetHB(idx).SetAsol( atof(infile_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadBsol()
int Parm_Amber::ReadBsol(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_BSOL, values_[NPHB], FMT)) return 1;
  for (int idx = 0; idx != values_[NPHB]; idx++)
    TopIn.SetNonbond().SetHB(idx).SetBsol( atof(infile_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadHBcut()
int Parm_Amber::ReadHBcut(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_HBCUT, values_[NPHB], FMT)) return 1;
  for (int idx = 0; idx != values_[NPHB]; idx++)
    TopIn.SetNonbond().SetHB(idx).SetHBcut( atof(infile_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadAtomTypes()
int Parm_Amber::ReadAtomTypes(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_TYPES, values_[NATOM], FMT)) return 1;
  for (int idx = 0; idx != values_[NATOM]; idx++)
    TopIn.SetAtom(idx).SetTypeName( NameType(infile_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadItree()
int Parm_Amber::ReadItree(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_ITREE, values_[NATOM], FMT)) return 1;
  for (int idx = 0; idx != values_[NATOM]; idx++)
    TopIn.SetExtraAtomInfo(idx).SetItree( NameType(infile_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadJoin()
int Parm_Amber::ReadJoin(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_JOIN, values_[NATOM], FMT)) return 1;
  for (int idx = 0; idx != values_[NATOM]; idx++)
    TopIn.SetExtraAtomInfo(idx).SetJoin( atoi(infile_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadIrotat()
int Parm_Amber::ReadIrotat(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_IROTAT, values_[NATOM], FMT)) return 1;
  for (int idx = 0; idx != values_[NATOM]; idx++)
    TopIn.SetExtraAtomInfo(idx).SetIrotat( atoi(infile_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadBox()
int Parm_Amber::ReadBox(FortranData const& FMT) {
  if (SetupBuffer(F_PARMBOX, 4, FMT)) return 1;
  double beta = atof(infile_.NextElement());
  double bx = atof(infile_.NextElement());
  double by = atof(infile_.NextElement());
  double bz = atof(infile_.NextElement());
  parmbox_.SetBetaLengths( beta, bx, by, bz );
  return 0;
}

// Parm_Amber::ReadCapInfo()
int Parm_Amber::ReadCapInfo(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_CAP_INFO, 1, FMT)) return 1;
  TopIn.SetCap().SetNatcap( atoi(infile_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadCapInfo2()
int Parm_Amber::ReadCapInfo2(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_CAP_INFO2, 4, FMT)) return 1;
  TopIn.SetCap().SetCutCap( atof(infile_.NextElement()) );
  TopIn.SetCap().SetXcap( atof(infile_.NextElement()) );
  TopIn.SetCap().SetYcap( atof(infile_.NextElement()) );
  TopIn.SetCap().SetZcap( atof(infile_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadGBradiiSet()
int Parm_Amber::ReadGBradiiSet(Topology& TopIn) {
  std::string radius_set = NoTrailingWhitespace(infile_.GetLine());
  mprintf("\tRadius Set: %s\n",radius_set.c_str());
  TopIn.SetGBradiiSet( radius_set );
  return 0;
}

// Parm_Amber::ReadGBradii()
int Parm_Amber::ReadGBradii(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_RADII, values_[NATOM], FMT)) return 1;
  for (int idx = 0; idx != values_[NATOM]; idx++)
    TopIn.SetAtom(idx).SetGBradius( atof(infile_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadGBscreen()
int Parm_Amber::ReadGBscreen(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_SCREEN, values_[NATOM], FMT)) return 1;
  for (int idx = 0; idx != values_[NATOM]; idx++)
    TopIn.SetAtom(idx).SetGBscreen( atof(infile_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadIpol()
int Parm_Amber::ReadIpol(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_IPOL, 1, FMT)) return 1;
  TopIn.SetIpol( atoi(infile_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadPolar()
int Parm_Amber::ReadPolar(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_POLAR, values_[NATOM], FMT)) return 1;
  for (int idx = 0; idx != values_[NATOM]; idx++)
    TopIn.SetAtom(idx).SetPolar( atof(infile_.NextElement()) );
  return 0;
}

// ----- Extra PDB Info --------------------------
int Parm_Amber::ReadPdbRes(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_PDB_RES, values_[NRES], FMT)) return 1;
  for (int idx = 0; idx != values_[NRES]; idx++)
    TopIn.SetRes(idx).SetOriginalNum( atoi(infile_.NextElement()) );
  return 0;
}

int Parm_Amber::ReadPdbChainID(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_PDB_CHAIN, values_[NRES], FMT)) return 1;
  for (int idx = 0; idx != values_[NRES]; idx++)
    TopIn.SetRes(idx).SetChainID( *(infile_.NextElement()) );
  return 0;
}

int Parm_Amber::ReadPdbIcode(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_PDB_ICODE, values_[NRES], FMT)) return 1;
  for (int idx = 0; idx != values_[NRES]; idx++)
    TopIn.SetRes(idx).SetIcode( *(infile_.NextElement()) );
  return 0;
}

int Parm_Amber::ReadPdbAlt(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_PDB_ALT, values_[NATOM], FMT)) return 1;
  for (int idx = 0; idx != values_[NATOM]; idx++)
    TopIn.SetExtraAtomInfo(idx).SetAltLoc( *(infile_.NextElement()) );
  return 0;
}

// ----- CHAMBER ---------------------------------
// ReadChamberFFtype(Topology& TopIn)
int Parm_Amber::ReadChamberFFtype(Topology& TopIn) {
  const char* ptr = infile_.NextLine();
  char ff_verstr[2];
  ff_verstr[0] = ptr[0];
  ff_verstr[1] = ptr[1];
  int ff_verno = atoi(ff_verstr);
  std::string fftype = NoTrailingWhitespace( ptr+2 );
  TopIn.SetChamber().SetVersion( ff_verno, fftype );
  mprintf("\tCHAMBER topology: %i: %s\n", ff_verno, fftype.c_str());
  TopIn.SetChamber().SetNLJ14terms( numLJparm_ );
  return 0;
}

// Parm_Amber::ReadChamberUBCount()
int Parm_Amber::ReadChamberUBCount(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_CHM_UBC, 2, FMT)) return 1;
  UB_count_[0] = atoi(infile_.NextElement()); // Number of bonds
  UB_count_[1] = atoi(infile_.NextElement()); // Number of parameters
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
    int a1 = atoi(infile_.NextElement()) - 1;
    int a2 = atoi(infile_.NextElement()) - 1;
    int bidx = atoi(infile_.NextElement()) - 1;
    TopIn.SetChamber().AddUBterm( BondType(a1, a2, bidx) );
  }
  return 0;
}

// Parm_Amber::ReadChamberUBFC()
int Parm_Amber::ReadChamberUBFC(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_CHM_UBFC, UB_count_[1], FMT)) return 1;
  for (int idx = 0; idx != UB_count_[1]; idx++)
    TopIn.SetChamber().SetUBparm(idx).SetRk( atof(infile_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadChamberUBEQ()
int Parm_Amber::ReadChamberUBEQ(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_CHM_UBEQ, UB_count_[1], FMT)) return 1;
  for (int idx = 0; idx != UB_count_[1]; idx++)
    TopIn.SetChamber().SetUBparm(idx).SetReq( atof(infile_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadChamberNumImpropers()
int Parm_Amber::ReadChamberNumImpropers(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_CHM_NIMP, 1, FMT)) return 1;
  N_impropers_ = atoi(infile_.NextElement());
  TopIn.SetChamber().ReserveImproperTerms( N_impropers_ );
  N_impropers_ *= 5; // Number of improper terms (impropers * 5)
  return 0;
}

// Parm_Amber::ReadChamberNumImpTerms()
int Parm_Amber::ReadChamberNumImpTerms(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_CHM_NIMPT, 1, FMT)) return 1;
  N_impTerms_ = atoi(infile_.NextElement());
  TopIn.SetChamber().ResizeImproperParm( N_impTerms_ );
  return 0;
}

// Parm_Amber::ReadChamberImpropers()
int Parm_Amber::ReadChamberImpropers(Topology& TopIn, FortranData const& FMT) {
  // NOTE: N_impropers_ already multiplied by 5
  if (SetupBuffer(F_CHM_IMP, N_impropers_, FMT)) return 1;
  for (int idx = 0; idx != N_impropers_; idx +=5) {
    int a1 = atoi(infile_.NextElement()) - 1;
    int a2 = atoi(infile_.NextElement()) - 1;
    int a3 = atoi(infile_.NextElement()) - 1;
    int a4 = atoi(infile_.NextElement()) - 1;
    int didx = atoi(infile_.NextElement()) - 1;
    TopIn.SetChamber().AddImproperTerm( DihedralType(a1, a2, a3, a4, didx) );
  }
  return 0;
}

// Parm_Amber::ReadChamberImpFC()
int Parm_Amber::ReadChamberImpFC(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_CHM_IMPFC, N_impTerms_, FMT)) return 1;
  for (int idx = 0; idx != N_impTerms_; idx++)
    TopIn.SetChamber().SetImproperParm(idx).SetPk( atof(infile_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadChamberImpPHASE()
int Parm_Amber::ReadChamberImpPHASE(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_CHM_IMPP, N_impTerms_, FMT)) return 1;
  for (int idx = 0; idx != N_impTerms_; idx++)
    TopIn.SetChamber().SetImproperParm(idx).SetPhase( atof(infile_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadChamberLJ14A()
int Parm_Amber::ReadChamberLJ14A(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_LJ14A, numLJparm_, FMT)) return 1;
  for (int idx = 0; idx != numLJparm_; idx++)
    TopIn.SetChamber().SetLJ14(idx).SetA( atof(infile_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadChamberLJ14B()
int Parm_Amber::ReadChamberLJ14B(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_LJ14B, numLJparm_, FMT)) return 1;
  for (int idx = 0; idx != numLJparm_; idx++)
    TopIn.SetChamber().SetLJ14(idx).SetB( atof(infile_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadChamberCmapCounts()
int Parm_Amber::ReadChamberCmapCounts(FortranData const& FMT) {
  if (SetupBuffer(F_CHM_CMAPC, 2, FMT)) return 1;
  n_cmap_terms_ = atoi( infile_.NextElement() );
  n_cmap_grids_ = atoi( infile_.NextElement() );
  return 0;
}

// Parm_Amber::ReadChamberCmapRes()
/** Get CMAP resolutions for each grid and allocate grids. */
int Parm_Amber::ReadChamberCmapRes(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_CHM_CMAPR, n_cmap_grids_, FMT)) return 1;
  for (int i = 0; i != n_cmap_grids_; i++)
    TopIn.SetChamber().AddCmapGrid( CmapGridType( atoi(infile_.NextElement()) ) );
  return 0;
}

// Parm_Amber::ReadChamberCmapGrid()
/** Read CMAP grid. */
int Parm_Amber::ReadChamberCmapGrid(const char* CmapFlag, Topology& TopIn, FortranData const& FMT)
{
  // Figure out which grid this is.
  // 012345678901234567890123
  // CHARMM_CMAP_PARAMETER_XX
  int gridnum = convertToInteger( std::string( CmapFlag+22 ) ) - 1;
  if (gridnum < 0 || gridnum >= (int)TopIn.Chamber().CmapGrid().size()) {
    mprinterr("Error: CMAP grid '%s' out of range.\n", CmapFlag);
    return 1;
  }
  CmapGridType& GRID = TopIn.SetChamber().SetCmapGrid( gridnum );
  if (SetupBuffer(F_CHM_CMAPP, GRID.Size(), FMT)) return 1;
  for (int idx = 0; idx != GRID.Size(); idx++)
    GRID.SetGridPt( idx, atof(infile_.NextElement()) );
  return 0;
}

// Parm_Amber::ReadChamberCmapTerms()
int Parm_Amber::ReadChamberCmapTerms(Topology& TopIn, FortranData const& FMT) {
  int nvals = n_cmap_terms_ * 6;
  if (SetupBuffer(F_CHM_CMAPI, nvals, FMT)) return 1;
  for (int idx = 0; idx != nvals; idx += 6) {
    int a1 = atoi(infile_.NextElement()) - 1;
    int a2 = atoi(infile_.NextElement()) - 1;
    int a3 = atoi(infile_.NextElement()) - 1;
    int a4 = atoi(infile_.NextElement()) - 1;
    int a5 = atoi(infile_.NextElement()) - 1;
    int gidx = atoi(infile_.NextElement()) - 1;
    TopIn.SetChamber().AddCmapTerm( CmapType(a1, a2, a3, a4, a5, gidx) );
  }
  return 0;
}

// ----- LES -------------------------------------
int Parm_Amber::ReadLESntyp(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_LES_NTYP, 1, FMT)) return 1;
  nlestyp_ = atoi(infile_.NextElement());
  // This allocates the FAC array (nlestyp*nlestyp) and reserves LES array (natoms)
  TopIn.SetLES().Allocate( values_[NATOM], nlestyp_ );
  return 0;
}

int Parm_Amber::ReadLESfac(Topology& TopIn, FortranData const& FMT) {
  int nvals = nlestyp_ * nlestyp_;
  if (SetupBuffer(F_LES_FAC, nvals, FMT)) return 1;
  for (int idx = 0; idx != nvals; idx++)
    TopIn.SetLES().SetFAC( idx, atof(infile_.NextElement()) );
  return 0;
}

int Parm_Amber::ReadLEStypes(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_LES_TYPE, values_[NATOM], FMT)) return 1;
  for (int idx = 0; idx != values_[NATOM]; idx++)
    TopIn.SetLES().SetType( idx, atoi(infile_.NextElement()) );
  return 0;
}

int Parm_Amber::ReadLEScnum(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_LES_CNUM, values_[NATOM], FMT)) return 1;
  for (int idx = 0; idx != values_[NATOM]; idx++)
    TopIn.SetLES().SetCopy( idx, atoi(infile_.NextElement()) );
  return 0;
}

int Parm_Amber::ReadLESid(Topology& TopIn, FortranData const& FMT) {
  if (SetupBuffer(F_LES_ID, values_[NATOM], FMT)) return 1;
  for (int idx = 0; idx != values_[NATOM]; idx++)
    TopIn.SetLES().SetID( idx, atoi(infile_.NextElement()) );
  return 0;
}

// -----------------------------------------------------------------------------
void Parm_Amber::WriteHelp() {
  mprintf("\tnochamber: Do not write CHAMBER information to topology (useful for e.g. using"
          " topology for visualization with VMD).\n");
}

int Parm_Amber::processWriteArgs(ArgList& argIn) {
  nochamber_ = argIn.hasKey("nochamber");
  return 0;
}

int Parm_Amber::WriteParm(FileName const& fname, Topology const& parmIn) {
  return 1;
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
  if (ptrIn == 0) return 1;
  std::string fformat( NoTrailingWhitespace( ptrIn ) );
  if ( fformat.empty() ) return 1;
  //mprintf("DEBUG: Fortran format: %s\n", fformat.c_str());
  // Make sure characters are upper case.
  for (std::string::iterator p = fformat.begin(); p != fformat.end(); p++)
    toupper(*p);
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


