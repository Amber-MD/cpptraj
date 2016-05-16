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
  NTOTRS;   total number of residues
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
  ptype_(OLDPARM)
{}

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

// Parm_Amber::ReadParm()
int Parm_Amber::ReadParm(FileName const& fname, Topology &TopIn ) {
  if (infile_.OpenRead( fname )) return 1;
  if (ptype_ == OLDPARM)
    return ReadOldParm( TopIn );
  else
    return ReadNewParm( TopIn );
}

// Parm_Amber::ReadOldParm()
int Parm_Amber::ReadOldParm(Topology& TopIn) {
  mprintf("\tReading old (<v7) Amber Topology file.\n");
  std::string title = NoTrailingWhitespace( infile_.GetLine() );
  int Npointers = 30; // No NEXTRA etc
  if ( ReadPointers( Npointers, FortranData(FINT, 12, 6, 0) ) ) return 1;

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
    if ( ptr[0] == '%' ) {
      if (ptr[1] == 'V' && ptr[2] == 'E' && ptr[3] == 'R') {
        // %VERSION line. Skip it.
        //if (debug_ > 0)
        mprintf("DEBUG: Version: %s\n", NoTrailingWhitespace(ptr).c_str());
      } else if (IsFLAG(ptr)) {
        // %FLAG <type> line. Determine the flag type.
        std::string flagType = NoTrailingWhitespace(ptr+6);
        int flagIdx = -1;
        for (ParmPtr P = FLAGS_; P->Flag != 0; ++P) {
          if (flagType.compare(P->Flag) == 0) {
            flagIdx = (int)(P - FLAGS_);
            break;
          }
        }
        // Read the format line. Do this even if FLAG is not recognized.
        if (ReadFormatLine(FMT)) return 1;
        // Process the FLAG
        if (flagIdx == -1) {
          mprintf("Warning: Amber topology flag '%s' is unrecognized and will be skipped.\n",
                  flagType);
        } else {
          int err = 0;
          switch ((AmberParmFlagType)flagIdx) {
            case F_CTITLE: ptype_ = CHAMBER; // Fall through to F_TITLE
            case F_TITLE: ReadTitle(TopIn); break;
            case F_POINTERS: err = ReadPointers(AMBERPOINTERS_, FMT); break;
            default: return 1; // SANITY CHECK
          }
          if (err != 0) return 1;
        }
      }  
    }
  }

  infile_.CloseFile();
  return 0;
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
void Parm_Amber::ReadTitle(Topology& TopIn) {
  std::string title = NoTrailingWhitespace( infile_.GetLine() );
  if (debug_>0) mprintf("\tAmberParm Title: \"%s\"\n",title.c_str());
  TopIn.SetParmName( title, infile_.Filename() );
}

int Parm_Amber::ReadPointers(int Npointers, FortranData const& FMT) {
  infile_.SetupFrameBuffer( Npointers, FMT.Width(), FMT.Ncols() );
  
  return 0;
}

// -----------------------------------------------------------------------------
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
    mprintf("[%s]: cols=%i type=%i width=%i precision=%i\n",fformat.c_str(),
            fncols_,(int)ftype_,fwidth_,fprecision_);

  return 0;
}


