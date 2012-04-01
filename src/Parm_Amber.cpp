#include <cstdio> // sscanf
#include <cstring> // strlen, strncmp
#include <locale> // isspace
#include <cstdlib> // atoi, atof
#include "Parm_Amber.h"
#include "CpptrajStdio.h"
#include "Box.h"
#include "Constants.h" // ELECTOAMBER, AMBERTOELEC

// CONSTRUCTOR
Parm_Amber::Parm_Amber() : 
  newParm_(false),
  ftype_(UNKNOWN_FTYPE),
  fncols_(0),
  fprecision_(0),
  fwidth_(0),
  error_count_(0),
  buffer_(NULL)
{ }

// DESTRUCTOR
Parm_Amber::~Parm_Amber() {
  if (buffer_!=NULL) delete[] buffer_;
}

bool Parm_Amber::ID_ParmFormat() {
  int iamber[12];
  // Assumes already set up for READ
  if (OpenFile()) return false;
  IO->Gets(lineBuffer_, BUF_SIZE);
  // Check for %VERSION
  if (strncmp(lineBuffer_,"%VERSION",8)==0) {
    IO->Gets(lineBuffer_, BUF_SIZE);
    // Check for %FLAG
    if (strncmp(lineBuffer_,"%FLAG",5)==0) {
      if (debug_>0) mprintf("  AMBER TOPOLOGY file\n");
      newParm_ = true;
      CloseFile();
      return true;
    }
  } else {
    // Since %VERSION not present, If first line is 81 bytes and the second 
    // line has 12 numbers in 12I6 format, assume old-style Amber topology
    // NOTE: Could also be less than 81? Only look for 12 numbers?
    int line1size = (int)strlen(lineBuffer_);
    if (line1size == (81 + isDos_)) {
      IO->Gets(lineBuffer_, BUF_SIZE);
      if ( sscanf(lineBuffer_,"%6i%6i%6i%6i%6i%6i%6i%6i%6i%6i%6i%6i", 
                  iamber,   iamber+1, iamber+2,  iamber+3, 
                  iamber+4, iamber+5, iamber+6,  iamber+7, 
                  iamber+8, iamber+9, iamber+10, iamber+11) == 12 )
      {
        if (debug_>0) mprintf("  AMBER TOPOLOGY, OLD FORMAT\n");
        CloseFile();
        return true;
      }
    }
  }
  CloseFile();
  return false;
}

int Parm_Amber::ReadParm( Topology &TopIn ) {
  int err;
  if (OpenFile()) return 1;
  if (!newParm_)
    err = ReadParmOldAmber( TopIn );
  else
    err = ReadParmAmber( TopIn );
  CloseFile();
  if (err != 0) return 1;

  return 0;
}

// ---------- Constants and Enumerated types -----------------------------------
/// Combined size of %FLAG and %FORMAT lines (81 * 2)
const size_t Parm_Amber::FFSIZE=162;
/// Enumerated type for FLAG_POINTERS section
const int Parm_Amber::AMBERPOINTERS=31;
enum topValues {
//0       1       2      3       4       5       6       7      8       9
  NATOM,  NTYPES, NBONH, MBONA,  NTHETH, MTHETA, NPHIH,  MPHIA, NHPARM, NPARM,
  NNB,    NRES,   NBONA, NTHETA, NPHIA,  NUMBND, NUMANG, NPTRA, NATYP,  NPHB,
  IFPERT, NBPER,  NGPER, NDPER,  MBPER,  MGPER,  MDPER,  IFBOX, NMXRS, IFCAP,
  NEXTRA
};
/// Number of unique amber parm FLAGs
const int Parm_Amber::NUMAMBERPARMFLAGS=44;
/// Constant strings for fortran formats corresponding to Amber parm flags
const char Parm_Amber::AmberParmFmt[NUMAMBERPARMFLAGS][16] = {
"%FORMAT(10I8)",   "%FORMAT(20a4)",   "%FORMAT(5E16.8)", "%FORMAT(5E16.8)", "%FORMAT(20a4)",
"%FORMAT(10I8)",   "%FORMAT(20a4)",   "%FORMAT(10I8)",   "%FORMAT(10I8)",   "%FORMAT(3I8)",
"%FORMAT(10I8)",   "%FORMAT(5E16.8)", "%FORMAT(10I8)",   "%FORMAT(10I8)",   "%FORMAT(10I8)",
"%FORMAT(5E16.8)", "%FORMAT(5E16.8)", "%FORMAT(10I8)",   "%FORMAT(5E16.8)", "%FORMAT(5E16.8)",
"%FORMAT(5E16.8)", "%FORMAT(5E16.8)", "%FORMAT(5E16.8)", "%FORMAT(5E16.8)", "%FORMAT(5E16.8)",
"%FORMAT(5E16.8)", "%FORMAT(5E16.8)", "%FORMAT(5E16.8)", "%FORMAT(5E16.8)", "%FORMAT(5E16.8)",
"%FORMAT(10I8)",   "%FORMAT(10I8)",   "%FORMAT(10I8)",   "%FORMAT(10I8)",   "%FORMAT(5E16.8)",
"%FORMAT(5E16.8)", "%FORMAT(5E16.8)", "%FORMAT(20a4)",   "%FORMAT(10I8)",   "%FORMAT(10I8)",
"%FORMAT(10I8)",   "%FORMAT(20a4)",   "%FORMAT(20a4)",   "%FORMAT(1a80)"
};
/// Constant strings for Amber parm flags
const char Parm_Amber::AmberParmFlag[NUMAMBERPARMFLAGS][27] = {
  "POINTERS",
  "ATOM_NAME",
  "CHARGE",
  "MASS",
  "RESIDUE_LABEL",
  "RESIDUE_POINTER",
  "AMBER_ATOM_TYPE",
  "BONDS_INC_HYDROGEN",
  "BONDS_WITHOUT_HYDROGEN",
  "SOLVENT_POINTERS",
  "ATOMS_PER_MOLECULE",
  "BOX_DIMENSIONS",
  "ATOM_TYPE_INDEX",
  "NUMBER_EXCLUDED_ATOMS",
  "NONBONDED_PARM_INDEX",
  "LENNARD_JONES_ACOEF",
  "LENNARD_JONES_BCOEF",
  "EXCLUDED_ATOMS_LIST",
  "RADII",
  "SCREEN",
  "BOND_FORCE_CONSTANT",
  "BOND_EQUIL_VALUE",
  "ANGLE_FORCE_CONSTANT",
  "ANGLE_EQUIL_VALUE",
  "DIHEDRAL_FORCE_CONSTANT",
  "DIHEDRAL_PERIODICITY",
  "DIHEDRAL_PHASE",
  "SCEE_SCALE_FACTOR",
  "SCNB_SCALE_FACTOR",
  "SOLTY",
  "ANGLES_INC_HYDROGEN",
  "ANGLES_WITHOUT_HYDROGEN",
  "DIHEDRALS_INC_HYDROGEN",
  "DIHEDRALS_WITHOUT_HYDROGEN",
  "HBOND_ACOEF",
  "HBOND_BCOEF",
  "HBCUT",
  "TREE_CHAIN_CLASSIFICATION",
  "JOIN_ARRAY",
  "IROTAT",
  "ATOMIC_NUMBER",
  "TITLE",
  "CTITLE",
  "RADIUS_SET"
};
// -----------------------------------------------------------------------------
// Parm_Amber::ReadParmOldAmber()
int Parm_Amber::ReadParmOldAmber( Topology &TopIn ) {
  return 1;
}

// Parm_Amber::ReadParmAmber()
int Parm_Amber::ReadParmAmber( Topology &TopIn ) {
  std::vector<int> atomsPerMol;
  int finalSoluteRes = -1;
  int firstSolvMol = -1;
  Box parmbox;
  bool chamber; // True if this top is a chamber-created top file.
  std::string title;
  if (debug_>0) mprintf("Reading Amber Topology file %s\n",BaseName());

  // Title. If not found check for CTITLE (chamber)
  if (PositionFileAtFlag(F_TITLE)) {
    title = GetLine();
  } else {
    if (PositionFileAtFlag(F_CTITLE)) {
      title = GetLine();
      chamber = true;
    } else {
      // No TITLE or CTITLE, weird, but dont crash out yet.
      mprintf("Warning: [%s] No TITLE in Amber Parm.\n",BaseName());
    }
  }
  mprintf("\tAmberParm Title: [%s]\n",title.c_str());
  if (!title.empty())
    TopIn.SetParmName( title );
  else
    TopIn.SetParmName( BaseName() );

  std::vector<int> values = GetFlagInteger(F_POINTERS, AMBERPOINTERS);
  if (values.empty()) {
    mprinterr("Error: [%s] Could not get POINTERS from Amber Topology.\n",BaseName());
    return 1;
  }
  //for (std::vector<int>::iterator v = values.begin(); v != values.end(); v++)
  //  mprintf("\t%i\n",*v);

  // Read parm variables
  std::vector<NameType> names = GetFlagName(F_NAMES, values[NATOM]);
  std::vector<double> charge = GetFlagDouble(F_CHARGE,values[NATOM]);
  std::vector<int> at_num = GetFlagInteger(F_ATOMICNUM,values[NATOM]);
  std::vector<double> mass = GetFlagDouble(F_MASS,values[NATOM]);
  std::vector<int> atype_index = GetFlagInteger(F_ATYPEIDX,values[NATOM]);
  //parmOut.numex = GetFlagInteger(F_NUMEX,values[NATOM]);
  //parmOut.NB_index = GetFlagInteger(F_NB_INDEX,values[NTYPES]*values[NTYPES]);
  std::vector<NameType> resnames = GetFlagName(F_RESNAMES,values[NRES]);
  std::vector<int> resnums = GetFlagInteger(F_RESNUMS,values[NRES]);
  /*parmOut.bond_rk = GetFlagDouble(F_BONDRK, values[NUMBND]);
  parmOut.bond_req = GetFlagDouble(F_BONDREQ, values[NUMBND]);
  parmOut.angle_tk = GetFlagDouble(F_ANGLETK, values[NUMANG]);
  parmOut.angle_teq = GetFlagDouble(F_ANGLETEQ, values[NUMANG]);
  parmOut.dihedral_pk = GetFlagDouble(F_DIHPK, values[NPTRA]);
  parmOut.dihedral_pn = GetFlagDouble(F_DIHPN, values[NPTRA]);
  parmOut.dihedral_phase = GetFlagDouble(F_DIHPHASE, values[NPTRA]);
  parmOut.scee_scale = GetFlagDouble(F_SCEE,values[NPTRA]);
  parmOut.scnb_scale = GetFlagDouble(F_SCNB,values[NPTRA]);
  // SOLTY: currently unused
  parmOut.solty = GetFlagDouble(F_SOLTY, values[NATYP]);
  int nlj = values[NTYPES] * (values[NTYPES]+1) / 2;
  parmOut.LJ_A = GetFlagDouble(F_LJ_A,nlj);
  parmOut.LJ_B = GetFlagDouble(F_LJ_B,nlj);*/
  std::vector<int> bondsh = GetFlagInteger(F_BONDSH,values[NBONH]*3);
  std::vector<int> bonds = GetFlagInteger(F_BONDS,values[NBONA]*3);
  /*parmOut.anglesh = GetFlagInteger(F_ANGLESH, values[NTHETH]*4);
  parmOut.angles = GetFlagInteger(F_ANGLES, values[NTHETA]*4);
  parmOut.dihedralsh = GetFlagInteger(F_DIHH, values[NPHIH]*5);
  parmOut.dihedrals = GetFlagInteger(F_DIH, values[NPHIA]*5);
  parmOut.excludedAtoms = GetFlagInteger(F_EXCLUDE, values[NNB]);
  parmOut.asol = GetFlagDouble(F_ASOL,values[NPHB]);
  parmOut.bsol = GetFlagDouble(F_BSOL,values[NPHB]);
  parmOut.hbcut = GetFlagDouble(F_HBCUT,values[NPHB]);*/
  std::vector<NameType> types = GetFlagName(F_TYPES,values[NATOM]);
  /*parmOut.itree = GetFlagName(F_ITREE,values[NATOM]);
  parmOut.join_array = GetFlagInteger(F_JOIN,values[NATOM]);
  parmOut.irotat = GetFlagInteger(F_IROTAT,values[NATOM]);*/
  // Get solvent info if IFBOX>0
  if (values[IFBOX]>0) {
    std::vector<int> solvent_pointer = GetFlagInteger(F_SOLVENT_POINTER,3);
    if (solvent_pointer.empty()) {
      mprintf("Could not get %s from Amber Topology file.\n",AmberParmFlag[F_SOLVENT_POINTER]);
      return 1;
    }
    finalSoluteRes = solvent_pointer[0] - 1;
    int molecules = solvent_pointer[1];
    firstSolvMol = solvent_pointer[2] - 1;
    atomsPerMol = GetFlagInteger(F_ATOMSPERMOL,molecules);
    // boxFromParm = {OLDBETA, BOX(1), BOX(2), BOX(3)}
    std::vector<double> boxFromParm = GetFlagDouble(F_PARMBOX,4);
    // If no box information present in the parm (such as with Chamber prmtops)
    // set the box info if ifbox = 2, otherwise default is NOBOX; the box info 
    // will eventually be set by angles from the first trajectory associated 
    // with this parm.
    if (boxFromParm.empty()) {
      if (not chamber) mprintf("Warning: Prmtop missing Box information.\n");
      // ifbox 2: truncated octahedron for certain
      if (values[IFBOX] == 2) 
        parmbox.SetTruncOct(); 
    // Determine box type, set Box angles and lengths from beta (boxFromParm[0])
    } else {
      parmbox.SetBetaLengths( boxFromParm );
    }
  }
  // GB parameters; radius set, radii, and screening parameters
  if (PositionFileAtFlag(F_RADSET)) {
    std::string radius_set = GetLine();
    mprintf("\tRadius Set: %s\n",radius_set.c_str());
  }
  std::vector<double> gb_radii = GetFlagDouble(F_RADII,values[NATOM]);
  std::vector<double> gb_screen = GetFlagDouble(F_SCREEN,values[NATOM]); 

  if (error_count_==0) {
    // Convert Amber Charge to elec
    for (std::vector<double>::iterator q = charge.begin(); q != charge.end(); q++)
      *q *= (AMBERTOELEC);
    // Shift atom #s in resnums by -1 so they start from 0
    for (std::vector<int>::iterator r = resnums.begin(); r != resnums.end(); r++)
      --(*r);
    // Shift atom #s in excludedAtoms by -1 so they start from 0 
    error_count_ += TopIn.CreateAtomArray( names, charge, at_num, mass, atype_index, 
                                           types, gb_radii, gb_screen,
                                           resnames, resnums );
    error_count_ += TopIn.SetBondInfo(bonds, bondsh);
    if (values[IFBOX]>0)
      error_count_ += TopIn.CreateMoleculeArray(atomsPerMol, parmbox, 
                                                finalSoluteRes, firstSolvMol);
  }

  return error_count_;
}

// Parm_Amber::GetFlagLine()
// NOTE: Useful?
std::string Parm_Amber::GetFlagLine(AmberParmFlagType flag) {
  std::string line;
  if (!PositionFileAtFlag(flag)) return line;
  return GetLine();
}

// Parm_Amber::GetLine()
std::string Parm_Amber::GetLine() {
  std::locale loc;

  IO->Gets(lineBuffer_, BUF_SIZE);
  std::string line( lineBuffer_ );
  // Remove any trailing whitespace
  std::string::iterator p = line.end();
  --p;
  for (; p != line.begin(); p--) 
    if (!isspace( *p, loc)) break;
  size_t lastSpace = (size_t)(p - line.begin()) + 1;
  //mprintf("lastSpace = %zu\n",lastSpace);
  if (lastSpace==1)
    line.clear();
  else
    line.resize( lastSpace );
  return line;
}

// Parm_Amber::GetInteger()
std::vector<int> Parm_Amber::GetInteger(int width, int ncols, int maxval) {
  std::vector<int> iarray;
  // Read prmtop section into buffer
  int err = AllocateAndRead( width, ncols, maxval );
  if (err == 0)
    return iarray;
  else if (err == -1) {
    mprinterr("Error in read of integer values from %s\n",BaseName());
    ++error_count_;
    return iarray;
  }
  // Reserve variable memory
  iarray.reserve( maxval );
  // Convert values in buffer to integer 
  char *ptrbegin = buffer_;
  char *ptrend = buffer_;
  for (int i = 0; i < maxval; i++) {
    // Advance past newlines / CR (dos)
    while (*ptrbegin=='\n' || *ptrbegin=='\r') ++ptrbegin;
    ptrend = ptrbegin + width;
    char lastchar = *ptrend;
    *ptrend = '\0';
    iarray.push_back( atoi(ptrbegin) );
    *ptrend = lastchar;
    ptrbegin = ptrend;
  }
  return iarray;
}

// Parm_Amber::GetDouble()
std::vector<double> Parm_Amber::GetDouble(int width, int ncols, int maxval) {
  std::vector<double> darray;
  // Read prmtop section into buffer
  int err = AllocateAndRead( width, ncols, maxval );
  if (err == 0)
    return darray;
  else if (err == -1) {
    mprinterr("Error in read of double values from %s\n",BaseName());
    ++error_count_;
    return darray;
  }
  // Reserve variable memory
  darray.reserve( maxval );
  // Convert values in buffer to integer 
  char *ptrbegin = buffer_;
  char *ptrend = buffer_;
  for (int i = 0; i < maxval; i++) {
    // Advance past newlines / CR (dos)
    while (*ptrbegin=='\n' || *ptrbegin=='\r') ++ptrbegin;
    ptrend = ptrbegin + width;
    char lastchar = *ptrend;
    *ptrend = '\0';
    darray.push_back( atof(ptrbegin) );
    *ptrend = lastchar;
    ptrbegin = ptrend;
  }
  return darray;
}

// Parm_Amber::GetName()
std::vector<NameType> Parm_Amber::GetName(int width, int ncols, int maxval) {
  std::vector<NameType> carray;
  // Read prmtop section into buffer
  int err = AllocateAndRead( width, ncols, maxval );
  if (err == 0)
    return carray;
  else if (err == -1) {
    mprinterr("Error in read of Name values from %s\n",BaseName());
    ++error_count_;
    return carray;
  }
  // Reserve variable memory
  carray.reserve( maxval );
  // Convert values in buffer to integer 
  char *ptrbegin = buffer_;
  char *ptrend = buffer_;
  for (int i = 0; i < maxval; i++) {
    // Advance past newlines / CR (dos)
    while (*ptrbegin=='\n' || *ptrbegin=='\r') ++ptrbegin;
    ptrend = ptrbegin + width;
    char lastchar = *ptrend;
    *ptrend = '\0';
    carray.push_back( ptrbegin );
    *ptrend = lastchar;
    ptrbegin = ptrend;
  }
  return carray;
}

// Parm_Amber::GetFlagInteger()
std::vector<int> Parm_Amber::GetFlagInteger(AmberParmFlagType fflag, int maxval) {
  std::vector<int> iarray;
  // Seek to flag and set up fncols, fwidth
  if (!SeekToFlag(fflag)) return iarray;
  // NOTE: Check that type matches?
  iarray = GetInteger( fwidth_, fncols_, maxval );
  return iarray;
}

// Parm_Amber::GetFlagDouble()
std::vector<double> Parm_Amber::GetFlagDouble(AmberParmFlagType fflag, int maxval) {
  std::vector<double> darray;
  // Seek to flag and set up fncols, fwidth
  if (!SeekToFlag(fflag)) return darray;
  // NOTE: Check that type matches?
  darray = GetDouble( fwidth_, fncols_, maxval );
  return darray;
}

// Parm_Amber::GetFlagName()
std::vector<NameType> Parm_Amber::GetFlagName(AmberParmFlagType fflag, int maxval) {
  std::vector<NameType> carray;
  // Seek to flag and set up fncols, fwidth
  if (!SeekToFlag(fflag)) return carray;
  // NOTE: Check that type matches?
  carray = GetName( fwidth_, fncols_, maxval );
  return carray;
}

bool Parm_Amber::SeekToFlag(AmberParmFlagType fflag) {
  // Find flag, set format line
  if (!PositionFileAtFlag(fflag)) return false;
  // Set up cols, width, etc from format
  if (!SetFortranType()) return false;
  return true;
}  

// Parm_Amber::GetFortranBufferSize()
/** Given number of columns and the width of each column, return the 
  * necessary char buffer size for N data elements.
  */
size_t Parm_Amber::GetFortranBufferSize(int width, int ncols, int N) {
  size_t bufferLines=0;
  size_t BufferSize=0;

  BufferSize = N * width;
  bufferLines = N / ncols;
  if ((N % ncols)!=0) ++bufferLines;
  // If DOS file there are CRs before Newlines
  if (isDos_) bufferLines *= 2;
  BufferSize += bufferLines;
  //if (debug_>0) 
  //  fprintf(stdout,"*** Buffer size is %i including %i newlines.\n",BufferSize,bufferLines);
  return BufferSize;
}

// Parm_Amber::AllocateAndRead()
int Parm_Amber::AllocateAndRead(int width, int ncols, int maxval) {
  char temp[3]; // Only for when maxval is 0, space for \n, \r, NULL
  // If # expected values is 0 there will still be a newline placeholder
  // in the parmtop. Read past that and return
  if (maxval==0) {
    IO->Gets(temp,2);
    return 0;
  }
  // Allocate buffer to read in entire section
  size_t BufferSize = GetFortranBufferSize(width, ncols, maxval);
  if (buffer_!=NULL) delete[] buffer_;
  buffer_ = new char[ BufferSize ];
  // Read section from file
  int err = IO->Read(buffer_,BufferSize);
  return err;
}

// Parm_Amber::PositionFileAtFlag()
bool Parm_Amber::PositionFileAtFlag(AmberParmFlagType flag) {
  const char *Key = AmberParmFlag[flag];
  char value[BUF_SIZE];
  bool searchFile = true;
  bool hasLooped = false;

  if (debug_ > 0) mprintf("Reading %s\n",Key);
  // Search for '%FLAG <Key>'
  while ( searchFile ) {
    while ( IO->Gets(lineBuffer_, BUF_SIZE) == 0 ) {
      if ( strncmp(lineBuffer_,"%FLAG",5)==0 ) {
        // Check flag Key
        sscanf(lineBuffer_, "%*s %s",value);
        if (strcmp(value,Key)==0) {
          if (debug_>1) mprintf("DEBUG: Found Flag Key [%s]\n",value);
          // Read next line; can be either a COMMENT or FORMAT. If COMMENT, 
          // read past until you get to the FORMAT line
          IO->Gets(lineBuffer_, BUF_SIZE); 
          while (strncmp(lineBuffer_, "%FORMAT",7)!=0)
            IO->Gets(lineBuffer_, BUF_SIZE);
          if (debug_>1) mprintf("DEBUG: Format line [%s]\n",lineBuffer_);
          // Set format
          // NOTE: Check against stored formats?
          fformat_.assign(lineBuffer_);
          return true;
        } // END found Key
      } // END found FLAG line
    } // END scan through file
    // If we havent yet tried to search from the beginning, try it now.
    // Otherwise the Key has not been found.
    if (!hasLooped) {
      IO->Rewind();
      hasLooped = true;
    } else
      searchFile = false;
  }

  // If we reached here Key was not found.
  if (debug_>0)
    mprintf("Warning: [%s] Could not find Key %s in file.\n",BaseName(),Key);
  fformat_.clear();
  return false;
}

// Parm_Amber::SetFortranType()
/** Given a fortran-type format string, set the corresponding fortran
  * type. Set fncols (if present), fwidth, and fprecision (if present).
  */
// 01234567
// %FORMAT([<cols>][(]<type><width>[<precision>][)])
bool Parm_Amber::SetFortranType() {
  std::locale loc;
  std::string arg;

  if ( fformat_.empty() ) return false;
  // Make sure characters are upper case.
  for (std::string::iterator p = fformat_.begin(); p != fformat_.end(); p++)
    toupper(*p, loc);
  // Advance past left parentheses
  std::string::iterator ptr = fformat_.begin() + 7;
  while (*ptr=='(') ++ptr;
  // If digit, have number of data columns
  fncols_ = 0;
  if (isdigit(*ptr, loc)) {
    while (ptr!=fformat_.end() && isdigit(*ptr, loc)) {
      arg += *ptr;
      ++ptr;
    }
    fncols_ = atoi( arg.c_str() );
  }
  // Advance past any more left parentheses
  while (ptr!=fformat_.end() && *ptr=='(') ++ptr;
  // Type
  if (ptr==fformat_.end()) {
    mprinterr("Error: Malformed fortran format string (%s) in Amber Topology %s\n",
              fformat_.c_str(), BaseName());
    return false;
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
  while (isdigit(*ptr,loc)) {
    arg += *ptr;
    ++ptr;
  }
  fwidth_ = atoi( arg.c_str() );
  // Precision
  fprecision_ = 0;
  if (*ptr == '.') {
    ++ptr;
    arg.clear();
    while (isdigit(*ptr,loc)) {
      arg += *ptr;
      ++ptr;
    }
    fprecision_ = atoi( arg.c_str() );
  }
  if (debug_ > 2)
    mprintf("[%s]: cols=%i type=%i width=%i precision=%i\n",fformat_.c_str(),
            fncols_,(int)ftype_,fwidth_,fprecision_);

  return true;
}
