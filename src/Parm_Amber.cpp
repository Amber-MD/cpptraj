// Parm_Amber.cpp
// NOTE: Eventually no memcpy!
#include <cstring> // memcpy, strcpy
#include <cstdlib> // atoi, atof
#include <cstdio>  // sscanf
#include <ctime>   // for writing time/date to amber parmtop (necessary?)
#include "Parm_Amber.h"
#include "CpptrajStdio.h"
#include "Constants.h" // ELECTOAMBER
// For file write
#include "CharBuffer.h"
// Compiler Defines:
// - USE_CHARBUFFER: Use CharBuffer to buffer entire file

// ---------- Constants and Enumerated types -----------------------------------
/// Combined size of %FLAG and %FORMAT lines (81 * 2)
const size_t AmberParmFile::FFSIZE=162;
/// Enumerated type for FLAG_POINTERS section
const int AmberParmFile::AMBERPOINTERS=31;
enum topValues {
//0       1       2      3       4       5       6       7      8       9
  NATOM,  NTYPES, NBONH, MBONA,  NTHETH, MTHETA, NPHIH,  MPHIA, NHPARM, NPARM,
  NNB,    NRES,   NBONA, NTHETA, NPHIA,  NUMBND, NUMANG, NPTRA, NATYP,  NPHB,
  IFPERT, NBPER,  NGPER, NDPER,  MBPER,  MGPER,  MDPER,  IFBOX, NMXRS, IFCAP,
  NEXTRA
};
/// Number of unique amber parm FLAGs
const int AmberParmFile::NUMAMBERPARMFLAGS=41;
/// Constant strings for fortran formats corresponding to Amber parm flags
const char AmberParmFile::AmberParmFmt[AmberParmFile::NUMAMBERPARMFLAGS][16] = {
"%FORMAT(10I8)",   "%FORMAT(20a4)",   "%FORMAT(5E16.8)", "%FORMAT(5E16.8)", "%FORMAT(20a4)",
"%FORMAT(10I8)",   "%FORMAT(20a4)",   "%FORMAT(10I8)",   "%FORMAT(10I8)",   "%FORMAT(3I8)",
"%FORMAT(10I8)",   "%FORMAT(5E16.8)", "%FORMAT(10I8)",   "%FORMAT(10I8)",   "%FORMAT(10I8)",
"%FORMAT(5E16.8)", "%FORMAT(5E16.8)", "%FORMAT(10I8)",   "%FORMAT(5E16.8)", "%FORMAT(5E16.8)",
"%FORMAT(5E16.8)", "%FORMAT(5E16.8)", "%FORMAT(5E16.8)", "%FORMAT(5E16.8)", "%FORMAT(5E16.8)",
"%FORMAT(5E16.8)", "%FORMAT(5E16.8)", "%FORMAT(5E16.8)", "%FORMAT(5E16.8)", "%FORMAT(5E16.8)",
"%FORMAT(10I8)",   "%FORMAT(10I8)",   "%FORMAT(10I8)",   "%FORMAT(10I8)",   "%FORMAT(5E16.8)",
"%FORMAT(5E16.8)", "%FORMAT(5E16.8)", "%FORMAT(20a4)",   "%FORMAT(10I8)",   "%FORMAT(10I8)",
"%FORMAT(10I8)"
};
/// Constant strings for Amber parm flags
const char AmberParmFile::AmberParmFlag[NUMAMBERPARMFLAGS][27] = {
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
  "ATOMIC_NUMBER"
};
// -----------------------------------------------------------------------------

// ---------- INTERNAL FUNCTIONS -----------------------------------------------
// GetFortranBufferSize()
/** Given number of columns and the width of each column, return the 
  * necessary char buffer size for N data elements.
  */
static size_t GetFortranBufferSize(int N, int isDos, int width, int ncols) {
  size_t bufferLines=0;
  size_t BufferSize=0;

  BufferSize = N * width;
  bufferLines = N / ncols;
  if ((N % ncols)!=0) ++bufferLines;
  // If DOS file there are CRs before Newlines
  if (isDos) bufferLines *= 2;
  BufferSize += bufferLines;
  //if (debug>0) 
  //  fprintf(stdout,"*** Buffer size is %i including %i newlines.\n",BufferSize,bufferLines);
  return BufferSize;
}

// PositionFileAtFlag()
/// Position the given file at the given flag. Set the corresponding format.
/** If fformat is NULL just report whether flag was found. This routine will
  * first attempt to search for the flag from the current File position. If
  * the flag is not found at first the routine will rewind and attempt to
  * search from the beginning. This can speed up reading of parm files
  * when this routine is called for Keys that are in the same order as
  * in the input file.
  * \return true if flag was found, false if not.
  */
static bool PositionFileAtFlag(CpptrajFile &File, const char *Key, char *fformat, int debug) {
  char lineBuffer[BUFFER_SIZE]; // Hold flag/format line from parmfile
  char value[83];
  bool hasLooped = false;
  bool searchFile = true;

  if (debug>0) mprintf("Reading %s\n",Key);
  // First, rewind the input file.
  //File->IO->Rewind();
  // Search for %FLAG <Key>
  while ( searchFile ) {
#   ifdef USE_CHARBUFFER
    while ( File.Gets(lineBuffer,BUFFER_SIZE) == 0)
#   else
    while ( File.IO->Gets(lineBuffer,BUFFER_SIZE) == 0)
#   endif
    {
      if ( strncmp(lineBuffer,"%FLAG",5)==0 ) {
        sscanf(lineBuffer,"%*s %s",value);
        if (strcmp(value,Key)==0) {
          if (debug>1) mprintf("DEBUG: Found Flag Key [%s]\n",value);
          // Read next line; can be either a COMMENT or FORMAT. If COMMENT, 
          // read past until you get to the FORMAT line
#         ifdef USE_CHARBUFFER
          File.Gets(lineBuffer,BUFFER_SIZE);
#         else
          File.IO->Gets(lineBuffer,BUFFER_SIZE);
#         endif
          while (strncmp(lineBuffer,"%FORMAT",7)!=0)
#           ifdef USE_CHARBUFFER
            File.Gets(lineBuffer,BUFFER_SIZE);
#           else
            File.IO->Gets(lineBuffer,BUFFER_SIZE);
#           endif
          if (debug>1) mprintf("DEBUG: Format line [%s]\n",lineBuffer);
          // Set format
          if (fformat!=NULL) strcpy(fformat, lineBuffer);
          return true;
        } // END found Key
      } // END found FLAG line
    } // END scan through file
    // If we havent yet tried to search from the beginning, try it now.
    // Otherwise the Key has not been found.
    if (!hasLooped) {
#     ifdef USE_CHARBUFFER
      File.Rewind();
#     else
      File.IO->Rewind();
#     endif
      hasLooped = true;
    } else
      searchFile = false;
  }

  // If we reached here Key was not found.
  if (debug>0)
    mprintf("Warning: [%s] Could not find Key %s in file.\n",File.filename,Key);
  if (fformat!=NULL) strcpy(fformat,"");
  return false;
}

// RemoveWhitespace()
/// Remove terminal whitespace from string (including newline + CR) 
static void RemoveWhitespace(char *bufferIn) {
  char *ptr = NULL;
  char *end;
  if (bufferIn==NULL) return;
  // Position ptr at the last char of bufferIn
  end = bufferIn + strlen(bufferIn) - 1;
  for (ptr = end; ptr >= bufferIn; ptr--) {
    // Stop at first non-whitespace non-newline char.
    if (*ptr!=' ' && *ptr!='\n' && *ptr!='\r') break;
  }
  // Put a NULL just after the first non-whitespace char
  ++ptr;
  *ptr='\0';
}

// ---------------------------------------------------------

// AmberParmFile::GetFortranType()
/** Given a fortran-type format string, return the corresponding fortran
  * type. Set ncols (if present), width, and precision (if present).
  */
// 01234567
// %FORMAT([<cols>][(]<type><width>[<precision>][)])
AmberParmFile::FortranType AmberParmFile::GetFortranType(char *FormatIn, int *ncols, 
                                                 int *width, int *precision) 
{
  int idx;
  char temp[32];
  char Format[32];
  char *ptr;
  FortranType ftype;
  // Make sure characters are upper-case
  // NOTE: Maybe do some checking for parentheses etc here
  //       Not expecting format strings to be > 32
  ptr = FormatIn;
  idx = 0;
  while (*ptr!='\0') {
    Format[idx++] = toupper( *ptr );
    if (idx==32) break;
    ptr++;
  }
  // Advance past left parentheses
  ptr = Format + 7;
  while (*ptr=='(') ptr++;
  // If digit, have number of data columns
  *ncols = 0;
  if (isdigit(*ptr)) {
    idx = 0;
    while (isdigit(*ptr)) {temp[idx++] = *ptr; ptr++;}
    temp[idx]='\0';
    *ncols = atoi(temp);
  }
  // Advance past any more left parentheses
  while (*ptr=='(') ptr++;
  // Type
  switch (*ptr) {
    case 'I' : ftype = FINT;    break;
    case 'E' : ftype = FDOUBLE; break;
    case 'A' : ftype = FCHAR;   break;
    case 'F' : ftype = FFLOAT;  break;
    default  : ftype = UNKNOWN_FTYPE;
  }
  ptr++;
  // Width
  idx = 0;
  while (isdigit(*ptr)) {temp[idx++] = *ptr; ptr++;}
  temp[idx]='\0';
  *width = atoi(temp);
  // Precision
  *precision = 0;
  if (*ptr == '.') {
    ptr++;
    idx = 0;
    while (isdigit(*ptr)) {temp[idx++] = *ptr; ptr++;}
    temp[idx]='\0';
    *precision = atoi(temp);
  }
  //mprintf("[%s]: cols=%i type=%c width=%i precision=%i\n",Format,
  //        *ncols,(int)ftype,*width,*precision);

  return ftype;
}


// DataToFortranBuffer()
/** Write N data elements stored in I, D, or C to character buffer with given 
  * fortran format.
  */
int AmberParmFile::DataToFortranBuffer(CharBuffer &buffer, 
                                       AmberParmFile::AmberParmFlagType fFlag,
                                       int *I, double *D, NAME *C, int N)
{
  int coord, width, numCols, precision;
  FortranType fType;
  char FormatString[32];

  // Determine type, cols, width, and precision from format string
  fType = GetFortranType((char*)AmberParmFmt[fFlag], &numCols, &width, &precision);
  if (fType == UNKNOWN_FTYPE) {
    mprinterr("Error: DataToFortranBuffer: Unknown format string [%s]\n",AmberParmFmt[fFlag]);
    return 1;
  }

  // If called with N == 0, or all NULL, want the FLAG and FORMAT lines but 
  // no data, just a newline.
  if (N==0 || (I==NULL && D==NULL && C==NULL)) {
    buffer.IncreaseSize( FFSIZE + 1 ); // FLAG + FORMAT lines, + newline
    buffer.Sprintf("%%FLAG %-74s\n",AmberParmFlag[fFlag]);
    buffer.Sprintf("%-80s\n",AmberParmFmt[fFlag]);
    buffer.NewLine();
    return 0;
  }

  // Increase the buffer by the appropriate amount
  size_t delta = GetFortranBufferSize(N,0,width,numCols);
  delta += FFSIZE; // FFSIZE is Combined size of %FLAG and %FORMAT lines (81 * 2)
  buffer.IncreaseSize( delta );

  //fprintf(stdout,"*** Called DataToBuffer: N=%i, width=%i, numCols=%i, format=%s\n",
  //        N, width, numCols, FormatString);
  //fprintf(stdout,"*** Buffer address is %p\n",buffer);

  // Print FLAG and FORMAT lines
  // '%FLAG '
  //  012345
  buffer.Sprintf("%%FLAG %-74s\n",AmberParmFlag[fFlag]);
  //sprintf(ptr,"%-80s\n",FLAG);
  buffer.Sprintf("%-80s\n",AmberParmFmt[fFlag]);

  // Write Integer data
  // NOTE: move the coord+1 code?
  coord=0;
  if (fType==FINT) {
    if (I==NULL) {
      mprinterr("Error: DataToFortranBuffer: INT is NULL.\n");
      return 1;
    }
    sprintf(FormatString,"%%%ii",width);
    for (coord=0; coord < N; coord++) {
      buffer.WriteInteger(FormatString, I[coord]);
      if ( ((coord+1)%numCols)==0 ) buffer.NewLine();
    }

  // Write Double data
  } else if (fType==FDOUBLE) {
    if (D==NULL) {
      mprinterr("Error: DataToFortranBuffer: DOUBLE is NULL.\n");
      return 1;
    }
    sprintf(FormatString,"%%%i.%ilE",width,precision);
    for (coord=0; coord < N; coord++) {
      buffer.WriteDouble(FormatString,D[coord]);
      if ( ((coord+1)%numCols)==0 ) buffer.NewLine();
    }

  // Write Char data
  } else if (fType==FCHAR) {
    if (C==NULL) {
      mprinterr("Error: DataToFortranBuffer: CHAR is NULL.\n");
      return 1;
    }
    sprintf(FormatString,"%%%is",width);
    for (coord=0; coord < N; coord++) {
      buffer.WriteString(FormatString,C[coord]);
      if ( ((coord+1)%numCols)==0 ) buffer.NewLine();
    }

  // Write Float data
  } else if (fType==FFLOAT) {
    if (D==NULL) {
      mprinterr("Error: DataToFortranBuffer: FLOAT is NULL.\n");
      return 1;
    }
    sprintf(FormatString,"%%%i.%ilf",width,precision);
    for (coord=0; coord < N; coord++) {
      buffer.WriteDouble(FormatString,D[coord]);
      if ( ((coord+1)%numCols)==0 ) buffer.NewLine();
    }
  }

  // If the coord record didnt end on a newline, print one
  if ( (coord%numCols)!=0 ) buffer.NewLine();
  return 0;
}

// ---------- STL functions ---------------------------------------------------- 
// CONSTRUCTOR
AmberParmFile::AmberParmFile() {
  buffer = NULL;
  File = NULL;
  error_count = 0;
}

// DESTRUCTOR
AmberParmFile::~AmberParmFile() {
  if (buffer!=NULL) delete[] buffer;
}

// AmberParmFile::AllocateAndRead()
int AmberParmFile::AllocateAndRead(int &width, int &ncols, int &maxval) {
  char temp[3]; // Only for when maxval is 0, space for \n, \r, NULL
  int err;
  // If # expected values is 0 there will still be a newline placeholder
  // in the parmtop. Read past that and return
  if (maxval==0) {
#   ifdef USE_CHARBUFFER
    File->Gets(temp,2);
#   else
    File->IO->Gets(temp,2);
#   endif
    return 0;
  }
  // Allocate buffer to read in entire section
  size_t BufferSize = GetFortranBufferSize(maxval, File->isDos, width, ncols);
  if (buffer!=NULL) delete[] buffer;
  buffer = new char[ BufferSize ];
  // Read section from file
  # ifdef USE_CHARBUFFER
  err = File->Read(buffer,BufferSize);
# else
  err = File->IO->Read(buffer,sizeof(char),BufferSize);
# endif
  return err;
}

// AmberParmFile::GetDouble()
double *AmberParmFile::GetDouble(int width, int ncols, int maxval)
{
  double *D;
  int err;
  // Read prmtop section into buffer
  err = AllocateAndRead(width,ncols,maxval);
  if (err == 0)
    return NULL;
  else if ( err == -1) {
    mprinterr("Error in read of double values from %s\n",File->filename);
    ++error_count;
    return NULL;
  }
  // Reserve variable memory
  D = new double[ maxval ];
  // Convert values in buffer to double
  char *ptrbegin = buffer;
  char *ptrend = buffer;
  for (int i = 0; i < maxval; i++) {
    // Advance past newlines / CR (dos)
    while (*ptrbegin=='\n' || *ptrbegin=='\r') ++ptrbegin;
    ptrend = ptrbegin + width;
    char lastchar = *ptrend;
    *ptrend = '\0';
    D[i] = atof(ptrbegin);
    *ptrend = lastchar;
    ptrbegin = ptrend;
  }
  return D;
}

// AmberParmFile::GetInteger()
int *AmberParmFile::GetInteger(int width, int ncols, int maxval)
{
  int *I;
  int err;
  // Read prmtop section into buffer
  err = AllocateAndRead(width,ncols,maxval);
  if (err == 0)
    return NULL;
  else if ( err == -1) {
    mprinterr("Error in read of integer values from %s\n",File->filename);
    ++error_count;
    return NULL;
  }
  // Reserve variable memory
  I = new int[ maxval ]; 
  // Convert values in buffer to integer 
  char *ptrbegin = buffer;
  char *ptrend = buffer;
  for (int i = 0; i < maxval; i++) {
    // Advance past newlines / CR (dos)
    while (*ptrbegin=='\n' || *ptrbegin=='\r') ++ptrbegin;
    ptrend = ptrbegin + width;
    char lastchar = *ptrend;
    *ptrend = '\0';
    I[i] = atoi(ptrbegin);
    *ptrend = lastchar;
    ptrbegin = ptrend;
  }
  return I;
}

// AmberParmFile::GetName()
NAME *AmberParmFile::GetName(int width, int ncols, int maxval)
{
  NAME *C;
  int err;
  // Read prmtop section into buffer
  err = AllocateAndRead(width,ncols,maxval);
  if (err == 0)
    return NULL;
  else if ( err == -1) {
    mprinterr("Error in read of string values from %s\n",File->filename);
    ++error_count;
    return NULL;
  }
  // Reserve variable memory
  C = new NAME[ maxval ]; 
  // Convert values in buffer to NAME 
  char *ptrbegin = buffer;
  for (int i = 0; i < maxval; i++) {
    // Advance past newlines / CR (dos)
    while (*ptrbegin=='\n' || *ptrbegin=='\r') ++ptrbegin;
    // Copy width characters
    for (int j = 0; j < width; j++) {
      C[i][j] = *ptrbegin;
      ++ptrbegin;
    }
    C[i][width]='\0';
  }
  return C;
}

// AmberParmFile::GetLine()
char *AmberParmFile::GetLine() {
  char *lineBuffer = new char[83]; // 80 + newline + NULL ( + CR if dos)
# ifdef USE_CHARBUFFER
  File->Gets(lineBuffer,82);
# else
  File->IO->Gets(lineBuffer,82);
# endif
  RemoveWhitespace(lineBuffer);
  return lineBuffer;
}

// AmberParmFile::GetFlagLine()
char *AmberParmFile::GetFlagLine(const char* Key) {
  char fformat[83];
  // Get Flag Key
  // Find flag, not concerned with format
  if (!PositionFileAtFlag(*File,Key,fformat,debug)) return NULL;
  if (debug>1) 
    mprintf("DEBUG: Flag line [%s]\n",Key);
  return GetLine();
}

// AmberParmFile::SeekToFlag()
AmberParmFile::FortranType AmberParmFile::SeekToFlag(AmberParmFlagType fflag, 
                                                    int &ncols, int &width, int &precision)
{
  char fformat[83]; // Hold FORMAT line
  FortranType fType;
  // Get Flag Key
  char *Key = (char*) AmberParmFlag[fflag];
  // Find flag, get format
  if (!PositionFileAtFlag(*File,Key,fformat,debug)) return UNKNOWN_FTYPE;
  // Determine cols, width etc from format
  fType = GetFortranType(fformat, &ncols, &width, &precision);
  if (debug>1) 
    mprintf("DEBUG: Flag [%s] Type %i Format[%s]\n",Key,(int)fType, fformat);
  return fType; 
}

// AmberParmFile::GetFlagDouble()
double *AmberParmFile::GetFlagDouble(AmberParmFlagType fflag, int maxval)
{
  int ncols, width, precision;
  FortranType fType;
  
  fType = SeekToFlag(fflag, ncols, width, precision);
  if (fType == UNKNOWN_FTYPE) return NULL;
  // NOTE: Check that type matches?
  // Read double
  return GetDouble(width, ncols, maxval);
}

// AmberParmFile::GetFlagInteger()
int *AmberParmFile::GetFlagInteger(AmberParmFlagType fflag, int maxval)
{
  int ncols, width, precision;
  FortranType fType;
  
  fType = SeekToFlag(fflag, ncols, width, precision);
  if (fType == UNKNOWN_FTYPE) return NULL;
  // NOTE: Check that type matches?
  // Read integer 
  return GetInteger(width, ncols, maxval);
}

// AmberParmFile::GetFlagName()
NAME *AmberParmFile::GetFlagName(AmberParmFlagType fflag, int maxval)
{
  int ncols, width, precision;
  FortranType fType;
  
  fType = SeekToFlag(fflag, ncols, width, precision);
  if (fType == UNKNOWN_FTYPE) return NULL;
  // NOTE: Check that type matches?
  // Read NAME 
  return GetName(width, ncols, maxval);
}

// -----------------------------------------------------------------------------

// AmberParmFile::SetParmFromValues()
/** Used by ReadParmAmber and ReadParmOldAmber to set AmberParm variables
  * from the POINTERS section of the parmtop.
  */
void AmberParmFile::SetParmFromValues(AmberParm &parmOut, int *values, bool isOld) {
  //mprintf("DEBUG: VALUES\n");
  //for (int i = 0; i < AMBERPOINTERS; i++) mprintf("\t%i\n",values[i]);
  // Set some commonly used values
  parmOut.natom = values[NATOM];
  parmOut.nres = values[NRES];
  parmOut.NbondsWithH = values[NBONH];
  parmOut.NbondsWithoutH = values[NBONA];
  if (debug>0) {
    if (isOld)
      mprintf("\tOld Amber top");
    else
      mprintf("\tAmber top");
    mprintf(", contains %i atoms, %i residues.\n",parmOut.natom,parmOut.nres);
    mprintf("\t%i bonds to hydrogen, %i other bonds.\n",
            parmOut.NbondsWithH,parmOut.NbondsWithoutH);
  }
  // Other values
  parmOut.ntypes = values[NTYPES];
  parmOut.nnb = values[NNB];
  parmOut.numbnd = values[NUMBND];
  parmOut.numang = values[NUMANG];
  parmOut.numdih = values[NPTRA];
  parmOut.NanglesWithH = values[NTHETH];
  parmOut.NanglesWithoutH = values[NTHETA];
  parmOut.NdihedralsWithH = values[NPHIH];
  parmOut.NdihedralsWithoutH = values[NPHIA];
  parmOut.natyp = values[NATYP];
  parmOut.nphb = values[NPHB];
  // Check that NBONA == MBONA etc. If not print a warning
  if (values[MBONA] != values[NBONA])
    mprintf("\tWarning: [%s] Amber parm has constraint bonds, but they will be ignored.\n",
            parmOut.parmName);
  if (values[MTHETA] != values[NTHETA])
    mprintf("\tWarning: [%s] Amber parm has constraint angles, but they will be ignored.\n",
            parmOut.parmName);
  if (values[MPHIA] != values[NPHIA])
    mprintf("\tWarning: [%s] Amber parm has constraint dihedrals, but they will be ignored.\n",
            parmOut.parmName);
  // If parm contains IFCAP or IFPERT info, print a warning since cpptraj
  // currently does not read these in.
  if (values[IFCAP] > 0)
    mprintf("\tWarning: Parm [%s] contains CAP information, which Cpptraj ignores.\n",
            parmOut.parmName);
  if (values[IFPERT] > 0)
    mprintf("\tWarning: Parm [%s] contains PERT information, which Cpptraj ignores.\n",
            parmOut.parmName);
}

// AmberParmFile::ReadParm()
/** Read parameters from Amber Topology file. */
int AmberParmFile::ReadParm(AmberParm &parmOut, CpptrajFile &parmfile) {
  int err;
# ifdef USE_CHARBUFFER
  // TEST: Close and reopen buffered.
  parmfile.CloseFile();
  parmfile.OpenFileBuffered();
# endif
  File = &parmfile; // For new STL functions.
  if (parmfile.fileFormat == OLDAMBERPARM)
    err = ReadParmOldAmber(parmOut,parmfile);
  else 
    err = ReadParmAmber(parmOut,parmfile);
  if (err!=0) return 1;
  // Common setup
  // Convert charge to units of electron charge
  if (parmOut.charge!=NULL) {
    double convert = AMBERTOELEC;
    for (int atom = 0; atom < parmOut.natom; atom++)
      parmOut.charge[atom] *= convert;
  }
  // Shift atom #s in resnums by -1 so they start from 0
  if (parmOut.resnums!=NULL) {
    for (int res = 0; res < parmOut.nres; res++)
      --parmOut.resnums[res];
  }
  // Shift atom #s in excludedAtoms by -1 so they start from 0
  if (parmOut.excludedAtoms!=NULL) {
    for (int atom = 0; atom < parmOut.nnb; atom++)
      --parmOut.excludedAtoms[atom];
  }
  // Print Box info
  if (debug>0 && parmOut.firstSolvMol!=-1) {
    mprintf("\tAmber parm %s contains box info: %i mols, first solvent mol is %i\n",
            parmOut.parmName, parmOut.molecules, parmOut.firstSolvMol);
    mprintf("\tBOX: %lf %lf %lf | %lf %lf %lf\n",
            parmOut.Box[0],parmOut.Box[1],parmOut.Box[2],
            parmOut.Box[3],parmOut.Box[4],parmOut.Box[5]);
    if (parmOut.boxType==ORTHO)
     mprintf("\t     Box is orthogonal.\n");
    else if (parmOut.boxType==NONORTHO)
      mprintf("\t     Box is non-orthogonal.\n");
    else
      mprintf("\t     Box will be determined from first associated trajectory.\n");
  }
  return 0;
}

// AmberParmFile::ReadParmOldAmber()
int AmberParmFile::ReadParmOldAmber(AmberParm &parmOut, CpptrajFile &parmfile) {
  int values[30];

  if (debug>0) mprintf("Reading Old-style Amber Topology file %s\n",parmOut.parmName);
  char *title = GetLine();
  if (debug>0) mprintf("\tOld AmberParm Title: %s\n",title);
  delete[] title;
  // Pointers - same as new format except only 30 values, no NEXTRA
  int *tempvalues = GetInteger(6,12,30);
  if (tempvalues==NULL) {
    mprintf("Could not get pointers from old amber topology file %s.\n",parmOut.parmName);
    return 1;
  }
  memcpy(values, tempvalues, 30 * sizeof(int));
  delete[] tempvalues;
  // Set some commonly used values
  SetParmFromValues(parmOut, values, true);
  // Load the rest of the parm
  parmOut.names = GetName(4,20,values[NATOM]);
  parmOut.charge = GetDouble(16,5,values[NATOM]);
  parmOut.mass = GetDouble(16,5,values[NATOM]);
  parmOut.atype_index = GetInteger(6,12,values[NATOM]);
  parmOut.numex = GetInteger(6,12,values[NATOM]);
  parmOut.NB_index = GetInteger(6,12,values[NTYPES]*values[NTYPES]);
  parmOut.resnames = GetName(4,20,values[NRES]);
  parmOut.resnums = GetInteger(6,12,values[NRES]); 
  parmOut.bond_rk = GetDouble(16,5,values[NUMBND]);
  parmOut.bond_req = GetDouble(16,5,values[NUMBND]);
  parmOut.angle_tk = GetDouble(16,5,values[NUMANG]);
  parmOut.angle_teq = GetDouble(16,5,values[NUMANG]);
  parmOut.dihedral_pk = GetDouble(16,5,values[NPTRA]);
  parmOut.dihedral_pn = GetDouble(16,5,values[NPTRA]);
  parmOut.dihedral_phase = GetDouble(16,5,values[NPTRA]);
  parmOut.solty = GetDouble(16,5,values[NATYP]);
  int nlj = values[NTYPES] * (values[NTYPES] + 1) / 2;
  parmOut.LJ_A = GetDouble(16,5,nlj);
  parmOut.LJ_B = GetDouble(16,5,nlj);
  parmOut.bondsh = GetInteger(6,12,values[NBONH]*3);
  parmOut.bonds = GetInteger(6,12,values[NBONA]*3);
  parmOut.anglesh = GetInteger(6,12,values[NTHETH]*4);
  parmOut.angles = GetInteger(6,12,values[NTHETA]*4);
  parmOut.dihedralsh = GetInteger(6,12,values[NPHIH]*5);
  parmOut.dihedrals = GetInteger(6,12,values[NPHIA]*5);
  parmOut.excludedAtoms = GetInteger(6,12,values[NNB]);
  parmOut.asol = GetDouble(16,5,values[NPHB]);
  parmOut.bsol = GetDouble(16,5,values[NPHB]);
  parmOut.hbcut = GetDouble(16,5,values[NPHB]);
  parmOut.types = GetName(4,20,values[NATOM]);
  parmOut.itree = GetName(4,20,values[NATOM]);
  parmOut.join_array = GetInteger(6,12,values[NATOM]);
  parmOut.irotat = GetInteger(6,12,values[NATOM]);
  // Solvent/Box info
  if (values[IFBOX] > 0) {
    int *solvent_pointer = GetInteger(6,12,3);
    if (solvent_pointer==NULL) {
      mprintf("Error in solvent pointers.\n");
      return 1;
    } 
    parmOut.finalSoluteRes = solvent_pointer[0];
    parmOut.molecules = solvent_pointer[1];
    parmOut.firstSolvMol = solvent_pointer[2];
    delete[] solvent_pointer;
    parmOut.atomsPerMol = GetInteger(6,12,parmOut.molecules);
    // boxFromParm = {OLDBETA, BOX(1), BOX(2), BOX(3)}
    double *boxFromParm = GetDouble(16,5,4);
    if (boxFromParm==NULL) {mprintf("Error in box info.\n"); return 1;}
    parmOut.boxType = SetBoxInfo(boxFromParm,parmOut.Box,debug);
    delete[] boxFromParm;
  }
  return 0;
}

// AmberParmFile::ReadParmAmber()
int AmberParmFile::ReadParmAmber(AmberParm &parmOut, CpptrajFile &parmfile) {
  int values[AMBERPOINTERS];
  bool chamber; // Set to true if this top is a chamber-created topology file

  if (debug>0) mprintf("Reading Amber Topology file %s\n",parmfile.filename);
  // Title
  char *title = GetFlagLine("TITLE");
  // If title is NULL, check for CTITLE (chamber parm)
  if (title==NULL) {
    title = GetFlagLine("CTITLE");
    chamber = true;
  } else {
    chamber = false;
  }
  if (debug>0) mprintf("\tAmberParm Title: %s\n",title);
  delete[] title;
  // Pointers
  int *tempvalues = GetFlagInteger(F_POINTERS,AMBERPOINTERS);
  if (tempvalues==NULL) {
    mprintf("Could not get %s from Amber Topology file.\n",AmberParmFlag[F_POINTERS]);
    return 1;
  }
  memcpy(values, tempvalues, AMBERPOINTERS * sizeof(int));
  delete[] tempvalues;
  // Set some commonly used values
  SetParmFromValues(parmOut,values, false);
  // Get parm variables
  parmOut.names = GetFlagName(F_NAMES, values[NATOM]);
  parmOut.charge = GetFlagDouble(F_CHARGE,values[NATOM]);
  parmOut.at_num = GetFlagInteger(F_ATOMICNUM,values[NATOM]);
  parmOut.mass = GetFlagDouble(F_MASS,values[NATOM]);
  parmOut.atype_index = GetFlagInteger(F_ATYPEIDX,values[NATOM]);
  parmOut.numex = GetFlagInteger(F_NUMEX,values[NATOM]);
  parmOut.NB_index = GetFlagInteger(F_NB_INDEX,values[NTYPES]*values[NTYPES]);
  parmOut.resnames = GetFlagName(F_RESNAMES,values[NRES]);
  parmOut.resnums = GetFlagInteger(F_RESNUMS,values[NRES]);
  parmOut.bond_rk = GetFlagDouble(F_BONDRK, values[NUMBND]);
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
  parmOut.LJ_B = GetFlagDouble(F_LJ_B,nlj);
  parmOut.bondsh = GetFlagInteger(F_BONDSH,values[NBONH]*3);
  parmOut.bonds = GetFlagInteger(F_BONDS,values[NBONA]*3);
  parmOut.anglesh = GetFlagInteger(F_ANGLESH, values[NTHETH]*4);
  parmOut.angles = GetFlagInteger(F_ANGLES, values[NTHETA]*4);
  parmOut.dihedralsh = GetFlagInteger(F_DIHH, values[NPHIH]*5);
  parmOut.dihedrals = GetFlagInteger(F_DIH, values[NPHIA]*5);
  parmOut.excludedAtoms = GetFlagInteger(F_EXCLUDE, values[NNB]);
  parmOut.asol = GetFlagDouble(F_ASOL,values[NPHB]);
  parmOut.bsol = GetFlagDouble(F_BSOL,values[NPHB]);
  parmOut.hbcut = GetFlagDouble(F_HBCUT,values[NPHB]);
  parmOut.types = GetFlagName(F_TYPES,values[NATOM]);
  parmOut.itree = GetFlagName(F_ITREE,values[NATOM]);
  parmOut.join_array = GetFlagInteger(F_JOIN,values[NATOM]);
  parmOut.irotat = GetFlagInteger(F_IROTAT,values[NATOM]);
  // Get solvent info if IFBOX>0
  if (values[IFBOX]>0) {
    int *solvent_pointer = GetFlagInteger(F_SOLVENT_POINTER,3);
    if (solvent_pointer==NULL) {
      mprintf("Could not get %s from Amber Topology file.\n",AmberParmFlag[F_SOLVENT_POINTER]);
      return 1;
    } 
    parmOut.finalSoluteRes = solvent_pointer[0];
    parmOut.molecules = solvent_pointer[1];
    parmOut.firstSolvMol = solvent_pointer[2];
    delete[] solvent_pointer;
    parmOut.atomsPerMol = GetFlagInteger(F_ATOMSPERMOL,parmOut.molecules);
    // boxFromParm = {OLDBETA, BOX(1), BOX(2), BOX(3)}
    double *boxFromParm = GetFlagDouble(F_PARMBOX,4);
    // If no box information present in the parm (such as with Chamber prmtops)
    // set the box info if ifbox = 2, otherwise set to NOBOX; the box info will 
    // eventually be set by angles from the first trajectory associated with 
    // this parm.
    if (boxFromParm==NULL) {
      if (not chamber) mprintf("Warning: Prmtop missing Box information.\n");
      // ifbox 2: truncated octahedron for certain
      if (values[IFBOX] == 2) {
        parmOut.boxType = NONORTHO;
        parmOut.Box[0] = 0.0; 
        parmOut.Box[1] = 0.0; 
        parmOut.Box[2] = 0.0;
        parmOut.Box[3] = TRUNCOCTBETA;
        parmOut.Box[4] = TRUNCOCTBETA;
        parmOut.Box[5] = TRUNCOCTBETA;
      } else
        parmOut.boxType = NOBOX;
    // Determine box type, set Box angles and lengths from beta (boxFromParm[0])
    } else {
      parmOut.boxType = SetBoxInfo(boxFromParm,parmOut.Box,debug);
      delete[] boxFromParm;
    }
  }
  // GB parameters; radius set, radii, and screening parameters
  parmOut.radius_set = GetFlagLine("RADIUS_SET");
  if (debug>0) mprintf("\tRadius Set: %s\n",parmOut.radius_set);
  parmOut.gb_radii = GetFlagDouble(F_RADII,values[NATOM]);
  parmOut.gb_screen = GetFlagDouble(F_SCREEN,values[NATOM]);
  
  return error_count;
}

// AmberParmFile::WriteParm()
/** Write out information from current AmberParm to an Amber parm file */
int AmberParmFile::WriteParm(AmberParm &parmIn, CpptrajFile &outfile) {
  CharBuffer buffer;
  int solvent_pointer[3];
  int values[AMBERPOINTERS];
  double parmBox[4];
  // For date and time
  time_t rawtime;
  struct tm *timeinfo;
  int largestRes=0; // For determining nmxrs
  int *tempResnums = NULL;
  double *tempCharge = NULL;

  if (parmIn.parmName==NULL) return 1;

  if (!outfile.IsOpen()) return 1;

  // HEADER AND TITLE (4 lines, version, flag, format, title)
  buffer.Allocate( 324 ); // (81 * 4), no space for NULL needed since using NewLine() 
  time(&rawtime);
  timeinfo = localtime(&rawtime);
  // VERSION
  // NOTE: tm_mon ranges from 0
  buffer.Sprintf("%-44s%02i/%02i/%02i  %02i:%02i:%02i                  \n",
                     "%VERSION  VERSION_STAMP = V0001.000  DATE = ",
                     timeinfo->tm_mon+1,timeinfo->tm_mday,timeinfo->tm_year%100,
                     timeinfo->tm_hour,timeinfo->tm_min,timeinfo->tm_sec);
  // TITLE
  buffer.Sprintf("%-80s\n%-80s\n%-80s","%FLAG TITLE","%FORMAT(20a4)","");
  buffer.NewLine();
  //outfile.IO->Printf("%-80s\n",parmName);

  // Shift atom #s in resnums by +1 to be consistent with AMBER
  // Also determine # atoms in largest residue for nmxrs
  tempResnums = new int[ parmIn.nres ];
  memcpy(tempResnums, parmIn.resnums, parmIn.nres * sizeof(int));
  for (int res=0; res < parmIn.nres; res++) {
    int diff = parmIn.resnums[res+1] - parmIn.resnums[res];
    if (diff > largestRes) largestRes = diff;
    tempResnums[res] += 1;
  }
  // POINTERS
  memset(values, 0, AMBERPOINTERS * sizeof(int));
  values[NATOM]=parmIn.natom;
  values[NTYPES]=parmIn.ntypes;
  values[NBONH]=parmIn.NbondsWithH;
  values[MBONA]=parmIn.NbondsWithoutH;
  values[NTHETH]=parmIn.NanglesWithH;
  values[MTHETA]=parmIn.NanglesWithoutH;
  values[NPHIH]=parmIn.NdihedralsWithH;
  values[MPHIA]=parmIn.NdihedralsWithoutH;
  values[NNB]=parmIn.nnb;
  values[NRES]=parmIn.nres;
  //   NOTE: Assuming NBONA == MBONA etc
  values[NBONA]=parmIn.NbondsWithoutH;
  values[NTHETA]=parmIn.NanglesWithoutH;
  values[NPHIA]=parmIn.NdihedralsWithoutH;
  values[NUMBND]=parmIn.numbnd;
  values[NUMANG]=parmIn.numang;
  values[NPTRA]=parmIn.numdih;
  values[NATYP]=parmIn.natyp;
  values[NPHB]=parmIn.nphb;
  values[IFBOX]=AmberIfbox(parmIn.Box[4]);
  values[NMXRS]=largestRes;
  DataToFortranBuffer(buffer,F_POINTERS, values, NULL, NULL, AMBERPOINTERS);
  // ATOM NAMES
  DataToFortranBuffer(buffer,F_NAMES, NULL, NULL, parmIn.names, parmIn.natom);
  // CHARGE - might be null if read from pdb
  if (parmIn.charge!=NULL) {
    // Convert charges to AMBER charge units
    tempCharge = new double[ parmIn.natom ];
    memcpy(tempCharge, parmIn.charge, parmIn.natom * sizeof(double));
    for (int atom=0; atom<parmIn.natom; atom++)
      tempCharge[atom] *= (ELECTOAMBER);
    DataToFortranBuffer(buffer,F_CHARGE, NULL, tempCharge, NULL, parmIn.natom);
    delete[] tempCharge;
  }
  // ATOMIC NUMBER
  if (parmIn.at_num!=NULL)
    DataToFortranBuffer(buffer,F_ATOMICNUM, parmIn.at_num, NULL, NULL, parmIn.natom);
  // MASS - might be null if read from pdb
  if (parmIn.mass!=NULL)
    DataToFortranBuffer(buffer,F_MASS, NULL, parmIn.mass, NULL, parmIn.natom);
  // ATOM_TYPE_INDEX
  if (parmIn.atype_index!=NULL)
    DataToFortranBuffer(buffer,F_ATYPEIDX, parmIn.atype_index, NULL, NULL, parmIn.natom);
  // NUMBER_EXCLUDED_ATOMS
  if (parmIn.numex!=NULL)
    DataToFortranBuffer(buffer,F_NUMEX, parmIn.numex, NULL, NULL, parmIn.natom);
  // NONBONDED_PARM_INDEX
  if (parmIn.NB_index!=NULL)
    DataToFortranBuffer(buffer,F_NB_INDEX, parmIn.NB_index, NULL, NULL, parmIn.ntypes * parmIn.ntypes);
  // RESIDUE LABEL - resnames
  DataToFortranBuffer(buffer,F_RESNAMES, NULL, NULL, parmIn.resnames, parmIn.nres);
  // RESIDUE POINTER - resnums, IPRES; tempResnums is shifted +1, see above
  DataToFortranBuffer(buffer,F_RESNUMS, tempResnums, NULL, NULL, parmIn.nres);
  delete[] tempResnums;
  // BOND_FORCE_CONSTANT and EQUIL VALUES
  if (parmIn.bond_rk!=NULL)
    DataToFortranBuffer(buffer,F_BONDRK, NULL, parmIn.bond_rk, NULL, parmIn.numbnd);
  if (parmIn.bond_req!=NULL)
    DataToFortranBuffer(buffer,F_BONDREQ, NULL, parmIn.bond_req, NULL, parmIn.numbnd);
  // ANGLE FORCE CONSTANT AND EQUIL VALUES
  if (parmIn.angle_tk!=NULL)
    DataToFortranBuffer(buffer,F_ANGLETK, NULL, parmIn.angle_tk, NULL, parmIn.numang);
  if (parmIn.angle_teq!=NULL)
    DataToFortranBuffer(buffer,F_ANGLETEQ, NULL, parmIn.angle_teq, NULL, parmIn.numang);
  // DIHEDRAL CONSTANT, PERIODICITY, PHASE
  if (parmIn.dihedral_pk!=NULL)
    DataToFortranBuffer(buffer,F_DIHPK, NULL, parmIn.dihedral_pk, NULL, parmIn.numdih);
  if (parmIn.dihedral_pn!=NULL)
    DataToFortranBuffer(buffer,F_DIHPN, NULL, parmIn.dihedral_pn, NULL, parmIn.numdih);
  if (parmIn.dihedral_phase!=NULL)
    DataToFortranBuffer(buffer,F_DIHPHASE, NULL, parmIn.dihedral_phase, NULL, parmIn.numdih);
  // SCEE and SCNB scaling factors
  if (parmIn.scee_scale!=NULL)
    DataToFortranBuffer(buffer,F_SCEE, NULL, parmIn.scee_scale, NULL, parmIn.numdih);
  if (parmIn.scnb_scale!=NULL)
    DataToFortranBuffer(buffer,F_SCNB, NULL, parmIn.scnb_scale, NULL, parmIn.numdih);
  // SOLTY
  if (parmIn.solty!=NULL)
    DataToFortranBuffer(buffer,F_SOLTY, NULL, parmIn.solty, NULL, parmIn.natyp);
  // Lennard-Jones A/B
  if (parmIn.LJ_A!=NULL)
    DataToFortranBuffer(buffer,F_LJ_A, NULL, parmIn.LJ_A, NULL, parmIn.ntypes*(parmIn.ntypes+1)/2);
  if (parmIn.LJ_B!=NULL)
    DataToFortranBuffer(buffer,F_LJ_B, NULL, parmIn.LJ_B, NULL, parmIn.ntypes*(parmIn.ntypes+1)/2);
  // BONDS INCLUDING HYDROGEN - might be null if read from pdb
  if (parmIn.bondsh != NULL)
    DataToFortranBuffer(buffer,F_BONDSH, parmIn.bondsh, NULL, NULL, parmIn.NbondsWithH*3);
  // BONDS WITHOUT HYDROGEN - might be null if read from pdb
  if (parmIn.bonds!=NULL)
    DataToFortranBuffer(buffer,F_BONDS, parmIn.bonds, NULL, NULL, parmIn.NbondsWithoutH*3);
  // ANGLES INCLUDING HYDROGEN
  if (parmIn.anglesh!=NULL)
    DataToFortranBuffer(buffer,F_ANGLESH, parmIn.anglesh, NULL, NULL, parmIn.NanglesWithH*4);
  // ANGLES WITHOUT HYDROGEN
  if (parmIn.angles!=NULL)
    DataToFortranBuffer(buffer,F_ANGLES, parmIn.angles, NULL, NULL, parmIn.NanglesWithoutH*4);
  // DIHEDRALS INCLUDING HYDROGEN
  if (parmIn.dihedralsh!=NULL)
    DataToFortranBuffer(buffer,F_DIHH, parmIn.dihedralsh, NULL, NULL, parmIn.NdihedralsWithH*5);
  // DIHEDRALS WITHOUT H
  if (parmIn.dihedrals!=NULL)
    DataToFortranBuffer(buffer,F_DIH, parmIn.dihedrals, NULL, NULL, parmIn.NdihedralsWithoutH*5);
  // EXCLUDED ATOMS LIST
  // Shift atom #s in excludedAtoms by +1 to be consistent with AMBER
  if (parmIn.excludedAtoms!=NULL) {
    int *tempexclude = new int[ parmIn.nnb ];
    memcpy(tempexclude, parmIn.excludedAtoms, parmIn.nnb * sizeof(int));
    for (int atom = 0; atom < parmIn.nnb; atom++) tempexclude[atom] += 1;
    DataToFortranBuffer(buffer,F_EXCLUDE, tempexclude, NULL, NULL, parmIn.nnb);
    delete[] tempexclude;
  }
  // HBOND ACOEFF, BCOEFF, HBCUT
  // NAB requires that these be printed even if no values exist.
  DataToFortranBuffer(buffer,F_ASOL,NULL,parmIn.asol,NULL,parmIn.nphb);
  DataToFortranBuffer(buffer,F_BSOL,NULL,parmIn.bsol,NULL,parmIn.nphb);
  DataToFortranBuffer(buffer,F_HBCUT,NULL,parmIn.hbcut,NULL,parmIn.nphb);
  // AMBER ATOM TYPE - might be null if read from pdb
  if (parmIn.types!=NULL)
    DataToFortranBuffer(buffer,F_TYPES, NULL, NULL, parmIn.types, parmIn.natom);
  // TREE CHAIN CLASSIFICATION
  if (parmIn.itree!=NULL)
    DataToFortranBuffer(buffer,F_ITREE, NULL, NULL, parmIn.itree, parmIn.natom);
  // JOIN ARRAY
  if (parmIn.join_array!=NULL)
    DataToFortranBuffer(buffer,F_JOIN, parmIn.join_array, NULL, NULL, parmIn.natom);
  // IROTAT
  if (parmIn.irotat!=NULL)
    DataToFortranBuffer(buffer,F_IROTAT, parmIn.irotat, NULL, NULL, parmIn.natom);
  // Write solvent info if IFBOX>0
  if (values[IFBOX]>0) {
    // SOLVENT POINTERS
    if (parmIn.firstSolvMol!=-1) {
      solvent_pointer[0]=parmIn.finalSoluteRes;
      solvent_pointer[1]=parmIn.molecules;
      solvent_pointer[2]=parmIn.firstSolvMol;
    } else {
      solvent_pointer[0]=parmIn.nres;
      solvent_pointer[1]=parmIn.molecules;
      solvent_pointer[2]=parmIn.molecules+1;
    }
    DataToFortranBuffer(buffer,F_SOLVENT_POINTER, solvent_pointer, NULL, NULL, 3);
    // ATOMS PER MOLECULE
    if (parmIn.atomsPerMol!=NULL)
      DataToFortranBuffer(buffer,F_ATOMSPERMOL, parmIn.atomsPerMol, NULL, NULL, parmIn.molecules);
    // BOX DIMENSIONS
    parmBox[0] = parmIn.Box[4]; // beta
    parmBox[1] = parmIn.Box[0]; // boxX
    parmBox[2] = parmIn.Box[1]; // boxY
    parmBox[3] = parmIn.Box[2]; // boxZ
    DataToFortranBuffer(buffer,F_PARMBOX, NULL, parmBox, NULL, 4);
  }
  // GB RADIUS SET
  if (parmIn.radius_set!=NULL) {
    buffer.IncreaseSize(243); // 3 * 81
    buffer.Sprintf("%-80s\n%-80s\n%-80s","%FLAG RADIUS_SET","%FORMAT(1a80)",parmIn.radius_set);
    buffer.NewLine();
  }
  // GB RADII
  if (parmIn.gb_radii!=NULL)
    DataToFortranBuffer(buffer,F_RADII, NULL, parmIn.gb_radii, NULL, parmIn.natom);
  // GB SCREENING PARAMETERS
  if (parmIn.gb_screen!=NULL)
    DataToFortranBuffer(buffer,F_SCREEN, NULL, parmIn.gb_screen, NULL, parmIn.natom);

  // Write buffer to file
  outfile.IO->Write(buffer.Buffer(), sizeof(char), buffer.CurrentSize());
  outfile.CloseFile();

  return 0;
}

