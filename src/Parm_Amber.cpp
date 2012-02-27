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

// ---------- Defines and Enumerated types --------------------------------------
/* Compiler Defines:
 * - USE_CHARBUFFER: Use CharBuffer to buffer entire file
 */
/// Enumerated type for Amber Parmtop Flags
enum AmberParmFlagType {
  F_POINTERS = 0, F_NAMES,   F_CHARGE,  F_MASS,    F_RESNAMES,
  F_RESNUMS,      F_TYPES,   F_BONDSH,  F_BONDS,   F_SOLVENT_POINTER,
  F_ATOMSPERMOL,  F_PARMBOX, F_ATYPEIDX,F_NUMEX,   F_NB_INDEX,
  F_LJ_A,         F_LJ_B,    F_EXCLUDE, F_RADII,   F_SCREEN,
  F_BONDRK,       F_BONDREQ, F_ANGLETK, F_ANGLETEQ,F_DIHPK,
  F_DIHPN,        F_DIHPHASE,F_SCEE,    F_SCNB,    F_SOLTY,
  F_ANGLESH,      F_ANGLES,  F_DIHH,    F_DIH,     F_ASOL,
  F_BSOL,         F_HBCUT,   F_ITREE,   F_JOIN,    F_IROTAT
};
#define NUMAMBERPARMFLAGS 40
/// Enumerated type for FLAG_POINTERS section
enum topValues {
//0       1       2      3       4       5       6       7      8       9
  NATOM,  NTYPES, NBONH, MBONA,  NTHETH, MTHETA, NPHIH,  MPHIA, NHPARM, NPARM,
  NNB,    NRES,   NBONA, NTHETA, NPHIA,  NUMBND, NUMANG, NPTRA, NATYP,  NPHB,
  IFPERT, NBPER,  NGPER, NDPER,  MBPER,  MGPER,  MDPER,  IFBOX, NMXRS, IFCAP,
  NEXTRA
};
#define AMBERPOINTERS 31
/// Enumerated type for Fortran data type
enum FortranType {
  UNKNOWN_FTYPE, FINT, FDOUBLE, FCHAR, FFLOAT
};
/// Combined size of %FLAG and %FORMAT lines (81 * 2)
#define FFSIZE 162
/// Constant strings for fortran formats corresponding to Amber parm flags
static const char AmberParmFmt[NUMAMBERPARMFLAGS][16] = {
"%FORMAT(10I8)",   "%FORMAT(20a4)",   "%FORMAT(5E16.8)", "%FORMAT(5E16.8)", "%FORMAT(20a4)",
"%FORMAT(10I8)",   "%FORMAT(20a4)",   "%FORMAT(10I8)",   "%FORMAT(10I8)",   "%FORMAT(3I8)",
"%FORMAT(10I8)",   "%FORMAT(5E16.8)", "%FORMAT(10I8)",   "%FORMAT(10I8)",   "%FORMAT(10I8)",
"%FORMAT(5E16.8)", "%FORMAT(5E16.8)", "%FORMAT(10I8)",   "%FORMAT(5E16.8)", "%FORMAT(5E16.8)",
"%FORMAT(5E16.8)", "%FORMAT(5E16.8)", "%FORMAT(5E16.8)", "%FORMAT(5E16.8)", "%FORMAT(5E16.8)",
"%FORMAT(5E16.8)", "%FORMAT(5E16.8)", "%FORMAT(5E16.8)", "%FORMAT(5E16.8)", "%FORMAT(5E16.8)",
"%FORMAT(10I8)",   "%FORMAT(10I8)",   "%FORMAT(10I8)",   "%FORMAT(10I8)",   "%FORMAT(5E16.8)",
"%FORMAT(5E16.8)", "%FORMAT(5E16.8)", "%FORMAT(20a4)",   "%FORMAT(10I8)",   "%FORMAT(10I8)"
};
/// Constant strings for Amber parm flags
static const char AmberParmFlag[NUMAMBERPARMFLAGS][27] = {
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
  "IROTAT"
};
// -----------------------------------------------------------------------------

// ---------- INTERNAL FUNCTIONS -----------------------------------------------
// GetFortranBufferSize()
/** Given number of columns and the width of each column, return the 
  * necessary char buffer size for N data elements.
  */
// NOTE: Convert to size_t?
static int GetFortranBufferSize(int N, int isDos, int width, int ncols) {
  int bufferLines=0;
  int BufferSize=0;

  BufferSize= N * width;
  bufferLines = N / ncols;
  if ((N % ncols)!=0) bufferLines++;
  // If DOS file there are CR before Newlines
  if (isDos) bufferLines*=2;
  BufferSize+=bufferLines;
  //if (debug>0) 
  //  fprintf(stdout,"*** Buffer size is %i including %i newlines.\n",BufferSize,bufferLines);
  return BufferSize;
}

// GetFortranType()
/** Given a fortran-type format string, return the corresponding fortran
  * type. Set ncols (if present), width, and precision (if present).
  */
// 01234567
// %FORMAT([<cols>][(]<type><width>[<precision>][)])
static FortranType GetFortranType(char *FormatIn, int *ncols, int *width, int *precision) {
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
  ptr++;
  *ptr='\0';
}

// getFlagFileString()
/** Search for the FLAG specified by Key. Assume the next line is a string 
  * of max length 80 chars and return it.
  */
static char *getFlagFileString(CpptrajFile &File, const char *Key, int debug) {
  char *lineBuffer;

  // Find flag
  if (PositionFileAtFlag(File,Key,NULL,debug)) {
    // Read next line and return
    lineBuffer = new char[83]; // 80 + newline + NULL ( + CR if dos)
#   ifdef USE_CHARBUFFER
    File.Gets(lineBuffer,82);
#   else
    File.IO->Gets(lineBuffer,82);
#   endif
    RemoveWhitespace(lineBuffer);
    return lineBuffer;
  }

  return NULL;
}

// getFlagFileValues()
/** Search for the FLAG specified by Key and return the values. The values will 
  * be put into an array type according to the FORMAT string in the top file, 
  * but it is necessary to explictly type the returned array. maxval is used to 
  * allocate memory for the return array - only maxval values will be read.
  */
static void *getFlagFileValues(CpptrajFile &File, AmberParmFlagType fflag, int maxval, int debug){
  int ncols, width, precision;
  char fformat[83];    // Hold Format from FORMAT line
  char temp[17];       // Hold data element
  char *buffer,*ptr,*Key;
  NAME *C = NULL;
  int *I = NULL;
  double *D = NULL;
  FortranType fType;
  int BufferSize;

  // Get Flag Key
  Key = (char*) AmberParmFlag[fflag];

  // Find flag, get format
  if (!PositionFileAtFlag(File,Key,fformat,debug)) return NULL;
  if (debug>1) mprintf("DEBUG: maxval= %i\n",maxval);

  // Determine cols, width etc from format
  fType = GetFortranType(fformat, &ncols, &width, &precision);
  if (debug>1) mprintf("DEBUG: Format type %i\n",(int)fType);
  if (fType == UNKNOWN_FTYPE) {
    mprinterr("Error: Unrecognized fortran format [%s] for key [%s]\n",fformat,Key);
    return NULL;
  }
  // Allocate memory based on data type
  switch (fType) {
    case UNKNOWN_FTYPE : return NULL;
    case FINT   : I = new int[ maxval ];    break;
    case FDOUBLE: D = new double[ maxval ]; break;
    case FCHAR  : C = new NAME[ maxval ];   break;
    case FFLOAT : D = new double[ maxval ]; break;
  }
  // Allocate memory to read in entire section
  BufferSize = GetFortranBufferSize(maxval, File.isDos, width, ncols);
  buffer = new char[ BufferSize ];
# ifdef USE_CHARBUFFER
  if ( File.Read(buffer,BufferSize)==-1 )
# else
  if ( File.IO->Read(buffer,sizeof(char),BufferSize)==-1 )
# endif
  {
    rprinterr("ERROR in read of prmtop section %s\n",Key);
    delete[] buffer;
    if (I!=NULL) delete[] I;
    if (D!=NULL) delete[] D;
    if (C!=NULL) delete[] C;
    return NULL;
  }
  if (debug>3) mprintf("DEBUG: Buffer [%s]\n",buffer);

  // Convert values in buffer to their type
  ptr=buffer;
  temp[width]='\0';
  for (int i=0; i<maxval; i++) {
    // Advance past newlines - DOS newline is CR
    //fprintf(stdout,"0 %i: %c %i\n",i,ptr[0],ptr[0]);
    while (ptr[0]=='\n' || ptr[0]=='\r') { ptr++; }
    //fprintf(stdout,"1 %i: %c %i\n",i,ptr[0],ptr[0]);
    strncpy(temp,ptr,width);
    if (debug>3) mprintf("DEBUG:   %8i buffer %s\n",i,temp);
    // Sanity check: If we have hit another FLAG before we've read maxval
    // values this is bad.
    // NOTE: Could do addtional type checking too.
    if ( strncmp(ptr,"%FLAG",5)==0 ) {
      rprintf("Error: #values read (%i) < # expected values (%i).\n",i,maxval);
      if (I!=NULL) delete[] I;
      if (D!=NULL) delete[] D;
      if (C!=NULL) delete[] C;
      delete[] buffer;
      return NULL;
    }
    // Now convert temp to appropriate type
    if (I!=NULL) I[i]=atoi(temp);
    if (D!=NULL) D[i]=atof(temp);
    if (C!=NULL) strcpy(C[i],temp);
    ptr+=width;
  } // END loop over maxval
  // All done! 
  delete[] buffer;
  if (I!=NULL) return I;
  if (D!=NULL) return D;
  if (C!=NULL) return C;

  return NULL;
}

// F_load20a4()
static char *F_load20a4(CpptrajFile &File) {
  char *lineBuffer = new char[83]; // 80 + newline + NULL ( + CR if dos)
# ifdef USE_CHARBUFFER
  File.Gets(lineBuffer,82);
# else
  File.IO->Gets(lineBuffer,82);
# endif
  RemoveWhitespace(lineBuffer);
  return lineBuffer;
}

// F_loadFormat()
static void *F_loadFormat(CpptrajFile &File, FortranType fType, int width, int ncols,
                   int maxval, int debug) {
  int *I = NULL;
  double *D = NULL;
  NAME *C = NULL;
  int BufferSize;
  char *ptr,*buffer;
  char temp[17];       // Hold data element
  // If # expected values is 0 there will still be a newline placeholder
  // in the parmtop. Read past that and return NULL 
  if (maxval==0) {
#   ifdef USE_CHARBUFFER
    File.Gets(temp,16);
#   else
    File.IO->Gets(temp,16);
#   endif
    return NULL;
  }
  // Allocate memory based on data type
  switch (fType) {
    case UNKNOWN_FTYPE : return NULL;
    case FINT   : I = new int[ maxval ] ;   break;
    case FDOUBLE: D = new double[ maxval ]; break;
    case FCHAR  : C = new NAME[ maxval ];   break;
    case FFLOAT : D = new double[ maxval ]; break;
  }
  // Allocate memory to read in entire section
  BufferSize = GetFortranBufferSize(maxval, File.isDos, width, ncols);
  buffer = new char[ BufferSize ];
# ifdef USE_CHARBUFFER
  if ( File.Read(buffer,BufferSize)==-1 )
# else
  if ( File.IO->Read(buffer,sizeof(char),BufferSize)==-1 )
# endif
  {
    rprinterr("ERROR in read of prmtop section; width=%i ncols=%i maxval=%i\n",
              width,ncols,maxval);
    delete[] buffer;
    if (I!=NULL) delete[] I;
    if (D!=NULL) delete[] D;
    if (C!=NULL) delete[] C;
    return NULL;
  }
  //mprintf("DEBUG: fType=%i width=%i ncols=%i maxval=%i\n",(int)fType,
  //        width,ncols,maxval);
  if (debug>3) mprintf("DEBUG: Buffer [%s]\n",buffer);

  // Convert values in buffer to their type
  ptr=buffer;
  temp[width]='\0';
  for (int i=0; i<maxval; i++) {
    // Advance past newlines - DOS newline is CR
    //fprintf(stdout,"0 %i: %c %i\n",i,ptr[0],ptr[0]);
    while (ptr[0]=='\n' || ptr[0]=='\r') { ptr++; }
    //fprintf(stdout,"1 %i: %c %i\n",i,ptr[0],ptr[0]);
    strncpy(temp,ptr,width);
    if (debug>3) mprintf("DEBUG:   %8i buffer %s\n",i,temp);
    // Now convert temp to appropriate type
    if (I!=NULL) I[i]=atoi(temp);
    if (D!=NULL) D[i]=atof(temp);
    if (C!=NULL) strcpy(C[i],temp);
    ptr+=width;
  } // END loop over maxval
  // All done! 
  delete[] buffer;
  if (I!=NULL) return I;
  if (D!=NULL) return D;
  if (C!=NULL) return C;

  return NULL;
}

// DataToFortranBuffer()
/** Write N data elements stored in I, D, or C to character buffer with given 
  * fortran format.
  */
static int DataToFortranBuffer(CharBuffer &buffer, AmberParmFlagType fFlag,
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

// -----------------------------------------------------------------------------

// AmberParmFile::SetParmFromValues()
/** Used by ReadParmAmber and ReadParmOldAmber to set AmberParm variables
  * from the POINTERS section of the parmtop.
  */
void AmberParmFile::SetParmFromValues(AmberParm &parmOut, int *values, bool isOld) {
  // Set some commonly used values
  parmOut.natom=values[NATOM];
  parmOut.nres=values[NRES];
  //ifbox=values[IFBOX];
  parmOut.NbondsWithH=values[NBONH];
  parmOut.NbondsWithoutH=values[NBONA];
  if (debug>0) {
    if (isOld)
      mprintf("    Old Amber top");
    else
      mprintf("    Amber top");
    mprintf("contains %i atoms, %i residues.\n",parmOut.natom,parmOut.nres);
    mprintf("    %i bonds to hydrogen, %i other bonds.\n",parmOut.NbondsWithH,parmOut.NbondsWithoutH);
  }
  // Other values
  parmOut.ntypes = values[NTYPES];
  parmOut.nnb = values[NNB];
  parmOut.numbnd = values[NUMBND];
  parmOut.numang = values[NUMANG];
  parmOut.numdih = values[NPTRA];
  parmOut.NanglesWithH=values[NTHETH];
  parmOut.NanglesWithoutH=values[NTHETA];
  parmOut.NdihedralsWithH=values[NPHIH];
  parmOut.NdihedralsWithoutH=values[NPHIA];
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
  if (parmfile.fileFormat == OLDAMBERPARM)
    return ReadParmOldAmber(parmOut,parmfile);
  return ReadParmAmber(parmOut,parmfile);
}

// AmberParmFile::ReadParmOldAmber()
int AmberParmFile::ReadParmOldAmber(AmberParm &parmOut, CpptrajFile &parmfile) {
  char *title;
  int values[30], ifbox;
# ifdef USE_CHARBUFFER
  // TEST: Close and reopen buffered.
  parmfile.CloseFile();
  parmfile.OpenFileBuffered();
# endif
  if (debug>0) mprintf("Reading Old-style Amber Topology file %s\n",parmOut.parmName);
  title = F_load20a4(parmfile);
  if (debug>0) mprintf("\tOld AmberParm Title: %s\n",title);
  delete[] title;
  // Pointers - same as new format except only 30 values, no NEXTRA
  int *tempvalues = (int*) F_loadFormat(parmfile, FINT, 6, 12, 30, debug);
  if (tempvalues==NULL) {
    mprintf("Could not get values from topfile\n");
    return 1;
  }
  memcpy(values, tempvalues, 30 * sizeof(int));
  delete[] tempvalues;
  // Set some commonly used values
  SetParmFromValues(parmOut, values, true);
  ifbox=values[IFBOX];
  // Load the rest of the parm
  // NOTE: Add error checking!
  parmOut.names = (NAME*) F_loadFormat(parmfile, FCHAR, 4, 20, parmOut.natom, debug);
  parmOut.charge = (double*) F_loadFormat(parmfile, FDOUBLE, 16, 5, parmOut.natom, debug);
  parmOut.mass = (double*) F_loadFormat(parmfile, FDOUBLE, 16, 5, parmOut.natom, debug);
  parmOut.atype_index = (int*) F_loadFormat(parmfile,FINT, 6, 12, parmOut.natom, debug);
  parmOut.numex = (int*) F_loadFormat(parmfile,FINT, 6, 12, parmOut.natom, debug);
  parmOut.NB_index = (int*) F_loadFormat(parmfile,FINT, 6, 12, parmOut.ntypes*parmOut.ntypes, debug);
  parmOut.resnames = (NAME*) F_loadFormat(parmfile, FCHAR, 4, 20, parmOut.nres, debug);
  parmOut.resnums = (int*) F_loadFormat(parmfile,FINT, 6, 12, parmOut.nres, debug);
  // Atom #s in resnums are currently shifted +1. Shift back to be consistent
  // with the rest of cpptraj.
  for (int atom=0; atom < parmOut.nres; atom++)
    parmOut.resnums[atom] -= 1;
  // Bond, angle, dihedral constants and values
  parmOut.bond_rk = (double*) F_loadFormat(parmfile,FDOUBLE,16,5,values[NUMBND],debug);
  parmOut.bond_req = (double*) F_loadFormat(parmfile,FDOUBLE,16,5,values[NUMBND],debug);
  parmOut.angle_tk = (double*) F_loadFormat(parmfile,FDOUBLE,16,5,values[NUMANG],debug);
  parmOut.angle_teq = (double*) F_loadFormat(parmfile,FDOUBLE,16,5,values[NUMANG],debug);
  parmOut.dihedral_pk = (double*) F_loadFormat(parmfile,FDOUBLE,16,5,values[NPTRA],debug);
  parmOut.dihedral_pn = (double*) F_loadFormat(parmfile,FDOUBLE,16,5,values[NPTRA],debug);
  parmOut.dihedral_phase = (double*) F_loadFormat(parmfile,FDOUBLE,16,5,values[NPTRA],debug);
  parmOut.solty = (double*) F_loadFormat(parmfile,FDOUBLE,16,5,values[NATYP],debug);
  // LJ params
  parmOut.LJ_A = (double*) F_loadFormat(parmfile,FDOUBLE,16,5,parmOut.ntypes*(parmOut.ntypes+1)/2,debug);
  parmOut.LJ_B = (double*) F_loadFormat(parmfile,FDOUBLE,16,5,parmOut.ntypes*(parmOut.ntypes+1)/2,debug);
  // Bonds, angles, dihedrals
  parmOut.bondsh = (int*) F_loadFormat(parmfile,FINT,6,12,values[NBONH]*3,debug);
  parmOut.bonds = (int*) F_loadFormat(parmfile,FINT,6,12,values[NBONA]*3,debug);
  parmOut.anglesh = (int*) F_loadFormat(parmfile,FINT,6,12,values[NTHETH]*4,debug);
  parmOut.angles = (int*) F_loadFormat(parmfile,FINT,6,12,values[NTHETA]*4,debug);
  parmOut.dihedralsh = (int*) F_loadFormat(parmfile,FINT,6,12,values[NPHIH]*5,debug);
  parmOut.dihedrals = (int*) F_loadFormat(parmfile,FINT,6,12,values[NPHIA]*5,debug);
  // Excluded atoms; shift by -1 so atom #s start from 0
  parmOut.excludedAtoms = (int*) F_loadFormat(parmfile,FINT,6,12,parmOut.nnb,debug);
  for (int atom=0; atom < parmOut.nnb; atom++)
    parmOut.excludedAtoms[atom] -= 1;
  // LJ 10-12 stuff 
  parmOut.asol = (double*) F_loadFormat(parmfile,FDOUBLE,16,5,values[NPHB],debug);
  parmOut.bsol = (double*) F_loadFormat(parmfile,FDOUBLE,16,5,values[NPHB],debug);
  parmOut.hbcut = (double*) F_loadFormat(parmfile,FDOUBLE,16,5,values[NPHB],debug);
  // Atom types
  parmOut.types = (NAME*) F_loadFormat(parmfile,FCHAR,4,20,parmOut.natom,debug);
  // Tree, join, irotat 
  parmOut.itree = (NAME*) F_loadFormat(parmfile,FCHAR,4,20,parmOut.natom,debug);
  parmOut.join_array = (int*) F_loadFormat(parmfile,FINT,6,12,parmOut.natom,debug);
  parmOut.irotat = (int*) F_loadFormat(parmfile,FINT,6,12,parmOut.natom,debug);
  // Solvent/Box info
  if (ifbox > 0) {
    int *solvent_pointer=(int*) F_loadFormat(parmfile,FINT,6,12,3,debug);
    if (solvent_pointer==NULL) {
      mprintf("Error in solvent pointers.\n");
      return 1;
    } else {
      parmOut.finalSoluteRes=solvent_pointer[0];
      parmOut.molecules=solvent_pointer[1];
      parmOut.firstSolvMol=solvent_pointer[2];
      delete[] solvent_pointer;
    }
    parmOut.atomsPerMol=(int*) F_loadFormat(parmfile,FINT,6,12,parmOut.molecules,debug);
    if (parmOut.atomsPerMol==NULL) {mprintf("Error in atoms per molecule.\n"); return 1;}
    // boxFromParm = {OLDBETA, BOX(1), BOX(2), BOX(3)}
    double *boxFromParm=(double*)  F_loadFormat(parmfile,FDOUBLE,16,5,4,debug);
    if (boxFromParm==NULL) {mprintf("Error in box info.\n"); return 1;}
    parmOut.boxType = SetBoxInfo(boxFromParm,parmOut.Box,debug);
    delete[] boxFromParm;
    if (debug>0) {
      mprintf("\t%s contains box info: %i mols, first solvent mol is %i\n",
              parmOut.parmName, parmOut.molecules, parmOut.firstSolvMol);
      mprintf("\tBOX: %lf %lf %lf | %lf %lf %lf\n",parmOut.Box[0],parmOut.Box[1],parmOut.Box[2],parmOut.Box[3],parmOut.Box[4],parmOut.Box[5]);
      if (parmOut.boxType==ORTHO)
        mprintf("\t     Box is orthogonal.\n");
      else if (parmOut.boxType==NONORTHO)
        mprintf("\t     Box is non-orthogonal.\n");
      else
        mprintf("\t     Box will be determined from first associated trajectory.\n");
    }
  }
  return 0;
}

// AmberParmFile::ReadParmAmber()
int AmberParmFile::ReadParmAmber(AmberParm &parmOut, CpptrajFile &parmfile) {
  int ifbox;
  int *solvent_pointer;
  double *boxFromParm;
  int values[AMBERPOINTERS];
  char *title;
  bool chamber; // true: This topology file is a chamber-created topology file
# ifdef USE_CHARBUFFER
  // TEST: Close and reopen buffered
  parmfile.CloseFile();
  parmfile.OpenFileBuffered();
# endif
  if (debug>0) mprintf("Reading Amber Topology file %s\n",parmfile.filename);
  // Title
  // NOTE: getFlagFileString uses 'new' operator.
  title = getFlagFileString(parmfile, "TITLE",debug);
  // If title is NULL, check for CTITLE (chamber parm)
  if (title==NULL) {
    title = getFlagFileString(parmfile,"CTITLE",debug);
    chamber = true;
  } else {
    chamber = false;
  }
  if (debug>0) mprintf("\tAmberParm Title: %s\n",title);
  delete[] title;
  // Pointers
  int *tempvalues=(int*) getFlagFileValues(parmfile,F_POINTERS,AMBERPOINTERS,debug);
  if (tempvalues==NULL) {
    mprintf("Could not get values from topfile\n");
    return 1;
  }
  memcpy(values, tempvalues, AMBERPOINTERS * sizeof(int));
  delete[] tempvalues;
  // Set some commonly used values
  SetParmFromValues(parmOut,values, false);
  ifbox=values[IFBOX];
  // Atom names
  parmOut.names=(NAME*) getFlagFileValues(parmfile,F_NAMES,parmOut.natom,debug);
  if (parmOut.names==NULL) {mprintf("Error in atom names.\n"); return 1;}
  // Charge; convert to units of electron charge
  parmOut.charge=(double*) getFlagFileValues(parmfile,F_CHARGE,parmOut.natom,debug);
  if (parmOut.charge==NULL) {mprintf("Error in charges.\n"); return 1;}
  for (int atom=0; atom < parmOut.natom; atom++) parmOut.charge[atom] *= (AMBERTOELEC);
  // Mass
  parmOut.mass=(double*) getFlagFileValues(parmfile,F_MASS,parmOut.natom,debug);
  if (parmOut.mass==NULL) {mprintf("Error in masses.\n"); return 1;}
  // Atom type index
  parmOut.atype_index = (int*) getFlagFileValues(parmfile,F_ATYPEIDX,parmOut.natom,debug);
  if (parmOut.atype_index==NULL) {mprintf("Error in atom type index.\n"); return 1;}
  // Number of excluded atoms
  parmOut.numex = (int*) getFlagFileValues(parmfile,F_NUMEX,parmOut.natom,debug);
  if (parmOut.numex==NULL) {mprintf("Error in number of excluded atoms.\n"); return 1;}
  // Nonbonded parm index
  parmOut.NB_index = (int*) getFlagFileValues(parmfile,F_NB_INDEX,parmOut.ntypes*parmOut.ntypes,debug);
  if (parmOut.NB_index==NULL) {mprintf("Error in nonbonded parameter index.\n"); return 1;}
  // Residue names
  parmOut.resnames=(NAME*) getFlagFileValues(parmfile,F_RESNAMES,parmOut.nres,debug);
  if (parmOut.resnames==NULL) {mprintf("Error in residue names.\n"); return 1;}
  // Residue atom #s; shift by -1 so that atom #s start from 0
  parmOut.resnums=(int*) getFlagFileValues(parmfile,F_RESNUMS,parmOut.nres,debug);
  if (parmOut.resnums==NULL) {mprintf("Error in residue numbers.\n"); return 1;}
  for (int res=0; res < parmOut.nres; res++) parmOut.resnums[res] -= 1;
  // Bond force constants and equilibrium values
  parmOut.bond_rk = (double*) getFlagFileValues(parmfile, F_BONDRK, values[NUMBND], debug);
  parmOut.bond_req = (double*) getFlagFileValues(parmfile, F_BONDREQ, values[NUMBND], debug);
  if (parmOut.bond_rk==NULL || parmOut.bond_req==NULL) {mprintf("Error in bond constants.\n"); return 1;}
  // Angle force constants and equilibrium values
  parmOut.angle_tk = (double*) getFlagFileValues(parmfile, F_ANGLETK, values[NUMANG], debug);
  parmOut.angle_teq = (double*) getFlagFileValues(parmfile, F_ANGLETEQ, values[NUMANG], debug);
  if (parmOut.angle_tk==NULL || parmOut.angle_teq==NULL) {mprintf("Error in angle constants.\n"); return 1;}
  // Dihedral force constants, periodicity, and phase values
  parmOut.dihedral_pk = (double*) getFlagFileValues(parmfile, F_DIHPK, values[NPTRA], debug);
  parmOut.dihedral_pn = (double*) getFlagFileValues(parmfile, F_DIHPN, values[NPTRA], debug);
  parmOut.dihedral_phase = (double*) getFlagFileValues(parmfile, F_DIHPHASE, values[NPTRA], debug);
  if (parmOut.dihedral_pk==NULL || parmOut.dihedral_pn==NULL || parmOut.dihedral_phase==NULL) {
    mprintf("Error in dihedral constants.\n"); return 1;
  }
  // SCEE and SCNB scale factors
  parmOut.scee_scale = (double*) getFlagFileValues(parmfile, F_SCEE, values[NPTRA], debug);
  parmOut.scnb_scale = (double*) getFlagFileValues(parmfile, F_SCNB, values[NPTRA], debug);
  // SOLTY: currently unused
  parmOut.solty = (double*) getFlagFileValues(parmfile,F_SOLTY,values[NATYP],debug);
  // Lennard-Jones A/B coefficient
  parmOut.LJ_A = (double*) getFlagFileValues(parmfile,F_LJ_A,parmOut.ntypes*(parmOut.ntypes+1)/2,debug);
  parmOut.LJ_B = (double*) getFlagFileValues(parmfile,F_LJ_B,parmOut.ntypes*(parmOut.ntypes+1)/2,debug);
  if (parmOut.LJ_A==NULL || parmOut.LJ_B==NULL) {mprintf("Error reading LJ parameters.\n"); return 1;}
  // Bond information
  parmOut.bondsh=(int*) getFlagFileValues(parmfile,F_BONDSH,parmOut.NbondsWithH*3,debug);
  parmOut.bonds=(int*) getFlagFileValues(parmfile,F_BONDS,parmOut.NbondsWithoutH*3,debug);
  if (parmOut.bondsh==NULL || parmOut.bonds==NULL) {mprintf("Error in bonds.\n"); return 1;}
  // Angle information
  parmOut.anglesh = (int*) getFlagFileValues(parmfile,F_ANGLESH, values[NTHETH]*4, debug);
  parmOut.angles  = (int*) getFlagFileValues(parmfile,F_ANGLES , values[NTHETA]*4, debug);
  if (parmOut.anglesh==NULL || parmOut.angles==NULL) {mprintf("Error in angles.\n"); return 1;}
  // Dihedral information
  parmOut.dihedralsh = (int*) getFlagFileValues(parmfile,F_DIHH, values[NPHIH]*5,  debug);
  parmOut.dihedrals  = (int*) getFlagFileValues(parmfile,F_DIH , values[NPHIA]*5,  debug);
  if (parmOut.dihedralsh==NULL || parmOut.dihedrals==NULL) {mprintf("Error in dihedrals.\n"); return 1;}
  // List of excluded atoms; shift by -1 so atom #s start from 0
  parmOut.excludedAtoms = (int*) getFlagFileValues(parmfile,F_EXCLUDE,parmOut.nnb,debug);
  if (parmOut.excludedAtoms==NULL) {mprintf("Error reading list of excluded atoms.\n"); return 1;}
  for (int atom=0; atom < parmOut.nnb; atom++) parmOut.excludedAtoms[atom] -= 1;
  // Hbond LJ 10-12 potential terms and cutoff
  parmOut.asol  = (double*) getFlagFileValues(parmfile,F_ASOL, values[NPHB],debug);
  parmOut.bsol  = (double*) getFlagFileValues(parmfile,F_BSOL, values[NPHB],debug);
  parmOut.hbcut = (double*) getFlagFileValues(parmfile,F_HBCUT,values[NPHB],debug);
  // Amber atom types
  parmOut.types=(NAME*) getFlagFileValues(parmfile,F_TYPES,parmOut.natom,debug);
  if (parmOut.types==NULL) {mprintf("Error in atom types.\n"); return 1;}
  // Tree chain classification and joining info 
  parmOut.itree = (NAME*) getFlagFileValues(parmfile,F_ITREE,parmOut.natom,debug);
  parmOut.join_array = (int*) getFlagFileValues(parmfile,F_JOIN,parmOut.natom,debug);
  // Last atom that would move if atom i was rotated; unused
  parmOut.irotat = (int*) getFlagFileValues(parmfile,F_IROTAT,parmOut.natom,debug);
  // Get solvent info if IFBOX>0
  if (values[IFBOX]>0) {
    solvent_pointer=(int*) getFlagFileValues(parmfile,F_SOLVENT_POINTER,3,debug);
    if (solvent_pointer==NULL) {
      mprintf("Error in solvent pointers.\n");
      return 1;
    } else {
      parmOut.finalSoluteRes=solvent_pointer[0];
      parmOut.molecules=solvent_pointer[1];
      parmOut.firstSolvMol=solvent_pointer[2];
      delete[] solvent_pointer;
    }
    parmOut.atomsPerMol=(int*) getFlagFileValues(parmfile,F_ATOMSPERMOL,parmOut.molecules,debug);
    if (parmOut.atomsPerMol==NULL) {mprintf("Error in atoms per molecule.\n"); return 1;}
    // boxFromParm = {OLDBETA, BOX(1), BOX(2), BOX(3)}
    boxFromParm=(double*) getFlagFileValues(parmfile,F_PARMBOX,4,debug);
    // If no box information present in the parm (such as with Chamber prmtops)
    // set the box info if ifbox = 2, otherwise set to NOBOX; the box info will 
    // eventually be set by angles from the first trajectory associated with 
    // this parm.
    if (boxFromParm==NULL) {
      if (not chamber) mprintf("Warning: Prmtop missing Box information.\n");
      // ifbox 2: truncated octahedron for certain
      if (ifbox == 2) {
        parmOut.boxType = NONORTHO;
        parmOut.Box[0] = 0.0; parmOut.Box[1] = 0.0; parmOut.Box[2] = 0.0;
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
    if (debug>0) {
      mprintf("\t%s contains box info: %i mols, first solvent mol is %i\n",
              parmOut.parmName, parmOut.molecules, parmOut.firstSolvMol);
      mprintf("\tBOX: %lf %lf %lf | %lf %lf %lf\n",parmOut.Box[0],parmOut.Box[1],parmOut.Box[2],parmOut.Box[3],parmOut.Box[4],parmOut.Box[5]);
      if (parmOut.boxType==ORTHO)
       mprintf("\t     Box is orthogonal.\n");
      else if (parmOut.boxType==NONORTHO)
        mprintf("\t     Box is non-orthogonal.\n");
      else
        mprintf("\t     Box will be determined from first associated trajectory.\n");
    }
  }
  // GB parameters; radius set, radii, and screening parameters
  parmOut.radius_set = getFlagFileString(parmfile,"RADIUS_SET",debug);
  if (parmOut.radius_set!=NULL) {
    /*radius_set.assign(title);
    delete[] title;
    // Remove whitespace from GB radius set
    radius_set.erase(std::remove(radius_set.begin(), radius_set.end(), ' '), radius_set.end());*/
    if (debug>0) mprintf("\tRadius Set: %s\n",parmOut.radius_set);
  }
  parmOut.gb_radii = (double*) getFlagFileValues(parmfile,F_RADII,parmOut.natom,debug);
  parmOut.gb_screen = (double*) getFlagFileValues(parmfile,F_SCREEN,parmOut.natom,debug);
  if (parmOut.gb_radii==NULL || parmOut.gb_screen==NULL) {mprintf("Error reading gb parameters.\n"); return 1;}


  return 0;
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
  if (parmIn.asol!=NULL)
    DataToFortranBuffer(buffer,F_ASOL,NULL,parmIn.asol,NULL,parmIn.nphb);
  if (parmIn.bsol!=NULL)
    DataToFortranBuffer(buffer,F_BSOL,NULL,parmIn.bsol,NULL,parmIn.nphb);
  if (parmIn.hbcut!=NULL)
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

