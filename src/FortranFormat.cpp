#include <cstring>
#include <cstdio> // sprintf
#include <cstdlib> // atof, atoi
#include <cctype> // toupper
#include "FortranFormat.h"
#include "CpptrajStdio.h"

// FFSIZE: Combined size of %FLAG and %FORMAT lines (81 * 2)
#define FFSIZE 162

//  F_POINTERS = 0, F_NAMES,   F_CHARGE,  F_MASS,    F_RESNAMES,             
//  F_RESNUMS,      F_TYPES,   F_BONDSH,  F_BONDS,   F_SOLVENT_POINTER, 
//  F_ATOMSPERMOL,  F_PARMBOX, F_ATYPEIDX,F_NUMEX,   F_NB_INDEX,
//  F_LJ_A,         F_LJ_B,    F_EXCLUDE, F_RADII,   F_SCREEN,
//  F_BONDRK,       F_BONDREQ, F_ANGLETK, F_ANGLETEQ,F_DIHPK,
//  F_DIHPN,        F_DIHPHASE,F_SCEE,    F_SCNB,    F_SOLTY
//  F_ANGLESH,      F_ANGLES,  F_DIHH,    F_DIH,     F_ASOL
//  F_BSOL,         F_HBCUT,   F_ITREE,   F_JOIN,    F_IROTAT

// Constant strings for fortran formats corresponding to Amber parm flags
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

// GetFortranBufferSize()
/** Given number of columns and the width of each column, return the 
  * necessary char buffer size for N data elements.
  */
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
/** If fformat is NULL just report whether flag was found.
  * \return true if flag was found, false if not.
  */
static bool PositionFileAtFlag(CpptrajFile *File, const char *Key, char *fformat, int debug) {
  char lineBuffer[BUFFER_SIZE]; // Hold flag/format line from parmfile
  char value[83];

  if (debug>0) mprintf("Reading %s\n",Key);
  // First, rewind the input file.
  File->IO->Rewind();
  // Search for %FLAG <Key>
  while ( File->IO->Gets(lineBuffer,BUFFER_SIZE) == 0) {
    if ( strncmp(lineBuffer,"%FLAG",5)==0 ) {
      sscanf(lineBuffer,"%*s %s",value);
      if (strcmp(value,Key)==0) {
        if (debug>1) mprintf("DEBUG: Found Flag Key [%s]\n",value);
        // Read next line; can be either a COMMENT or FORMAT. If COMMENT, 
        // read past until you get to the FORMAT line
        File->IO->Gets(lineBuffer,BUFFER_SIZE);
        while (strncmp(lineBuffer,"%FORMAT",7)!=0)
          File->IO->Gets(lineBuffer,BUFFER_SIZE);
        if (debug>1) mprintf("DEBUG: Format line [%s]\n",lineBuffer);
        // Set format
        if (fformat!=NULL) strcpy(fformat, lineBuffer);
        return true;
      } // END found Key
    } // END found FLAG line
  } // END scan through file

  // If we reached here Key was not found.
  if (debug>0)
    mprintf("Warning: [%s] Could not find Key %s in file.\n",File->filename,Key);
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
char *getFlagFileString(CpptrajFile *File, const char *Key, int debug) {
  char *lineBuffer;

  // Find flag
  if (PositionFileAtFlag(File,Key,NULL,debug)) {
    // Read next line and return
    lineBuffer = new char[83]; // 80 + newline + NULL ( + CR if dos)
    File->IO->Gets(lineBuffer,82);
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
void *getFlagFileValues(CpptrajFile *File, AmberParmFlagType fflag, int maxval, int debug){
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
  BufferSize = GetFortranBufferSize(maxval, File->isDos, width, ncols);
  buffer = new char[ BufferSize ];
  if ( File->IO->Read(buffer,sizeof(char),BufferSize)==-1 ) {
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
char *F_load20a4(CpptrajFile *File) {
  char *lineBuffer = new char[83]; // 80 + newline + NULL ( + CR if dos)
  File->IO->Gets(lineBuffer,82);
  RemoveWhitespace(lineBuffer);
  return lineBuffer;
}

// F_loadFormat()
void *F_loadFormat(CpptrajFile *File, FortranType fType, int width, int ncols, 
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
    File->IO->Gets(temp,16);
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
  BufferSize = GetFortranBufferSize(maxval, File->isDos, width, ncols);
  buffer = new char[ BufferSize ];
  if ( File->IO->Read(buffer,sizeof(char),BufferSize)==-1 ) {
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
int DataToFortranBuffer(CharBuffer &buffer, AmberParmFlagType fFlag,
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

