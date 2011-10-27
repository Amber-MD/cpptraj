#include "FortranFormat.h"
#include "CpptrajStdio.h"
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <cctype> // toupper

//  F_POINTERS = 0, F_NAMES,  F_CHARGE, F_MASS,   F_RESNAMES,             
//  F_RESNUMS,      F_TYPES,  F_BONDSH, F_BONDS,  F_SOLVENT_POINTER, 
//  F_ATOMSPERMOL,  F_PARMBOX

// Constant strings for fortran formats corresponding to Amber parm flags
static const char AmberParmFmt[NUMAMBERPARMFLAGS][16] = {
"%FORMAT(10I8)",  "%FORMAT(20a4)", "%FORMAT(5E16.8)", "%FORMAT(5E16.8)","%FORMAT(20a4)",
"%FORMAT(10I8)",  "%FORMAT(20a4)", "%FORMAT(10I8)",   "%FORMAT(10I8)",  "%FORMAT(3I8)",
"%FORMAT(10I8)",  "%FORMAT(5E16.8)"
}; 
static const char AmberParmFlag[NUMAMBERPARMFLAGS][23] = {
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
  "BOX_DIMENSIONS"
};

// Enumerated type for Fortran data type
enum FortranType {
  UNKNOWN_FTYPE, FINT, FDOUBLE, FCHAR, FFLOAT
};

/* GetFortranBufferSize()
 * Given number of columns and the width of each column, return the 
 * necessary char buffer size for N data elements.
 */
int GetFortranBufferSize(int N, int isDos, int width, int ncols) {
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

/* GetFortranType()
 * Given a fortran-type format string, return the corresponding fortran
 * type. Set ncols (if present), width, and precision (if present).
 * 01234567
 * %FORMAT([<cols>][(]<type><width>[<precision>][)])
 */
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

/* getFlagFileString()
 * Search for the FLAG specified by Key. Assume the next line is a string 
 * of max length 80 chars and return it.
 */
char *getFlagFileString(CpptrajFile *File, const char *Key, int debug) {
  char *lineBuffer = new char[83]; // 80 + newline + NULL ( + CR if dos)
  char value[83];

  if (debug>0) mprintf("Reading %s\n",Key);

  // First, rewind the input file.
  File->IO->Rewind();

  // Next, search for the required FLAG
  while ( File->IO->Gets(lineBuffer,82) == 0) {
    if ( strncmp(lineBuffer,"%FLAG",5)==0 ) {
      sscanf(lineBuffer,"%*s %s",value);
      if (strcmp(value,Key)==0) {
        if (debug>1) mprintf("DEBUG: Found Flag Key [%s]\n",value);
        // Read next line; can be either a COMMENT or FORMAT. If COMMENT, 
        // read past until you get to the FORMAT line
        File->IO->Gets(lineBuffer,82);
        while (strncmp(lineBuffer,"%FORMAT",7)!=0)
          File->IO->Gets(lineBuffer,82);
        if (debug>1) mprintf("DEBUG: Format line [%s]\n",lineBuffer);
        // Read next line and return
        File->IO->Gets(lineBuffer,82);
        return lineBuffer;
      }
    }
  }

  // If here, Key not found, could be bad news, but let the calling function
  // print the error message.
  if (debug>0)
    mprintf("Warning: [%s] Could not find Key %s in file.\n",File->filename,Key);
  delete[] lineBuffer;
  return NULL;
}

/* getFlagFileValues()
 * Search for the FLAG specified by Key and return the values. The values will 
 * be put into an array type according to the FORMAT string in the top file, 
 * but it is necessary to explictly type the returned array. maxval is used to 
 * allocate memory for the return array - only maxval values will be read.
 */
void *getFlagFileValues(CpptrajFile *File, const char *Key, int maxval, int debug){
  int i, ncols, width, precision; 
  char lineBuffer[BUFFER_SIZE]; // Hold flag/format line from parmfile
  char value[83];      // Hold Key from Flag line
  char temp[17];       // Hold data element
  char *buffer,*ptr;
  NAME *C;
  int *I;
  double *D;
  FortranType fType;
  int BufferSize;

  if (debug>0) {
    mprintf("Reading %s\n",Key);
    if (debug>1) mprintf("DEBUG: maxval= %i\n",maxval);
  }

  C=NULL; D=NULL; I=NULL;
  
  // First, rewind the input file.
  File->IO->Rewind();

  // Next, search for the required FLAG
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
        fType = GetFortranType(lineBuffer, &ncols, &width, &precision);
        if (debug>1) mprintf("DEBUG: Format type %i\n",(int)fType);
        if (fType == UNKNOWN_FTYPE) {
          mprinterr("Error: Unrecognized fortran format [%s] for key [%s]\n",lineBuffer,Key);
          return NULL;
        } 
        // Allocate memory based on data type
        switch (fType) {
          case UNKNOWN_FTYPE : return NULL;
          case FINT   : I=(int*)    malloc(maxval*sizeof(int)); break;
          case FDOUBLE: D=(double*) malloc(maxval*sizeof(double)); break;
          case FCHAR  : C=(NAME*)   malloc(maxval*sizeof(NAME)); break;
          case FFLOAT : D=(double*) malloc(maxval*sizeof(double)); break;
        }
        // Allocate memory to read in entire section
        BufferSize = GetFortranBufferSize(maxval, File->isDos, width, ncols);
        buffer=(char*) calloc(BufferSize,sizeof(char));
        if ( File->IO->Read(buffer,sizeof(char),BufferSize)==-1 ) {
          rprinterr("ERROR in read of prmtop section %s\n",Key);
          free(buffer);
          break; // Send us outside the while loop
        }
        if (debug>3) mprintf("DEBUG: Buffer [%s]\n",buffer);

        // Convert values in buffer to their type
        ptr=buffer;
        temp[width]='\0';          
        for (i=0; i<maxval; i++) {
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
            if (I!=NULL) free(I);
            if (D!=NULL) free(D);
            if (C!=NULL) free(C);
            free(buffer);
            return NULL;
          }
          /* Now convert temp to appropriate type */
          if (I!=NULL) I[i]=atoi(temp);
          if (D!=NULL) D[i]=atof(temp);
          if (C!=NULL) strcpy(C[i],temp);
          ptr+=width;
        }
        /* All done! */
        free(buffer);
        if (I!=NULL) return I;
        if (D!=NULL) return D;
        if (C!=NULL) return C;
      } // End if (strcmp(value,Key)==0)
    } // End if ( strncmp(lineBuffer,"%FLAG",5)==0 )
  } // End While loop
  // If here, Key not found, could be bad news, but let the calling function
  // print the error message.
  if (debug>0)
    mprintf("Warning: [%s] Could not find Key %s in file.\n",File->filename,Key); 
  if (I!=NULL) free(I);
  if (D!=NULL) free(D);
  if (C!=NULL) free(C);

  return NULL;
}

/* DataToFortranBuffer()
 * Write N data elements stored in I, D, or C to buffer with given 
 * fortran format.
 */
char *DataToFortranBuffer(char *bufferIn, AmberParmFlagType fFlag,
                          int *I, double *D, NAME *C, int N) 
{
  int coord, width, numCols, precision;
  char *ptr;
  FortranType fType;
  char FormatString[32];


  // Determine type, cols, width, and precision from format string
  fType = GetFortranType((char*)AmberParmFmt[fFlag], &numCols, &width, &precision);
  if (fType == UNKNOWN_FTYPE) {
    mprinterr("Error: DataToFortranBuffer: Unknown format string [%s]\n",AmberParmFmt[fFlag]);
    return NULL;
  }

  if (bufferIn==NULL) return NULL;
  ptr = bufferIn;

  //fprintf(stdout,"*** Called DataToBuffer: N=%i, width=%i, numCols=%i, format=%s\n",
  //        N, width, numCols, FormatString);
  //fprintf(stdout,"*** Buffer address is %p\n",buffer);

  // Print FLAG and FORMAT lines
  // '%FLAG '
  //  012345
  sprintf(ptr,"%%FLAG %-74s\n",AmberParmFlag[fFlag]);
  //sprintf(ptr,"%-80s\n",FLAG);
  ptr += 81;
  sprintf(ptr,"%-80s\n",AmberParmFmt[fFlag]);
  ptr += 81;

  // Write Integer data
  // NOTE: move the coord+1 code?
  coord=0;
  if (fType==FINT) {
    if (I==NULL) {
      mprinterr("Error: DataToFortranBuffer: INT is NULL.\n");
      return NULL;
    }
    sprintf(FormatString,"%%%ii",width);
    for (coord=0; coord < N; coord++) {
      sprintf(ptr,FormatString,I[coord]);
      ptr+=width;
      if ( ((coord+1)%numCols)==0 ) {
        sprintf(ptr,"\n");
        ptr++;
      }
    }

  // Write Double data
  } else if (fType==FDOUBLE) {
    if (D==NULL) {
      mprinterr("Error: DataToFortranBuffer: DOUBLE is NULL.\n");
      return NULL;
    }
    sprintf(FormatString,"%%%i.%ilE",width,precision);
    for (coord=0; coord < N; coord++) {
      sprintf(ptr,FormatString,D[coord]);
      ptr+=width;
      if ( ((coord+1)%numCols)==0 ) {
        sprintf(ptr,"\n");
        ptr++;
      }
    }

  // Write Char data
  } else if (fType==FCHAR) {
    if (C==NULL) {
      mprinterr("Error: DataToFortranBuffer: CHAR is NULL.\n");
      return NULL;
    }
    sprintf(FormatString,"%%%is",width);
    for (coord=0; coord < N; coord++) {
      sprintf(ptr,FormatString,C[coord]);
      ptr+=width;
      if ( ((coord+1)%numCols)==0 ) {
        sprintf(ptr,"\n");
        ptr++;
      }
    }

  // Write Float data
  } else if (fType==FFLOAT) {
    if (D==NULL) {
      mprinterr("Error: DataToFortranBuffer: FLOAT is NULL.\n");
      return NULL;
    }
    sprintf(FormatString,"%%%i.%ilf",width,precision);
    for (coord=0; coord < N; coord++) {
      sprintf(ptr,FormatString,D[coord]);
      ptr+=width;
      if ( ((coord+1)%numCols)==0 ) {
        sprintf(ptr,"\n");
        ptr++;
      }
    }
  }

  // If the coord record didnt end on a newline, print one
  if ( (coord%numCols)!=0 ) {
    sprintf(ptr,"\n");
    ptr++; // Only needed if more will be written
  }

  //fprintf(stdout,"*** Ptr ended on %p\n",ptr);
  return ptr;
}

