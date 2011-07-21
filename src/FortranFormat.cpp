#include "FortranFormat.h"
#include "CpptrajStdio.h"
#include <cstring>
#include <cstdio>
#include <cstdlib>
// Enumerated type for Fortran data type
enum FortranType {
  UNKNOWN_FTYPE, FINT, FDOUBLE, FCHAR
};
// Static arrays containing data for each Fortran Format
const FortranType FortranFType[NUMFORTRANFORMAT]  = {
  UNKNOWN_FTYPE, FINT, FDOUBLE, FCHAR, FINT, FINT
};
const int FortranFWidth[NUMFORTRANFORMAT] = {0,  8, 16,  4,  6, 8};
const int FortranFCols[NUMFORTRANFORMAT]  = {0, 10,  5, 20, 12, 3};
const char FortranFString[NUMFORTRANFORMAT][8] = {
  "\0", "%8i", "%16.8lE", "%4s", "%6i", "%8i"
};
const char FortranFOutString[NUMFORTRANFORMAT][16] = {
  "\0", "%FORMAT(10I8)", "%FORMAT(5E16.8)", "%FORMAT(20a4)",
  "%FORMAT(12I6)", "%FORMAT(3I8)"
};

/* GetFortranBufferSize()
 * Given a fortran format, return the necessary char buffer size for
 * N data elements of that format.
 */
int GetFortranBufferSize(FortranFormat fFormat, int N, int isDos) {
  int bufferLines=0;
  int BufferSize=0;

  BufferSize= N * FortranFWidth[fFormat];
  bufferLines = N / FortranFCols[fFormat];
  if ((N % FortranFCols[fFormat])!=0) bufferLines++;
  // If DOS file there are CR before Newlines
  if (isDos) bufferLines*=2;
  BufferSize+=bufferLines;
  //if (debug>0) 
  //  fprintf(stdout,"*** Buffer size is %i including %i newlines.\n",BufferSize,bufferLines);
  return BufferSize;
}

/* GetFortranFormat()
 * Given a fortran-type format string, return the corresponding fortran
 * type. 
 */
static FortranFormat GetFortranFormat(char *Format) {
  if        ( strncmp(Format,"%FORMAT(10I8)"  ,13)==0 ) {
    return F10I8;
  } else if ( strncmp(Format,"%FORMAT(5E16.8)",15)==0 ) {
    return F5E16_8;
  } else if ( strncmp(Format,"%FORMAT(20a4)"  ,13)==0 ) {
    return F20a4;
  } else if ( strncmp(Format,"%FORMAT(12I6)",  13)==0 ) {
    return F12I6;
  } else if ( strncmp(Format,"%FORMAT(3I8)"   ,12)==0 ) {
    return F3I8;
  } else {
    mprinterr("Error: Unrecognized fortran format [%s]\n",Format);
    return UNKNOWN_FFORMAT;
  }
}

/* getFlagFileValues()
 * Search for the FLAG specified by Key and return the values. The values will 
 * be put into an array type according to the FORMAT string in the top file, 
 * but it is necessary to explictly type the returned array. maxval is used to 
 * allocate memory for the return array - only maxval values will be read.
 */
void *getFlagFileValues(PtrajFile *File, const char *Key, int maxval, int debug){
  int i; 
  char lineBuffer[BUFFER_SIZE]; // Hold flag/format line from parmfile
  char value[83];      // Hold Key from Flag line
  char temp[17];       // Hold data element
  char *buffer,*ptr;
  NAME *C;
  int *I;
  double *D;

  FortranFormat fFormat;
  FortranType fType;
  int BufferSize;
  int width;

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
        fFormat = GetFortranFormat(lineBuffer);
        if (debug>1) mprintf("DEBUG: Format type %i\n",fFormat);
        if (fFormat == UNKNOWN_FFORMAT) return NULL;
        fType = FortranFType[fFormat];
        width = FortranFWidth[fFormat];
        // Allocate memory based on data type
        switch (fType) {
          case UNKNOWN_FTYPE : return NULL;
          case FINT   : I=(int*)    malloc(maxval*sizeof(int)); break;
          case FDOUBLE: D=(double*) malloc(maxval*sizeof(double)); break;
          case FCHAR  : C=(NAME*)   malloc(maxval*sizeof(NAME)); break;
        }
        // Allocate memory to read in entire section
        BufferSize = GetFortranBufferSize(fFormat, maxval, File->isDos);
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
  // If we have scanned through the input file and have not found Key, bad! 
  rprintf("Error: Could not find key %s in file.\n",Key);
  if (I!=NULL) free(I);
  if (D!=NULL) free(D);
  if (C!=NULL) free(C);

  return NULL;
}

/* DataToFortranBuffer()
 * Write N data elements stored in I, D, or C to buffer with given 
 * fortran format.
 */
char *DataToFortranBuffer(char *bufferIn, const char *FLAG, FortranFormat fFormat, 
                          int *I, double *D, NAME *C, int N) {
  int coord, width, numCols;
  char *ptr;
  FortranType fType;
  const char *FormatString;

  width = FortranFWidth[fFormat];
  numCols=FortranFCols[fFormat];
  fType = FortranFType[fFormat];
  FormatString = FortranFString[fFormat];

  if (bufferIn==NULL) return NULL;
  ptr = bufferIn;

  //fprintf(stdout,"*** Called DataToBuffer: N=%i, width=%i, numCols=%i, format=%s\n",
  //        N, width, numCols, FormatString);
  //fprintf(stdout,"*** Buffer address is %p\n",buffer);

  // Print FLAG and FORMAT lines
  sprintf(ptr,"%-80s\n",FLAG);
  ptr += 81;
  sprintf(ptr,"%-80s\n",FortranFOutString[fFormat]);
  ptr += 81;
/*
  for (coord=0; coord<N; coord++) {
    if      (I!=NULL) sprintf(ptr,FormatString,I[coord]);
    else if (D!=NULL) sprintf(ptr,FormatString,D[coord]);
    else if (C!=NULL) sprintf(ptr,FormatString,C[coord]);
    ptr+=width;
    if ( ((coord+1)%numCols)==0 ) {
      sprintf(ptr,"\n");
      ptr++;
    }
  } 
  // If the coord record didnt end on a newline, print one
  if ( (coord%numCols)!=0 ) {
    sprintf(ptr,"\n");
    ptr++; // Only needed if more will be written
  } 
  
  //fprintf(stdout,"*** Ptr ended on %p\n",ptr);
  return ptr;
*/ 
  // Write Integer data
  // NOTE: move the coord+1 code?
  coord=0;
  if (fType==FINT) {
    if (I==NULL) {
      mprinterr("Error: DataToFortranBuffer: INT is NULL.\n");
      return NULL;
    }
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
    for (coord=0; coord < N; coord++) {
      sprintf(ptr,FormatString,C[coord]);
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

