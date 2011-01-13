#include <cstring>
#include "Mol2FileRoutines.h"
#include "CpptrajStdio.h"
// Mol2FileRoutines

// Tripos Tags - must be in same order as enum type TRIPOSTAG
const char TRIPOSTAGTEXT[3][20]={"@<TRIPOS>MOLECULE\0",
                                 "@<TRIPOS>ATOM\0",
                                 "@<TRIPOS>BOND\0"
                                };

/*
 * Mol2ScanTo()
 * Scan to the specified TRIPOS section of file
 */
int Mol2ScanTo( PtrajFile *File, TRIPOSTAG tag ) {
  int tagSize;
  char buffer[MOL2BUFFERSIZE];

  tagSize = strlen(TRIPOSTAGTEXT[tag]);
  while ( File->IO->Gets(buffer,MOL2BUFFERSIZE)==0 ) {
    //mprintf("DEBUG: Line [%s]\n",buffer);
    //mprintf("DEBUG: Targ [%s]\n",TRIPOSTAGTEXT[tag]); 
    if (strncmp(buffer,TRIPOSTAGTEXT[tag],tagSize)==0) return 0;
  }

  mprintf("Error: Mol2File::ScanTo(): Could not find tag %s\n",TRIPOSTAGTEXT[tag]);

  return 1;
}
