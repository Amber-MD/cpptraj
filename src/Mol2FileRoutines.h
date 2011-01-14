// Mol2FileRoutines
#include "PtrajFile.h"
#define MOL2BUFFERSIZE 256

enum TRIPOSTAG { MOLECULE, ATOM, BOND, SUBSTRUCT };

int Mol2ScanTo( PtrajFile *, TRIPOSTAG );
