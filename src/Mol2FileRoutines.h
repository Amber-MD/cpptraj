// Mol2FileRoutines
#include "CpptrajFile.h"
#include "Name.h"
// Note: Just yse buffersize in CpptrajFile?
#define MOL2BUFFERSIZE 256

enum TRIPOSTAG { MOLECULE, ATOM, BOND, SUBSTRUCT };
#define NUMTRIPOSTAGS 4

int Mol2ScanTo( CpptrajFile *, TRIPOSTAG );
int Mol2AtomName(char *, NAME);
int Mol2XYZ(char *, double *);
int Mol2AtomType(char *, NAME);
int Mol2ResNumName(char *, int *, NAME);
double Mol2Charge(char *);
