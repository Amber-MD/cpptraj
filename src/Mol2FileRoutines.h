// Mol2FileRoutines
#include "CpptrajFile.h"
#define MOL2BUFFERSIZE 256

enum TRIPOSTAG { MOLECULE, ATOM, BOND, SUBSTRUCT };

int Mol2ScanTo( CpptrajFile *, TRIPOSTAG );
int Mol2AtomName(char *, char *);
int Mol2XYZ(char *, double *);
int Mol2AtomType(char *, char *);
int Mol2ResNumName(char *, int *, char *);
double Mol2Charge(char *);
