#ifndef INC_MOL2FILEROUTINES_H
#define INC_MOL2FILEROUTINES_H
/*! \file Mol2FileRoutines.h
    \brief Collection of routines used for accessing mol2 files.
 */
#include "FileIO.h"
#include "Name.h"
// NOTE: Just use buffersize in CpptrajFile?
#define MOL2BUFFERSIZE 256
/// Recognized TRIPOS tags
enum TRIPOSTAG { MOLECULE, ATOM, BOND, SUBSTRUCT };
#define NUMTRIPOSTAGS 4

bool IsMol2Keyword(char *);
int Mol2ScanTo( FileIO *, TRIPOSTAG );
int Mol2AtomName(char *, NAME);
int Mol2XYZ(char *, double *);
int Mol2AtomType(char *, NAME);
int Mol2ResNumName(char *, int *, NAME);
double Mol2Charge(char *);
#endif
