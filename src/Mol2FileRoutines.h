#ifndef INC_MOL2FILEROUTINES_H
#define INC_MOL2FILEROUTINES_H
/*! \file Mol2FileRoutines.h
    \brief Collection of routines used for accessing mol2 files.
 */
#include "CpptrajFile.h"
#include "Name.h"
// NOTE: Just yse buffersize in CpptrajFile?
#define MOL2BUFFERSIZE 256
/// Recognized TRIPOS tags
enum TRIPOSTAG { MOLECULE, ATOM, BOND, SUBSTRUCT };
#define NUMTRIPOSTAGS 4

int Mol2ScanTo( CpptrajFile *, TRIPOSTAG );
int Mol2AtomName(char *, NAME);
int Mol2XYZ(char *, double *);
int Mol2AtomType(char *, NAME);
int Mol2ResNumName(char *, int *, NAME);
double Mol2Charge(char *);
#endif
