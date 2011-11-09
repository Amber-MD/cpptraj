// FortranFormat
// This file contains routines that pertain to reading and writing files with
// fortran FLAG and FORMAT keywords (currently only Amber Topology files).
#include "CpptrajFile.h"
#include "Name.h"
// Enumerated type for Amber Parmtop Flags
enum AmberParmFlagType {
  F_POINTERS = 0, F_NAMES,  F_CHARGE,  F_MASS,   F_RESNAMES,    
  F_RESNUMS,      F_TYPES,  F_BONDSH,  F_BONDS,  F_SOLVENT_POINTER, 
  F_ATOMSPERMOL,  F_PARMBOX,F_ATYPEIDX,F_NUMEX,  F_NB_INDEX,
  F_LJ_A,         F_LJ_B,   F_EXCLUDE, F_RADII,  F_SCREEN
};
#define NUMAMBERPARMFLAGS 20 
// Enumerated type for Fortran data type
enum FortranType {
  UNKNOWN_FTYPE, FINT, FDOUBLE, FCHAR, FFLOAT
};
// FFSIZE: Combined size of %FLAG and %FORMAT lines (81 * 2)
#define FFSIZE 162
int GetFortranBufferSize(int,int,int,int);
char *getFlagFileString(CpptrajFile *, const char *, int);
void *getFlagFileValues(CpptrajFile *,AmberParmFlagType,int,int);
char *F_load20a4(CpptrajFile *);
void *F_loadFormat(CpptrajFile *, FortranType, int, int, int, int);
char *DataToFortranBuffer(char*,AmberParmFlagType,int*,double*,NAME*,int);
