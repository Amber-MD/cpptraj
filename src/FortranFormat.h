#ifndef INC_FORTRANFORMAT_H
#define INC_FORTRANFORMAT_H
/*! \file FortranFormat.h
    \brief Routines that pertain to reading and writing files with
           fortran FLAG and FORMAT keywords (currently only Amber Topology files).
 */
#include "CpptrajFile.h"
#include "Name.h"
#include "CharBuffer.h"
/// Enumerated type for Amber Parmtop Flags
enum AmberParmFlagType {
  F_POINTERS = 0, F_NAMES,   F_CHARGE,  F_MASS,    F_RESNAMES,    
  F_RESNUMS,      F_TYPES,   F_BONDSH,  F_BONDS,   F_SOLVENT_POINTER, 
  F_ATOMSPERMOL,  F_PARMBOX, F_ATYPEIDX,F_NUMEX,   F_NB_INDEX,
  F_LJ_A,         F_LJ_B,    F_EXCLUDE, F_RADII,   F_SCREEN,
  F_BONDRK,       F_BONDREQ, F_ANGLETK, F_ANGLETEQ,F_DIHPK,
  F_DIHPN,        F_DIHPHASE,F_SCEE,    F_SCNB,    F_SOLTY,
  F_ANGLESH,      F_ANGLES,  F_DIHH,    F_DIH,     F_ASOL,
  F_BSOL,         F_HBCUT,   F_ITREE,   F_JOIN,    F_IROTAT
};
#define NUMAMBERPARMFLAGS 40 
/// Enumerated type for Fortran data type
enum FortranType {
  UNKNOWN_FTYPE, FINT, FDOUBLE, FCHAR, FFLOAT
};
char *getFlagFileString(CpptrajFile &, const char *, int);
void *getFlagFileValues(CpptrajFile &,AmberParmFlagType,int,int);
char *F_load20a4(CpptrajFile &);
void *F_loadFormat(CpptrajFile &, FortranType, int, int, int, int);
int DataToFortranBuffer(CharBuffer&,AmberParmFlagType,int*,double*,NAME*,int);
#endif
