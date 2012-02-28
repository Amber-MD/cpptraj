#ifndef INC_PARM_AMBER_H
#define INC_PARM_AMBER_H
#include "ParmIO.h"
#include "CharBuffer.h" // NOTE: Get rid of this eventually?
class AmberParmFile : public ParmIO {
  public :
    AmberParmFile();
    ~AmberParmFile();
    int ReadParm(AmberParm&, CpptrajFile&);
    int WriteParm(AmberParm&, CpptrajFile&);

  private :
    /// Enumerated type for Fortran data type
    enum FortranType {
      UNKNOWN_FTYPE, FINT, FDOUBLE, FCHAR, FFLOAT
    };
    /// Enumerated type for Amber Parmtop Flags
    enum AmberParmFlagType {
      F_POINTERS = 0, F_NAMES,   F_CHARGE,  F_MASS,    F_RESNAMES,
      F_RESNUMS,      F_TYPES,   F_BONDSH,  F_BONDS,   F_SOLVENT_POINTER,
      F_ATOMSPERMOL,  F_PARMBOX, F_ATYPEIDX,F_NUMEX,   F_NB_INDEX,
      F_LJ_A,         F_LJ_B,    F_EXCLUDE, F_RADII,   F_SCREEN,
      F_BONDRK,       F_BONDREQ, F_ANGLETK, F_ANGLETEQ,F_DIHPK,
      F_DIHPN,        F_DIHPHASE,F_SCEE,    F_SCNB,    F_SOLTY,
      F_ANGLESH,      F_ANGLES,  F_DIHH,    F_DIH,     F_ASOL,
      F_BSOL,         F_HBCUT,   F_ITREE,   F_JOIN,    F_IROTAT,
      F_ATOMICNUM
    };
    static const int NUMAMBERPARMFLAGS;
    static const int AMBERPOINTERS;
    static const size_t FFSIZE;

    static const char AmberParmFmt[][16];
    static const char AmberParmFlag[][27];

    CpptrajFile *File;
    char *buffer;
    int error_count; 

    FortranType GetFortranType(char *, int *, int *, int *);
    int DataToFortranBuffer(CharBuffer &, AmberParmFlagType, int *, double *, NAME *, int);
    int AllocateAndRead(int &, int &, int &);
    double *GetDouble(int, int, int);
    int *GetInteger(int, int, int);
    NAME *GetName(int, int, int);
    char *GetLine();
    char *GetFlagLine(const char*);
    FortranType SeekToFlag(AmberParmFlagType fflag, int &, int &, int &);
    double *GetFlagDouble(AmberParmFlagType, int);
    int *GetFlagInteger(AmberParmFlagType, int);
    NAME *GetFlagName(AmberParmFlagType, int);

    void SetParmFromValues(AmberParm &, int *, bool);
    int ReadParmAmber(AmberParm&, CpptrajFile&);
    int ReadParmOldAmber(AmberParm&, CpptrajFile&);
};
#endif
