#ifndef INC_PARM_AMBER_H
#define INC_PARM_AMBER_H
#include "ParmIO.h"
class Parm_Amber : public ParmIO {
  public :
    Parm_Amber();
    ~Parm_Amber();
    bool ID_ParmFormat();
    int ReadParm(Topology&);
    int WriteParm(Topology&);
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
      F_ATOMICNUM,    F_TITLE,   F_CTITLE,  F_RADSET
    };
    static const int NUMAMBERPARMFLAGS;
    static const int AMBERPOINTERS;
    //static const size_t FFSIZE;
    static const char AmberParmFmt[][16];
    static const char AmberParmFlag[][27];

    static const size_t BUF_SIZE = 83;
    char lineBuffer_[BUF_SIZE];
    bool newParm_;
    std::string fformat_;
    FortranType ftype_;
    int fncols_;
    int fprecision_;
    int fwidth_;
    int error_count_;
    char *buffer_;
    size_t buffer_size_;
    size_t buffer_max_size_;

    int ReadParmOldAmber(Topology&);
    int ReadParmAmber(Topology&);

    //std::string GetFlagLine(AmberParmFlagType);
    std::string GetLine();
    std::vector<int> GetInteger(int,int,int);
    std::vector<double> GetDouble(int,int,int);
    std::vector<NameType> GetName(int,int,int);
    std::vector<int> GetFlagInteger(AmberParmFlagType,int);
    std::vector<double> GetFlagDouble(AmberParmFlagType,int);
    std::vector<NameType> GetFlagName(AmberParmFlagType,int);
    bool SeekToFlag(AmberParmFlagType);
    int AllocateAndRead(int,int,int);
    bool PositionFileAtFlag(AmberParmFlagType);

    int WriteSetup(AmberParmFlagType,size_t);
    int WriteInteger(AmberParmFlagType,std::vector<int>&);
    int WriteDouble(AmberParmFlagType,std::vector<double>&);
    int WriteName(AmberParmFlagType,std::vector<NameType>&);

    size_t GetFortranBufferSize(int,int,int);
    bool SetFortranType();
};
#endif
