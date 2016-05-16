#ifndef INC_PARM_AMBER_H
#define INC_PARM_AMBER_H
#include "ParmIO.h"
#include "BufferedFrame.h"
class Parm_Amber : public ParmIO {
  public :
    Parm_Amber();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new Parm_Amber(); }
    static void WriteHelp();
    // ----- Inherited functions -----------------
    bool ID_ParmFormat(CpptrajFile&);
    int processReadArgs(ArgList&) { return 0; }
    int ReadParm(FileName const&, Topology&);
    int processWriteArgs(ArgList&);
    int WriteParm(FileName const&, Topology const&);
  private :
    typedef std::vector<double> Darray;
    typedef std::vector<int> Iarray;
    /// Class for determining data size from Fortran format string.
    class FortranData;
    /// Enumerated type for Amber Parmtop Flags
    enum AmberParmFlagType {
      F_POINTERS = 0, F_NAMES,     F_CHARGE,    F_MASS,     F_RESNAMES,
      F_RESNUMS,      F_TYPES,     F_BONDSH,    F_BONDS,    F_SOLVENT_POINTER,
      F_ATOMSPERMOL,  F_PARMBOX,   F_ATYPEIDX,  F_NUMEX,    F_NB_INDEX,
      F_LJ_A,         F_LJ_B,      F_EXCLUDE,   F_RADII,    F_SCREEN,
      F_BONDRK,       F_BONDREQ,   F_ANGLETK,   F_ANGLETEQ, F_DIHPK,
      F_DIHPN,        F_DIHPHASE,  F_SCEE,      F_SCNB,     F_SOLTY,
      F_ANGLESH,      F_ANGLES,    F_DIHH,      F_DIH,      F_ASOL,
      F_BSOL,         F_HBCUT,     F_ITREE,     F_JOIN,     F_IROTAT,
      F_ATOMICNUM,    F_TITLE,     F_RADSET,    F_LES_NTYP, F_LES_TYPE,
      F_LES_FAC,      F_LES_CNUM,  F_LES_ID,    F_CAP_INFO, F_CAP_INFO2,
      F_IPOL,         F_POLAR,     F_CTITLE,    F_CHM_UBC,  F_CHM_UB,
      F_CHM_UBFC,     F_CHM_UBEQ,  F_CHM_NIMP,  F_CHM_IMP,  F_CHM_NIMPT,
      F_CHM_IMPFC,    F_CHM_IMPP,  F_LJ14A,     F_LJ14B,    F_CHM_CMAPC,
      F_CHM_CMAPR,    F_CHM_CMAPP, F_CHM_CMAPI, F_FF_TYPE,  F_PDB_RES,
      F_PDB_CHAIN,    F_PDB_ICODE, F_PDB_ALT
    };
    /// Used to hold %FLAG/FORMAT string pairs. Corresponds to AmberParmFlagType.
    struct ParmFlag {
      const char* Flag; ///< %FLAG name in topology.
      const char* Fmt;  ///< Fortran format string for writing.
    };
    typedef const ParmFlag* ParmPtr;
    // NOTE: Although amber topology files should only ever require 83 chars
    //       to read each line (80 chars + newline + CR (if dos) + NULL)
    //       increase the size to handle non-standard files.
    static const size_t BUF_SIZE = 256;
    /// Amber topology subtype
    enum ParmType { OLDPARM = 0, NEWPARM, CHAMBER };

    int ReadOldParm(Topology&);
    int ReadNewParm(Topology&);
    int ReadFormatLine(FortranData&);
    void ReadTitle(Topology&);
    int ReadPointers(int, FortranData&);
 
    static const int AMBERPOINTERS_;
    static const ParmFlag FLAGS_[];

    ParmType ptype_;
    BufferedFrame infile_;
};
// -----------------------------------------------------------------------------
class Parm_Amber::FortranData {
  public:
    /// Enumerated type for Fortran data type
    enum Type { UNKNOWN_FTYPE=0, FINT, FDOUBLE, FCHAR, FFLOAT };
    FortranData() : ftype_(UNKNOWN_FTYPE), fncols_(0), fwidth_(0), fprecision_(0) {}
    FortranData(const char*);
    int ParseFortranFormat(const char*);

    Type Ftype()    const { return ftype_; }
    int Ncols()     const { return fncols_; }
    int Width()     const { return fwidth_; }
    int Precision() const { return fprecision_; }
  private:
    Type ftype_;
    int fncols_;
    int fwidth_;
    int fprecision_;
};
#endif
