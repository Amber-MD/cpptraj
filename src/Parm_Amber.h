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
    /// Enumerated type for Fortran data type
    enum Type { UNKNOWN_FTYPE=0, FINT, FDOUBLE, FCHAR, FFLOAT };
    /// Enumerated type for Amber Parmtop Flags. KEEP IN SYNC WITH FLAGS_ ARRAY
    enum FlagType {
      F_POINTERS = 0, F_NAMES,     F_CHARGE,    F_MASS,      F_RESNAMES,
      F_RESNUMS,      F_TYPES,     F_BONDSH,    F_BONDS,     F_SOLVENT_POINTER,
      F_ATOMSPERMOL,  F_PARMBOX,   F_ATYPEIDX,  F_NUMEX,     F_NB_INDEX,
      F_LJ_A,         F_LJ_B,      F_LJ_C,      F_EXCLUDE,   F_RADII,
      F_SCREEN,       F_BONDRK,    F_BONDREQ,   F_ANGLETK,   F_ANGLETEQ,
      F_DIHPK,        F_DIHPN,     F_DIHPHASE,  F_SCEE,      F_SCNB,
      F_SOLTY,        F_ANGLESH,   F_ANGLES,    F_DIHH,      F_DIH,
      F_ASOL,         F_BSOL,      F_HBCUT,     F_ITREE,     F_JOIN,
      F_IROTAT,       F_ATOMICNUM, F_TITLE,     F_RADSET,    F_LES_NTYP,
      F_LES_TYPE,     F_LES_FAC,   F_LES_CNUM,  F_LES_ID,    F_CAP_INFO,
      F_CAP_INFO2,    F_IPOL,      F_POLAR,     F_CTITLE,    F_CHM_UBC,
      F_CHM_UB,       F_CHM_UBFC,  F_CHM_UBEQ,  F_CHM_NIMP,  F_CHM_IMP,
      F_CHM_NIMPT,    F_CHM_IMPFC, F_CHM_IMPP,  F_LJ14A,     F_LJ14B,
      F_CHM_CMAPC,    F_CHM_CMAPR, F_CHM_CMAPP, F_CHM_CMAPI, F_FF_TYPE,
      F_PDB_RES,      F_PDB_CHAIN, F_PDB_ICODE, F_PDB_ALT,   F_PDB_BFAC,
      F_PDB_OCC,      F_PDB_NUM,   F_CMAPC,     F_CMAPR,     F_CMAPP,
      F_CMAPI
    };
    /// Used to hold %FLAG/FORMAT string pairs. Corresponds to FlagType.
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
    inline const char* SkipToNextFlag();

    void ResetFileToFlag(FlagType);
    void ProblemFlagWarning(FlagType, unsigned int, unsigned int);
    double FileBufferToDouble(FlagType, unsigned int, unsigned int);

    int ReadTitle(Topology&);
    int ReadPointers(int, Topology&, FortranData const&);
    inline int SetupBuffer(FlagType, int, FortranData const&);
    int ReadAtomNames(Topology&, FortranData const&);
    int ReadAtomCharges(Topology&, FortranData const&);
    int ReadAtomicNum(FortranData const&);
    int ReadAtomicMass(Topology&, FortranData const&);
    int ReadAtomTypeIndex(Topology&, FortranData const&);
    int ReadNonbondIndices(Topology&, FortranData const&);
    int ReadResidueNames(Topology&, FortranData const&);
    int ReadResidueAtomNums(Topology&, FortranData const&);
    int ReadBondRK(Topology&, FortranData const&);
    int ReadBondREQ(Topology&, FortranData const&);
    int ReadAngleTK(Topology&, FortranData const&);
    int ReadAngleTEQ(Topology&, FortranData const&);
    int ReadDihedralPK(Topology&, FortranData const&);
    int ReadDihedralPN(Topology&, FortranData const&);
    int ReadDihedralPHASE(Topology&, FortranData const&);
    int ReadDihedralSCEE(Topology&, FortranData const&);
    int ReadDihedralSCNB(Topology&, FortranData const&);
    int ReadLJA(Topology&, FortranData const&);
    int ReadLJB(Topology&, FortranData const&);
    int ReadLJC(Topology&, FortranData const&);
    inline BondType GetBond();
    int ReadBondsH(Topology&, FortranData const&);
    int ReadBonds(Topology&, FortranData const&);
    inline AngleType GetAngle();
    int ReadAnglesH(Topology&, FortranData const&);
    int ReadAngles(Topology&, FortranData const&);
    inline DihedralType GetDihedral();
    int ReadDihedralsH(Topology&, FortranData const&);
    int ReadDihedrals(Topology&, FortranData const&);
    int ReadAsol(Topology&, FortranData const&);
    int ReadBsol(Topology&, FortranData const&);
    int ReadHBcut(Topology&, FortranData const&);
    int ReadAtomTypes(Topology&, FortranData const&);
    int ReadItree(Topology&, FortranData const&);
    int ReadJoin(Topology&, FortranData const&);
    int ReadIrotat(Topology&, FortranData const&);
    int ReadBox(FortranData const&);
    int ReadCapInfo(Topology&, FortranData const&);
    int ReadCapInfo2(Topology&, FortranData const&);
    int ReadGBradiiSet(Topology&);
    int ReadGBradii(Topology&, FortranData const&);
    int ReadGBscreen(Topology&, FortranData const&);
    int ReadIpol(Topology&, FortranData const&);
    int ReadPolar(Topology&, FortranData const&);
    // Extra PDB Info
    int ReadPdbRes(Topology&, FortranData const&);
    int ReadPdbChainID(Topology&, FortranData const&);
    int ReadPdbIcode(Topology&, FortranData const&);
    int ReadPdbAlt(Topology&, FortranData const&);
    int ReadPdbBfactor(Topology&, FortranData const&);
    int ReadPdbOccupancy(Topology&, FortranData const&);
    int ReadPdbNumbers(Topology&, FortranData const&);
    // CHAMBER
    int ReadChamberFFtype(Topology&, FortranData const&);
    int ReadChamberUBCount(Topology&, FortranData const&);
    int ReadChamberUBTerms(Topology&, FortranData const&);
    int ReadChamberUBFC(Topology&, FortranData const&);
    int ReadChamberUBEQ(Topology&, FortranData const&);
    int ReadChamberNumImpropers(Topology&, FortranData const&);
    int ReadChamberNumImpTerms(Topology&, FortranData const&);
    int ReadChamberImpropers(Topology&, FortranData const&);
    int ReadChamberImpFC(Topology&, FortranData const&);
    int ReadChamberImpPHASE(Topology&, FortranData const&);
    int ReadChamberLJ14A(Topology&, FortranData const&);
    int ReadChamberLJ14B(Topology&, FortranData const&);
    int ReadCmapCounts(FlagType,FortranData const&);
    int ReadCmapRes(FlagType,Topology&, FortranData const&);
    int ReadCmapGrid(FlagType,const char*, Topology&, FortranData const&);
    int ReadCmapTerms(FlagType,Topology&, FortranData const&);
    // LES
    int ReadLESntyp(Topology&, FortranData const&);
    int ReadLESfac(Topology&, FortranData const&);
    int ReadLEStypes(Topology&, FortranData const&);
    int ReadLEScnum(Topology&, FortranData const&);
    int ReadLESid(Topology&, FortranData const&);

    // ----- Write -------------------------------
    FortranData WriteFormat(FlagType) const;
    int BufferAlloc(FlagType, FortranData const&, int, int);
    int BufferAlloc(FlagType, int, int);
    int BufferAlloc(FlagType f, int n) { return BufferAlloc(f, n, -1); }
    int WriteLJ(FlagType, FlagType, NonbondArray const&);
    int WriteBondParm(FlagType, FlagType, BondParmArray const&);
    int WriteBonds(FlagType, BondArray const&);
    int WriteAngles(FlagType, AngleArray const&);
    int WriteDihedrals(FlagType, DihedralArray const&);
    void WriteLine(FlagType, std::string const&);
    int WriteTreeChainClassification(std::vector<NameType> const&);
    int WriteIjoin(std::vector<int> const&);
    int WriteIrotat(std::vector<int> const&);
    int WriteExtra(Topology const&, int);
 
    static const int AMBERPOINTERS_;
    /// Contain topology flags enumerated by FlagType
    static const ParmFlag FLAGS_[];

    ParmType ptype_;
    BufferedFrame file_;
    double elec_to_parm_; ///< Convert elec to topology units
    double parm_to_elec_; ///< Convert topology units to elec.

    // Read variables
    Iarray values_; ///< Values read in from POINTERS
    Iarray atomicNums_; ///< Set to atomic numbers if ATOMIC_NUMBER section found.
    Box parmbox_; ///< Box coords/type, set from beta, x, y, and z.
    int numLJparm_; ///< Number of LJ parameters
    bool SCEE_set_; ///< True if SCEE section found
    bool SCNB_set_; ///< True if SCNB section found
    bool atProblemFlag_; ///< True if a problematic flag was encountered and needs to be skipped.

    // CHAMBER variables
    static const double ELECTOCHAMBER_;
    static const double CHAMBERTOELEC_;
    int UB_count_[2];   ///< Urey-Bradley count: # bonds (x3), # parameters
    int N_impropers_;   ///< Number of impropers (x5)
    int N_impTerms_;    ///< Number of improper terms
    int n_cmap_terms_;  ///< Number of CMAP terms
    int n_cmap_grids_;  ///< Number of CMAP grids

    // LES variables
    int nlestyp_; ///< Number of LES types

    // Write variables
    bool writeChamber_;     ///< If true write CHAMBER info
    bool writeEmptyArrays_; ///< If true try to write TREE, IROTATE, JOIN even if not present 
    bool writePdbInfo_;     ///< If true write chain IDs etc
};
// -----------------------------------------------------------------------------
class Parm_Amber::FortranData {
  public:
    FortranData() : fstr_(0), ftype_(UNKNOWN_FTYPE), fncols_(0), fwidth_(0), fprecision_(0) {}
    /// CONSTRUCTOR - from format string
    FortranData(const char*);
    /// CONSTRUCTOR - does not set format string 
    FortranData(Type tIn, int colsIn, int widthIn, int precIn) :
      ftype_(tIn), fncols_(colsIn), fwidth_(widthIn), fprecision_(precIn) {}
    int ParseFortranFormat(const char*);

    const char* Fstr() const { return fstr_; }
    Type Ftype()       const { return ftype_; }
    int Ncols()        const { return fncols_; }
    int Width()        const { return fwidth_; }
    int Precision()    const { return fprecision_; }
  private:
    const char* fstr_;
    Type ftype_;
    int fncols_;
    int fwidth_;
    int fprecision_;
};
#endif
