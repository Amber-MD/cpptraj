#ifndef INC_TOPINFO_H
#include "CpptrajFile.h"
#include "DataSet_Coords.h"
/// Class for printing formatted topology info to a file.
class TopInfo {
  public:
    TopInfo () : outfile_(0), parm_(0), toStdout_(false) {}
    ~TopInfo();
    TopInfo(Topology*);
    int SetupTopInfo(CpptrajFile*, Topology*, DataSet_Coords*);
    int SetupTopInfo(Topology* p, DataSet_Coords* c) { return SetupTopInfo(0, p, c); }
    int PrintAtomInfo(std::string const&) const;
    int PrintShortResInfo(std::string const&, int) const;
    int PrintResidueInfo(std::string const&) const;
    int PrintMoleculeInfo(std::string const&) const;
    int PrintShortMolInfo(std::string const&) const;
    int PrintBondInfo(std::string const&, std::string const&) const;
    int PrintAngleInfo(std::string const&, std::string const&, std::string const&) const;
    int PrintDihedralInfo(std::string const&, std::string const&,
                          std::string const&, std::string const&) const;
    int PrintChargeInfo(std::string const&) const;
    int PrintMassInfo(std::string const&) const;
  private:
    inline int SetupMask(CharMask&) const;
    inline int SetupMask(std::string const&, CharMask&) const;
    void PrintBonds(BondArray const&, BondParmArray const&,
                    CharMask const&, CharMask const&, int, int&) const;
    void PrintAngles(AngleArray const&, AngleParmArray const&,
                     CharMask const&, CharMask const&, CharMask const&,
                     int, int&) const;
    void PrintDihedrals(DihedralArray const&, DihedralParmArray const&,
                        CharMask const&, CharMask const&,
                        CharMask const&, CharMask const&, int, int&) const;

    CpptrajFile* outfile_;
    Topology* parm_;
    Frame coords_;
    int Awidth_; ///< Max width of field for holding atom numbers.
    int Rwidth_; ///< Max width of AtomMaskName for topology
    int max_type_len_; ///< Max width of atom type name in topology
    bool toStdout_;
};
#endif
