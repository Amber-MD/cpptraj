#ifndef INC_TOPINFO_H
#include "CpptrajFile.h"
#include "Topology.h"
/// Class for printing formatted topology info to a file.
class TopInfo {
  public:
    TopInfo () : outfile_(0), parm_(0), toStdout_(false) {}
    ~TopInfo();
    TopInfo(Topology*);
    int SetupTopInfo(CpptrajFile*, Topology*);
    int SetupTopInfo(Topology* p) { return SetupTopInfo(0, p); }
    int PrintAtomInfo(std::string const&) const;
    int PrintShortResInfo(std::string const&, int) const;
    int PrintResidueInfo(std::string const&) const;
    int PrintBondInfo(std::string const&, std::string const&) const;
    int PrintAngleInfo(std::string const&, std::string const&, std::string const&) const;
    int PrintDihedralInfo(std::string const&, std::string const&,
                          std::string const&, std::string const&) const;
  private:
    inline int SetupMask(CharMask&) const;
    inline int SetupMask(std::string const&, CharMask&) const;
    void PrintBonds(BondArray const&, BondParmArray const&,
                    CharMask const&, CharMask const&, int&) const;
    void PrintAngles(AngleArray const&, AngleParmArray const&,
                     CharMask const&, CharMask const&, CharMask const&,
                     int&) const;
    void PrintDihedrals(DihedralArray const&, DihedralParmArray const&,
                        CharMask const&, CharMask const&,
                        CharMask const&, CharMask const&, int&) const;

    CpptrajFile* outfile_;
    Topology* parm_;
    bool toStdout_;
};
#endif
