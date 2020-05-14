#ifndef INC_TOPINFO_H
#define INC_TOPINFO_H
#include "ParameterTypes.h"
#include "Frame.h"
class CpptrajFile;
class DataSet_Coords;
class Topology;
class CharMask;
/// Class for printing formatted topology info to a file.
class TopInfo {
  public:
    /// CONSTRUCTOR
    TopInfo();
    /// CONSTRUCTOR - Take pointer to topology, output to STDOUT
    TopInfo(Topology const*);
    ~TopInfo();
    void SetNoIntraRes(bool b) { noIntraRes_ = b; }
    int SetupTopInfo(CpptrajFile*, Topology const*, DataSet_Coords*);
    int SetupTopInfo(Topology const* p, DataSet_Coords* c) { return SetupTopInfo(0, p, c); }

    int PrintAtomInfo(std::string const&) const;
    int PrintShortResInfo(std::string const&, int) const;
    int PrintResidueInfo(std::string const&) const;
    int PrintMoleculeInfo(std::string const&) const;
    int PrintShortMolInfo(std::string const&) const;
    int PrintBondInfo(std::string const&, std::string const&, bool) const;
    int PrintAngleInfo(std::string const&, std::string const&, std::string const&) const;
    int PrintDihedralInfo(std::string const&, std::string const&,
                          std::string const&, std::string const&,bool) const;
    int PrintChargeInfo(std::string const&, double&) const;
    int PrintMassInfo(std::string const&, double&) const;
  private:
    inline int SetupMask(CharMask&) const;
    inline int SetupMask(std::string const&, CharMask&) const;
    /// \return the longest atom/type name selected by mask.
    int maxAtomNamesWidth(AtomMask const&) const;
    /// \return the longest residue name in given list of residue #s
    int maxResNameWidth(std::vector<int> const&) const;
    /// \return the longest molecule name in the given list of molecle #s
    int maxMolNameWidth(std::vector<int> const&) const;
    void PrintBonds(BondArray const&, BondParmArray const&,
                    CharMask const&, CharMask const&, int, int&) const;
    void PrintAngles(AngleArray const&, AngleParmArray const&,
                     CharMask const&, CharMask const&, CharMask const&,
                     int, int&) const;
    void PrintDihedrals(DihedralArray const&, DihedralParmArray const&,
                        CharMask const&, CharMask const&,
                        CharMask const&, CharMask const&, int, int&) const;

    CpptrajFile* outfile_;
    Topology const* parm_;
    Frame coords_;
    int Awidth_;        ///< Max width of field for holding atom numbers.
    int amn_width_;     ///< Max width of AtomMaskName (:<resname>@<atomname>) for topology
    int max_type_len_;  ///< Max width of atom type name in topology
    int max_aname_len_; ///< Max width of atom name in topology
    bool toStdout_;
    bool noIntraRes_;   ///< If true, ignore intra-residue bonds/angles/dihedrals etc
};
#endif
