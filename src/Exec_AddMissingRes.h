#ifndef INC_EXEC_ADDMISSINGRES_H
#define INC_EXEC_ADDMISSINGRES_H
#include "Exec.h"
#include "TrajectoryFile.h"
// Forward declare
class DataSet_Coords_CRD;
/// Attempt to add missing residues in a PDB 
class Exec_AddMissingRes : public Exec {
  public:
    Exec_AddMissingRes() : Exec(GENERAL) {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_AddMissingRes(); }
    RetType Execute(CpptrajState&, ArgList&);
  private:
    typedef std::vector<std::string> Sarray;
    // Forward declaration of Gap class
    class Gap;
    typedef std::vector<Gap> Garray;

    int FindGaps(Garray&, CpptrajFile&, std::string const&) const;
    int Minimize(Topology const&, Frame&, CharMask const&) const;
    int WriteStructure(std::string const&, Topology*, Frame const&, TrajectoryFile::TrajFormatType) const;
    static void GenerateLinearCoords(int, int, Frame&);
    static void GenerateTerminalCoords(int,int,int,int, Frame&);
    int AssignLinearCoords(Topology const&, CharMask const&, Frame&) const;
    static int CalcFvecAtIdx(Vec3&, Vec3&, int, Topology const&, Frame const&, CharMask const&);
    int CoordSearch(int, int, int, Topology const&, CharMask&, Frame&) const;
    int AssignCoordsBySearch(Topology const&, Frame const&, Topology const&, Frame&,
                             Garray const&, CharMask const&) const;

    int AddMissingResidues(DataSet_Coords_CRD*, Topology const&, Frame const&, Garray const&);

    int debug_;
};

// ----- Gap class -------------------------------------------------------------
/// Record a gap in PDB
class Exec_AddMissingRes::Gap {
  public:
    Gap() {}
    /// Start gap with res name, res num, and chain
    Gap(std::string const& startNameIn, int startResIn, std::string const& startChainIn) :
      resNames_(1, startNameIn), startRes_(startResIn), stopRes_(-1), chainId_(startChainIn[0])
      {}
    /// Start gap with res num and chain
    Gap(int startResIn, std::string const& startChainIn) :
      startRes_(startResIn), stopRes_(-1), chainId_(startChainIn[0])
      {}

    std::string const& FirstName() const { return resNames_.front(); }
    std::string const& LastName()  const { return resNames_.back(); }
    int StartRes()                 const { return startRes_; }
    int StopRes()                  const { return stopRes_; }
    char Chain()                   const { return chainId_; }
    unsigned int Nres()            const { return resNames_.size(); }
    typedef Sarray::const_iterator name_iterator;
    name_iterator nameBegin()      const { return resNames_.begin(); }
    name_iterator nameEnd()        const { return resNames_.end(); }

    void SetStopRes(int s) { stopRes_ = s; }
    void AddGapRes(std::string const& r) { resNames_.push_back( r ); }
  private:
    Sarray resNames_; ///< Residue names in the Gap
    int startRes_;    ///< pdb start residue number for the gap 
    int stopRes_;     ///< pdb stop residue number for the gap
    char chainId_;    ///< chain ID of the gap
};

#endif
