#ifndef INC_EXEC_ADDMISSINGRES_H
#define INC_EXEC_ADDMISSINGRES_H
#include "Exec.h"
#include "TrajectoryFile.h"
// Forward declare
class DataSet_Coords_CRD;
/// Attempt to add missing residues in a PDB 
class Exec_AddMissingRes : public Exec {
  public:
    Exec_AddMissingRes();
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_AddMissingRes(); }
    RetType Execute(CpptrajState&, ArgList&);
  private:
    typedef std::vector<std::string> Sarray;
    // Forward declaration of Pres class
    class Pres;
    typedef std::vector<Pres> ResArray;
    /// Hold "gaps"
    typedef std::vector<ResArray> Garray;
    typedef std::vector<int> Iarray;

    int FindGaps(Garray&, CpptrajFile&, std::string const&) const;
    int Minimize(Topology const&, Frame&, CharMask const&) const;
    int WriteStructure(std::string const&, Topology*, Frame const&, TrajectoryFile::TrajFormatType) const;
    static void GenerateLinearGapCoords(int, int, Frame&);
    //static void GenerateLinearTerminalCoords(int,int,int,int, Frame&);
    //int AssignLinearCoords(Topology const&, CharMask const&, Frame&) const;
    static int CalcFvecAtIdx(Vec3&, Vec3&, int, Topology const&, Frame const&, CharMask const&);
    int CoordSearchGap(int,int,int,int,Topology const&, CharMask&, Frame&) const;
    int CoordSearchTerminal(int, int, int, Topology const&, CharMask&, Frame&) const;
    int AssignCoordsBySearch(Topology const&, Frame const&, Topology const&, Frame&,
                             CharMask const&) const;

    int AddMissingResidues(DataSet_Coords_CRD*, Topology const&, Frame const&, Garray const&) const;

    int debug_;     ///< Debug level
    int nMinSteps_; ///< Number of minimization steps.
    bool optimize_; ///< If true, try to optimize coordinates.
};
#endif
