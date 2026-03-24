#ifndef INC_SOLVATE_H
#define INC_SOLVATE_H
#include <string>
#include <vector>
class ArgList;
class DataSet_Coords;
class DataSetList;
class Frame;
class Topology;
class Vec3;
namespace Cpptraj {
namespace Parm {
class ParameterSet;
}
namespace Structure {
/// Used to add solvent to a frame/topology pair.
class Solvate {
  public:
    /// CONSTRUCTOR
    Solvate();
    /// Initialize
    int InitSolvate(ArgList&, bool, int);
    /// Initialize for setbox
    int InitSetbox(ArgList&, int);
    /// Set VDW bounding box
    int SetVdwBoundingBox(Topology&, Frame&, Cpptraj::Parm::ParameterSet const& set0);
    /// Solvate with box
    int SolvateBox(Topology&, Frame&, Cpptraj::Parm::ParameterSet const&, DataSetList const&);

    static const char* SetboxKeywords();
    static const char* SolvateKeywords1();
    static const char* SolvateKeywords2();

    /// \return Solvent box name
    std::string const& SolventBoxName() const { return solventBoxName_; }
    /// Print solvate info to stdout
    void PrintSolvateInfo() const;
  private:
    /// \return Solvent unit selected from given DataSetList
    DataSet_Coords* GetSolventUnit(DataSetList const&) const;
    /// Adjust box widths
    void adjustBoxWidths(double&, double&, double&, double&, double&, double&) const;
    /// Calculate density
    static void calculateDensity(Topology const&, Frame const&, Cpptraj::Parm::ParameterSet const&);
    /// \return number of boxes needed to create layer
    static unsigned int nBoxesInLayer(int);
    /// \return number of boxes needed to create cube 
    static unsigned int nBoxesInCube(int);
    /// Solvate with box and specified number of solvent
    int SolvateBoxWithExactNumber(Topology&, Frame&, Cpptraj::Parm::ParameterSet const&, DataSetList const&);
    /// Get buffer arguments
    int getBufferArg(ArgList&, double);
    /// Get radii for atoms in topology
    std::vector<double> getAtomRadii(double&, Topology const&,
                                     Cpptraj::Parm::ParameterSet const&) const;
    /// Set vdW bounding box
    int setVdwBoundingBox(double&, double&, double&, std::vector<double> const&, Frame&, bool) const;
    /// Find solute atoms within a solvent box at given center
    int findCloseSoluteAtoms(std::vector<int>&, double, int, Frame const&, Vec3 const&, double, double, double
#                            ifdef CPPTRAJ_DEBUG_SOLVATE
                             , std::vector<double> const&
#                            endif
                            ) const;
    /// Determine which solvent residues do not clash with given solute atoms
    int determineValidSolventResidues(std::vector<int>&, std::vector<int> const&,
                                      Frame const&, Topology const&, Frame const&,
                                      std::vector<double> const&, std::vector<double> const&) const;
    /// Scale up buffer sizes for oct box if needed
    void octBoxCheck(Frame const&,double,double,double,std::vector<double> const&);
    /// Ewald rotate for trun. oct system
    static void ewald_rotate(Frame&, double&);
    /// Add solvent unit boxes
    int addSolventUnits(int, int, int, double, double, double, double, double, double, double,
                        Frame&, Topology const&, Frame&, Topology&,
                        std::vector<double> const&, std::vector<double> const&) const;

    int debug_;
    unsigned int nsolvent_;
    double bufferX_;
    double bufferY_;
    double bufferZ_;
    double bufferD_;
    double closeness_;
    double clipX_;
    double clipY_;
    double clipZ_;
    bool isotropic_;
    bool clip_;
    bool center_;
    bool doTruncatedOct_;
    std::string solventBoxName_;
    static const double ATOM_DEFAULT_RADIUS_; ///< Atom default radius from LEaP
    static const double CLOSENESSMODIFIER_;   ///< Overlap closeness modifier from LEaP
};
}
}
#endif
