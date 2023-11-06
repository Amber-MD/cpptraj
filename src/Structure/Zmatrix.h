#ifndef INC_STRUCTURE_ZMATRIX_H
#define INC_STRUCTURE_ZMATRIX_H
#include "InternalCoords.h"
#include "StructureEnum.h"
#include "../Vec3.h"
class Frame;
class Topology;
class Molecule;
namespace Cpptraj {
namespace Structure {
class BuildAtom;
/// Hold internal coordinates for a system.
class Zmatrix {
    typedef std::vector<InternalCoords> ICarray;
  public:
    typedef std::vector<bool> Barray;
    /// CONSTRUCTOR
    Zmatrix();
    /// Set debug level
    void SetDebug(int d) { debug_ = d; }
    /// Print to stdout
    void print(Topology*) const;
    /// Print to stdout, no atom names
    void print() const { print(0); }

    /// \reserve space for # of internal coords TODO zero out seeds?
    void reserve(unsigned int n) { IC_.reserve( n ); }
    /// Clear the Zmatrix
    void clear();
    /// Add internal coordinate
    int AddIC(InternalCoords const&);
    /// Set specified IC
    void SetIC(unsigned int, InternalCoords const&);

    /// Set seed atoms from frame/top
    int SetSeedPositions(Frame const&, Topology const&, int, int, int);

    /// Convert specifed molecule of Frame/Topology to internal coordinates array
    int SetFromFrame(Frame const&, Topology const&, int);
    /// Convert molecule 0 of Frame/Topology to internal coordinates array
    int SetFromFrame(Frame const&, Topology const&);
    /// Get internal coordinates around bond in one direction.
    int SetupICsAroundBond(int, int, Frame const&, Topology const&,
                           std::vector<bool> const&, BuildAtom const&, BuildAtom const&);

    /// Set Frame from internal coords
    int SetToFrame(Frame&) const;
    /// Set Frame from internal coords with some positions already set
    int SetToFrame(Frame&, Barray&) const;

    typedef ICarray::const_iterator const_iterator;
    const_iterator begin() const { return IC_.begin(); }
    const_iterator end()   const { return IC_.end(); }

    /// \return Specified internal coord
    InternalCoords const& operator[](int idx) const { return IC_[idx]; }
    /// \return Array containing indices of ICs involving specified atom I
    std::vector<int> AtomI_indices(int) const;
    /// \return number of internal coords
    unsigned int N_IC() const { return IC_.size(); }
    /// \return memory usage in bytes
    unsigned int sizeInBytes() const { return (7*sizeof(int)) +
                                              (9*sizeof(double)) + // 3 Vec3
                                              (IC_.size() * InternalCoords::sizeInBytes()); }

    /// \return XYZ position of atom I for given internal coordinate
    static Vec3 AtomIposition(InternalCoords const&, Frame const&);
  private:
    typedef std::vector<int> Iarray;
    /// Simple version of auto set seeds based on connectivity only
    int autoSetSeeds_simple(Frame const&, Topology const&, Molecule const&);
    /// Calculate and add an internal coordinate given indices and Cartesian coords.
    void addIc(int,int,int,int,const double*,const double*,const double*,const double*);
    /// Add internal coordiantes by tracing a molecule
    int traceMol(int, int, int, Frame const&, Topology const&, unsigned int, unsigned int&, Barray&);
    /// Convert from Cartesian to minimal Zmatrix by tracing a molecule
    int SetFromFrame_Trace(Frame const&, Topology const&, int);
    /// \return True if IC seeds are set
    //bool HasICSeeds() const;
    /// \return True if Cartesian seeds are set
    bool HasCartSeeds() const;

    int debug_;     ///< Print debug info
    ICarray IC_;    ///< Hold internal coordinates for all atoms
    int icseed0_;   ///< Index into IC_ of first seed
    int icseed1_;   ///< Index into IC_ of second seed
    int icseed2_;   ///< Index into IC_ of third seed
    Vec3 seed0Pos_; ///< Seed 0 xyz
    Vec3 seed1Pos_; ///< Seed 1 xyz
    Vec3 seed2Pos_; ///< Seed 2 xyz
    int seedAt0_;   ///< Seed 0 topology index if seed0Pos_ is set
    int seedAt1_;   ///< Seed 1 topology index if seed1Pos_ is set
    int seedAt2_;   ///< Seed 2 topology index if seed2Pos_ is set
};
}
}
#endif
