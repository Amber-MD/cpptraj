#ifndef INC_STRUCTURE_ZMATRIX_H
#define INC_STRUCTURE_ZMATRIX_H
#include "InternalCoords.h"
#include "../Vec3.h"
class Frame;
class Topology;
namespace Cpptraj {
namespace Structure {
/// Hold internal coordinates for a system.
class Zmatrix {
    typedef std::vector<InternalCoords> ICarray;
    typedef std::vector<int> Iarray;
  public:
    /// CONSTRUCTOR
    Zmatrix();
    /// Set debug level
    void SetDebug(int d) { debug_ = d; }
    /// Add internal coordinate
    void AddIC(InternalCoords const&);
    /// Add internal coordinate as next available IC seed
    int AddICseed(InternalCoords const&);
    /// Set seed atoms from frame/top
    int SetSeedPositions(Frame const&, Topology const&, int, int, int);
    /// Convert Frame/Topology to internal coordinates array
    int SetFromFrame(Frame const&, Topology const&, int);

    /// Set Frame from internal coords
    int SetToFrame(Frame&) const;
    /// Print to stdout
    void print() const;
    /// \return True if IC seeds are set
    bool HasICSeeds() const;
    /// \return True if Cartesian seeds are set
    bool HasCartSeeds() const;

    typedef ICarray::const_iterator const_iterator;
    const_iterator begin() const { return IC_.begin(); }
    const_iterator end()   const { return IC_.end(); }
  private:
    //static const int DUMMY0; ///< Used to denote a dummy atom 0
    //static const int DUMMY1; ///< Used to denote a dummy atom 1
    //static const int DUMMY2; ///< Used to denote a dummy atom 2

    void addIc(int,int,int,int,const double*,const double*,const double*,const double*);

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
