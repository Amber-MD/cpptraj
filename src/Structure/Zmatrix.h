#ifndef INC_STRUCTURE_ZMATRIX_H
#define INC_STRUCTURE_ZMATRIX_H
#include "InternalCoords.h"
#include "../Vec3.h"
class Frame;
class Topology;
class Molecule;
namespace Cpptraj {
namespace Structure {
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
    void print(Topology const*) const;
    /// Print to stdout, no atom names
    void print() const { print(0); }

    /// \reserve space for # of internal coords TODO zero out seeds?
    void reserve(unsigned int n) { IC_.reserve( n ); }
    /// Clear the Zmatrix
    void clear();
    /// Sort ICs by atom index (I < L < K < J)
    void sort();
    /// \return the last added IC
    InternalCoords const& back() const { return IC_.back(); }
    /// Add internal coordinate
    int AddIC(InternalCoords const&);
    /// Calculate and add internal coordinate for specified atoms
    int AddIC(int,int,int,int,Frame const&);
    /// Set specified IC
    void SetIC(unsigned int, InternalCoords const&);
    /// Offset IC indices by given value
    void OffsetIcIndices(int);

    /// Set seed atoms from frame/top
    int SetSeedPositions(Frame const&, Topology const&, int, int, int);
    /// Set seed atoms as 3 consecutive atoms with known positions for specified residue # TODO deprecate?
    int AutoSetSeedsWithPositions(Frame const&, Topology const&, int, Barray const&);

    /// Try to generate complete ICs in same manner as LEaP
    int GenerateInternals(Frame const&, Topology const&);
    /// Try to generate complete ICs from atom connectivity
    int SetFromFrameAndConnect(Frame const&, Topology const&);//, int);
    /// Convert specifed molecule of Frame/Topology to internal coordinates array
    int SetFromFrame(Frame const&, Topology const&, int);
    /// Convert molecule 0 of Frame/Topology to internal coordinates array
    int SetFromFrame(Frame const&, Topology const&);
    /// \return Modifable reference to specified IC
    InternalCoords& ModifyIC(unsigned int idx) { return IC_[idx]; }

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
    /// \return XYZ position of atom I bonded to C using positions of atoms A-C-B, C-I distance, and 2 angles
    static Vec3 PosFromBondTwoAnglesOrientation(Vec3 const&, Vec3 const&, Vec3 const&,
                                                double, double, double, double);
    /// \return XYZ position of atom I using positions of atoms J/K and distance (ang)/angle(rad)
    static Vec3 AtomIposition(Vec3 const&, Vec3 const&, double, double);
    /// \return XYZ position of atom I using position of atom J and distance (ang)
    static Vec3 AtomIposition(Vec3 const&, double);
    /// \return XYZ position of atom I using positions of atoms J/K/L and distance(ang)/angle/torsion (deg)
    static Vec3 AtomIposition(Vec3 const&, Vec3 const&, Vec3 const&, double, double, double);
    /// \return XYZ position of atom I for given internal coordinate
    static Vec3 AtomIposition(InternalCoords const&, Frame const&);
  private:
    typedef std::vector<int> Iarray;
    /// Set seeds as 3 consecutive atoms with known positions
    int autoSetSeeds_withPositions(Frame const&, Topology const&, int, int, Barray const&);
    /// Simple version of auto set seeds based on connectivity only
    int autoSetSeeds_simple(Frame const&, Topology const&, Molecule const&);
    /// Calculate and add an internal coordinate given indices and Cartesian coords.
    //inline void addIc(int,int,int,int,const double*,const double*,const double*,const double*);
    /// Calculate and add an internal coordinate given indices and Cartesian coords.
    inline void addIc(int,int,int,int,Frame const&);
    /// Add internal coordinates by tracing a molecule
    int traceMol(int, int, int, Frame const&, Topology const&, unsigned int, unsigned int&, Barray&);
    /// Add internal coordinate for given atom
    int addInternalCoordForAtom(int, Frame const&, Topology const&);
    /// Calculate coordinates for an atom based on three angles and a bond length
    static Vec3 calculatePositionFromAngles(double, double, double, double);

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
