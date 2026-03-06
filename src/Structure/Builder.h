#ifndef INC_STRUCTURE_BUILDER_H
#define INC_STRUCTURE_BUILDER_H
#include <vector>
#include "../AtomType.h"
#ifdef TIMER
# include "../Timer.h" // DEBUG
#endif
class Atom;
class Topology;
class Frame;
class Residue;
class Vec3;
namespace Cpptraj {
namespace Parm {
class ParameterSet;
}
namespace Structure {
class Zmatrix;
class InternalCoords;
/// Used to attach different topology/frame combos using internal coordinates
class Builder {
  public:
    typedef std::vector<bool> Barray;
    /// Build type
    enum BuildType {
      BUILD = 0, ///< For cases where some external coordinates may already be known.
      SEQUENCE   ///< For cases where no external coordinates are known.
    };
    /// CONSTRUCTOR
    Builder();
    /// Set debug level
    void SetDebug(int d) { debug_ = d; }
    /// Set optional parameter set
    void SetParameters(Cpptraj::Parm::ParameterSet const*);
    /// Set atom chirality
    void SetAtomChirality(int, double);
    /// Update the internal coordinates in given Zmatrix with values from Frame/Parameters TODO combine into BuildFromInternals?
    int UpdateICsFromFrame(Frame const&, Topology const&, Barray const&);
    /// Generate internal coordinates in the same manner as LEaP
    int GenerateInternals(Frame const&, Topology const&, Barray const&);
    /// Generate internal coordinates around a link between residues in same manner as LEaP
    int GenerateInternalsAroundLink(int, int, Frame const&, Topology const&, Barray const&, BuildType, Barray const&);
    int GenerateInternalsAroundLink(int, int, Frame const&, Topology const&, Barray const&, BuildType);
    /// Update existing indices with given offset
    void UpdateIndicesWithOffset(int);

    /// \return LEaP chirality value around given atom
    static double DetermineChiralityAroundAtom(int, Frame const&, Topology const&);
    /// Build position from internals for any atom with unset position; some positions may be known.
    int BuildFromInternals(Frame&, Topology const&, Barray&) const;
    /// Build position from internals for any atom with unset position
    int BuildSequenceFromInternals(Frame&, Topology const&, Barray&, int, int) const;
    /// \return Zmatrix from current internals
    //int GetZmatrixFromInternals(Zmatrix&, Topology const&) const;
    /// Adjust torsion around a bond so that atoms with longest 'depth' are trans
    int AdjustIcAroundLink(int, int, Frame const&, Topology const&);

    static void PrintTiming(int, double);
  private:
    typedef std::vector<int> Iarray;

    /// Used to hold parameters for modeling a torsion
    class TorsionModel;
    /// Hold torsion
    class InternalTorsion;
    /// Hold angle
    class InternalAngle;
    /// Hold bond
    class InternalBond;
    /// Hold chirality
    class InternalChirality;

    class AtomIC;

    typedef std::vector<InternalTorsion> Tarray;
    typedef std::vector<InternalAngle> Aarray;
    typedef std::vector<InternalBond> Larray;
    typedef std::vector<InternalChirality> Carray;

    /// Get length parameter for atoms
    int getLengthParam(double&, int, int, Topology const&) const;
    /// Get angle parameter for atoms.
    int getAngleParam(double&, int, int, int, Topology const&) const;
    /// Calculte phi value for a torsion in a TorsionModel
    void ModelTorsion(TorsionModel const&, unsigned int, unsigned int, double);
    /// Get angle parameter/model an angle value
    double ModelBondAngle(int, int, int, Topology const&) const;
    /// Get bond parameter/model a bond length
    double ModelBondLength(int, int, Topology const&) const;
    /// Create torsion around SP3-SP3 linkage
    void createSp3Sp3Torsions(TorsionModel const&);
    /// Create torsion around SP3-SP2 linkage
    void createSp3Sp2Torsions(TorsionModel const&);
    /// Create torsion around SP2-SP2 linkage
    void createSp2Sp2Torsions(TorsionModel const&);
    /// Dtermine atom hybridization in the same way as leap
    AtomType::HybridizationType getAtomHybridization(Atom const&, std::vector<Atom> const&) const;
    /// Model torsions around a bond in the same manner as LEaP
    int assignTorsionsAroundBond(int, int, Frame const&, Topology const&, Barray const&, int);
    /// Get any existing internal torsion indices around specified atoms
    Iarray getExistingTorsionIdxs(int, int) const;
    /// Get specific internal torsion
    int getExistingTorsionIdx(int, int, int, int) const;
    /// Get specific internal angle
    int getExistingAngleIdx(int, int, int) const;
    /// Get specific internal bond
    int getExistingBondIdx(int, int) const;
    /// Get specific chirality
    int getExistingChiralityIdx(int) const;
    /// Build angle internal
    void buildAngleInternal(int, int, int, Frame const&, Topology const&, Barray const&);
    /// Build bond internal
    void buildBondInternal(int, int, Frame const&, Topology const&, Barray const&);

    /// Determine chirality around an atom
    static int determineChirality(double&, int, Frame const&, Topology const&, Barray const&);

    /// Generate internal coords for a given atom
    int generateAtomInternals(int, Frame const&, Topology const&, Barray const&);
    /// Get complete internal coords that can be used to construct specified atom
    int getIcFromInternals(InternalCoords&, int, Barray const&) const;
    /// Get two angles with ai-aj in common and all positions but ai known
    int getTwoAnglesFromInternals(InternalAngle&, InternalAngle&, InternalBond&, int, Barray const&) const;
    /// Get angle with all positions but ai known
    int getAngleFromInternals(InternalAngle&, InternalBond&, int, Barray const&) const;
    /// Get bond with aj position known
    int getBondFromInternals(InternalBond&, int, Barray const&) const;
    /// For debug, print all valid internals associated with an atom
    void printAllInternalsForAtom(int, Topology const&, Barray const&) const;
    /// \\return index of atom with longest 'depth' bonded to a given atom (ignoring one bonded atom).
    int get_depths_around_atom(int, int, Topology const&) const;
    /// Get any complete internal coords for specified atom
    AtomIC getInternalCoordsForAtom(int, int, Barray const&, Topology const&) const;
    /// Build XYZ coords for an atom from 2 angles and 1 bond
    int buildCoordsFromTwoAngles(Vec3&, int, InternalAngle const&, InternalAngle const&, InternalBond const&,
                                 Frame const&, Topology const&, Barray const&) const;
    /// Build XYZ coords for an atom if internals are available
    int buildCoordsForAtom(int, Frame&, Topology const&, Barray const&) const;
    /// Build atoms based on build priority; for cases where known atoms are sparse
    int sparseBuildFromInternals(std::vector<AtomIC>&, Frame&, Topology const&, Barray&) const;
    /// \return array containing all residues with atoms that need positions
    std::vector<Residue> residuesThatNeedPositions(Topology const&, Barray const&) const;

    int debug_;
    Cpptraj::Parm::ParameterSet const* params_;

    Topology const* currentTop_; ///< Topology for the createSpXSpXTorsions/ModelTorsion routines
    Tarray internalTorsions_;
    Aarray internalAngles_;
    Larray internalBonds_;
    Carray internalChirality_;
    Iarray Rnums_;               ///< Hold residue indices pertaining to current internals
#   ifdef TIMER
    static Timer timeg_builder_IALsetup_;
    static Timer timeg_builder_IALspan_;
    static Timer timeg_builder_IALgen_;
    static Timer timeg_builder_IALgen_internals_;
    static Timer timeg_builder_IALgen_missing_;
    static Timer timeg_builder_torsions_;
    static Timer timeg_builder_angles_;
    static Timer timeg_builder_bonds_;
    static Timer timeg_builder_chirality_;
    static Timer timeg_ang_search_;
    static Timer timeg_ang_model_;
    static Timer timeg_MBA_getAngleParam_;
    static Timer timeg_MBA_getAtomHybrid_;
#   endif
};
/// ----- Hold torsion internal ------------------
class Cpptraj::Structure::Builder::InternalTorsion {
  public:
    /// CONSTRUCTOR
    InternalTorsion() : ai_(-1), aj_(-1), ak_(-1), al_(-1), phi_(0) {}
    /// CONSTRUCTOR
    InternalTorsion(int i, int j, int k, int l, double p) :
      ai_(i), aj_(j), ak_(k), al_(l), phi_(p) {}
    /// Set the phi value in radians
    void SetPhiVal(double p) { phi_ = p; }
    /// Offset indices by given value
    void OffsetIndices(int o) {
      ai_ += o;
      aj_ += o;
      ak_ += o;
      al_ += o;
    }

    int AtI() const { return ai_; }
    int AtJ() const { return aj_; }
    int AtK() const { return ak_; }
    int AtL() const { return al_; }
    double PhiVal() const { return phi_; }
  private:
    int ai_;
    int aj_;
    int ak_;
    int al_;
    double phi_;
};
// ----- Hold angle internal ---------------------
class Cpptraj::Structure::Builder::InternalAngle {
  public:
    /// CONSTRUCTOR
    InternalAngle() : ai_(-1), aj_(-1), ak_(-1), theta_(0) {}
    /// CONSTRUCTOR
    InternalAngle(int i, int j, int k, double t) :
      ai_(i), aj_(j), ak_(k), theta_(t) {}
    /// Set the phi value in radians
    void SetThetaVal(double t) { theta_ = t; }
    /// Offset indices by given value
    void OffsetIndices(int o) {
      ai_ += o;
      aj_ += o;
      ak_ += o;
    }

    int AtI() const { return ai_; }
    int AtJ() const { return aj_; }
    int AtK() const { return ak_; }
    double ThetaVal() const { return theta_; }
  private:
    int ai_;
    int aj_;
    int ak_;
    double theta_;
};
// ----- Hold bond internal ----------------------
class Cpptraj::Structure::Builder::InternalBond {
  public:
    /// CONSTRUCTOR
    InternalBond() : ai_(-1), aj_(-1), dist_(0) {}
    /// CONSTRUCTOR
    InternalBond(int i, int j, double d) :
      ai_(i), aj_(j), dist_(d) {}
    /// Set the distance value in angstroms
    void SetDistVal(double d) { dist_ = d; }
    /// Offset indices by given value
    void OffsetIndices(int o) {
      ai_ += o;
      aj_ += o;
    }

    int AtI() const { return ai_; }
    int AtJ() const { return aj_; }
    double DistVal() const { return dist_; }
  private:
    int ai_;
    int aj_;
    double dist_;
};
// ----- Hold chirality value --------------------
class Cpptraj::Structure::Builder::InternalChirality {
  public:
    /// CONSTRUCTOR
    InternalChirality() : ai_(-1), dChi_(0) {}
    /// CONSTRUCTOR
    InternalChirality(int i, double d) : ai_(i), dChi_(d) {}
    /// Set the chirality value TODO use enum?
    void SetChiralVal(double d) { dChi_ = d; }
    /// Offset indices by given value
    void OffsetIndices(int o) { ai_ += o; }

    /// \return true if given chirality matches this one
    bool ChiralityMatches(double d) const {
      int thisC = 0;
      if (dChi_ > 0)
        thisC = 1;
      else if (dChi_ < 0)
        thisC = -1;
      int otherC = 0;
      if (d > 0)
        otherC = 1;
      else if (d < 0)
        otherC = -1;
      return (thisC == otherC);
    }

    int AtI() const { return ai_; }
    double ChiralVal() const { return dChi_; }
  private:
    int ai_;
    double dChi_;
};
// -----------------------------------------------
}
}
#endif
