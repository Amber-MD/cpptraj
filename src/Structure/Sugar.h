#ifndef INC_STRUCTURE_SUGAR_H
#define INC_STRUCTURE_SUGAR_H
#include <vector>
// Forward declares
class Topology;

namespace Cpptraj {
namespace Structure {
/// Hold indices for a sugar residue
class Sugar {
  public:
    typedef std::vector<int> Iarray;
    /// Sugar status. Keep synced with StatTypeStr_
    enum StatType { SETUP_OK = 0,    ///< Regular sugar, setup complete
                    MISSING_O,       ///< Could not find ring oxygen
                    MULTIPLE_O,      ///< Multiple potential ring oxygen atoms
                    MISSING_CHAIN,   ///< Could not find all chain carbons
                    MISSING_ANO_REF, ///< Missing anomeric reference atom
                    MISSING_CONFIG,  ///< Missing configurational carbon
                    MISSING_C1X      ///< Missing C1 X substituent
                  };
    /// CONSTRUCTOR - Status and residue first atom, incomplete setup
    Sugar(StatType, int);
    /// CONSTRUCTOR - Status, ring O, anomeric, ring atoms, chain atoms, incomplete setup
    Sugar(StatType, int, int, Iarray const&,Iarray const&);
    /// CONSTRUCTOR - ring O, Anomeric, Anomeric Ref, Highest Sterocenter, ring atoms, chain atoms, isMissingAtoms
    Sugar(int,int,int,int,Iarray const&,Iarray const&, bool);

    /// Base ring type
    enum RingTypeEnum { PYRANOSE = 0,  ///< Ring is 5 carbons, 1 oxygen
                        FURANOSE,      ///< Ring is 4 carbons, 1 oxygen
                        UNKNOWN_RING   ///< Some unknown ring type
                      };
    /// Keep synced with RingTypeEnum
    static const char* ringstr_[];


    inline int ResNum(Topology const&) const;
    StatType Status()          const { return stat_; }
    int RingOxygenAtom()       const { return ring_oxygen_atom_; }
    int AnomericAtom()         const { return anomeric_atom_; }
    int AnomericRefAtom()      const { return ano_ref_atom_; }
    int HighestStereocenter()  const { return highest_stereocenter_; }
    RingTypeEnum RingType()    const { return ringType_; }
    bool IsMissingAtoms()      const { return isMissingAtoms_; }
    Iarray const& RingAtoms()  const { return ring_atoms_; }
    int RingEndAtom()          const { return ring_atoms_.back(); }
    Iarray const& ChainAtoms() const { return chain_atoms_; }
//    typedef std::vector<int>::const_iterator const_iterator;
//    const_iterator ringbegin() const { return ring_atoms_.begin(); }
//    const_iterator ringend()   const { return ring_atoms_.end(); }

    bool NotSet() const { return (stat_ != SETUP_OK); }
    /// \return Number of ring atoms
    unsigned int NumRingAtoms() const;
    void PrintInfo(Topology const&) const;
    /// Remap internal indices according to given atom map.
    void RemapIndices(Iarray const&, int, int);

    void SetStatus(StatType s) { stat_ = s; }
  private:
    /// Strings corresponding to StatType
    static const char* StatTypeStr_[];

    StatType stat_;            ///< Setup status
    int ring_oxygen_atom_;     ///< Index of the ring oxygen atom
    int anomeric_atom_;        ///< Index of the anomeric C atom (ring start)
    int ano_ref_atom_;         ///< Index of the anomeric reference C atom
    int highest_stereocenter_; ///< Index of the highest stereocenter in the carbon chain
    RingTypeEnum ringType_;    ///< Will be set to ring type
    bool isMissingAtoms_;      ///< True if original PDB indicates sugar has missing atoms.
    Iarray ring_atoms_;        ///< Index of all non-oxygen ring atoms
    Iarray chain_atoms_;       ///< Index of all chain carbon atoms (from anomeric carbon).
};

} // end Structure
} // end Cpptraj
#endif
