#ifndef INC_EXEC_GRAFT_H
#define INC_EXEC_GRAFT_H
#include "Exec.h"
/// Graft part of one COORDS to another COORDS 
class Exec_Graft : public Exec {
  public:
    Exec_Graft();
    ~Exec_Graft();
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_Graft(); }
    RetType Execute(CpptrajState&, ArgList&);
  private:
    typedef std::vector<int> Iarray;
    typedef std::vector<std::string> Sarray;
    /// \return Array of connect atoms from associated data
    static Iarray getConnectAtoms(AssociatedData*);
    /// print connect atoms to stdout
    static void print_connect(const char*, Iarray const&, Topology const&);
    /// Select bond index from expression in given topology
    static int select_bond_idx(std::string const&, Topology const&); 
    /// Modify topology and frame according to mask expression
    static Topology* modify_top(Topology const&, AtomMask const&, Frame&);
    /// Get COORDS data set
    static DataSet_Coords* get_crd(ArgList&, DataSetList const&, const char*, const char*, Frame&, const char*);
    /// Graft using internal coordinates
    int graft_ic(DataSet_Coords*, Topology const&, Frame const&, Topology const&, Frame const&, Sarray const&, Sarray const&) const;
    /// Graft assuming structures have been rms fit
    int graft_rms(DataSet_Coords*, Topology const&, Frame const&, Topology const&, Frame const&, Sarray const&, Sarray const&) const;

    int get_original_orientations(Topology const&, Frame const&, Topology const&, Frame const&,
                                  AtomMask const&, AtomMask const&,
                                  Sarray const&, Sarray const&);
    int debug_;
    int verbose_;          ///< Parameter assign verbosity
    Topology* newMol0Top_; ///< Hold target topology if modified.
    Topology* newMol1Top_; ///< Hold source topology if modified.
    bool hasOrient0_;
    bool hasOrient1_;
//    double orient0_; ///< When using IC, record original orientation around bonding atom 0
//    double orient1_; ///< When using IC, record original orientation around bonding atom 1
    double chi0_;    ///< When using IC, record original chirality around bonding atom 0
    double chi1_;    ///< When using IC, record original chirality around bonding atom 1
};
#endif
