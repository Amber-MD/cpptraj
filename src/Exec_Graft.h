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

    int debug_;
    Topology* newMol0Top_; ///< Hold target topology if modified.
    Topology* newMol1Top_; ///< Hold source topology if modified.
};
#endif
