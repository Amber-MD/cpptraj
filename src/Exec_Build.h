#ifndef INC_EXEC_BUILD_H
#define INC_EXEC_BUILD_H
#include "Exec.h"
#include "Parm/GB_Params.h"
namespace Cpptraj {
namespace Structure {
class Creator;
}
}
/// Used to build a structure 
class Exec_Build : public Exec {
  public:
    Exec_Build();
    ~Exec_Build() {}
    void Help() const;
    DispatchObject* Alloc() const { return (DispatchObject*)new Exec_Build(); }
    RetType Execute(CpptrajState&, ArgList&);

    /// Stand-alone execution - only build, do not parameterize/solvate
    RetType BuildStructure(DataSet*, DataSetList&, int, std::string const&, std::string const&, bool);
    /// \return Output COORDS set
    DataSet* OutCrdPtr() const { return outCrdPtr_; }
    /// Print Timing data
    void PrintTiming() const;
  private:
    typedef std::vector<int> Iarray;
    // Keep track of which residues are connected to each other
    typedef std::vector<Iarray> ResConnectArray;
    // For holding bonded atom pairs.
    typedef std::pair<int,int> Ipair;
    typedef std::vector<Ipair> IParray;

    /// \return true if given IParray has the given Ipair
    static inline bool hasBondingPair(IParray const&, Ipair const&);
    /// \return true if given array of residue connections has target residue index.
    static inline bool resIsConnected(Iarray const&, int);
    /// Create new topology/frame using templates
    int FillAtomsWithTemplates(Topology&, Frame&, Topology const&, Frame const&,
                               Cpptraj::Structure::Creator const&,
                               std::vector<BondType> const&);
    /// Map atoms in topology to template
    static std::vector<int> MapAtomsToTemplate(Topology const&, int, DataSet_Coords*, Cpptraj::Structure::Creator const&, std::vector<NameType>&, int&);
    /// Transfer bonds from old topology to new topology
    int transfer_bonds(Topology&, Topology const&, std::vector<BondType> const&) const;

    /// Set up output COORDS (outCrdPtr_) and Topology
    int setupOutputCoords(DataSet*, std::string const&, std::string const&, DataSetList&);
    /// Build, parameterize, solvate, and check
    RetType BuildAndParmStructure(DataSet*, std::string const&, DataSetList&, int, ArgList&, Cpptraj::Parm::GB_RadiiType);
    /// Do input structure clean/prep and fill in missing atoms from templates
    int StructurePrepAndFillTemplates(ArgList&, Topology&, Frame&, Topology&, Frame&, std::string const&, Cpptraj::Structure::Creator const&); // NOTE not const bc of timers

    int debug_;
    int check_box_natom_;         ///< Systems larger than this will have box added so PL check can be used
    bool check_structure_;        ///< If true check the resulting structure
    bool keepMissingSourceAtoms_; ///< If true, attempt to keep input atoms that are missing from templates
    bool requireAllInputAtoms_;   ///< If true, require all input atoms to be found in templates.
    bool doHisDetect_;
    bool doDisulfide_;
    bool doSugar_;
    DataSet* outCrdPtr_;          ///< Hold built (output) COORDS

    Timer t_total_;
    Timer t_hisDetect_;
    Timer t_clean_;
    Timer t_get_templates_;
    Timer t_disulfide_;
    Timer t_sugar_;
    Timer t_fill_;
    Timer t_fill_template_;
    Timer t_fill_build_;
    Timer t_fill_build_internals_;
    Timer t_fill_build_build_;
    Timer t_fill_build_link_;
    Timer t_fill_build_link_bond_;
    Timer t_assign_;
    Timer t_check_;
    Timer t_check_bonds_;
    Timer t_check_overlaps_;
    Timer t_check_rings_;
    Timer t_check_setup_;
    Timer t_solvate_;
};
#endif
