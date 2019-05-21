#ifndef INC_TRAJ_PDBFILE_H
#define INC_TRAJ_PDBFILE_H
#include "TrajectoryIO.h"
#include "PDBfile.h"
// Forward declarations
class DataSet;
/// TrajectoryIO class for reading coordinates from PDB files.
class Traj_PDBfile: public TrajectoryIO {
  public:
    // NOTE: PDBWRITEMODE must remain public for pytraj.
    /** PDBWRITEMODE: Indicate how the pdb should be written.
      *  SINGLE: Writing only a single frame.
      *  MODEL: Multiple frames written to the same file separated with 
      *         the MODEL keyword
      *  MULTI: Each frame written to a different file with name filename.frame
      */
    enum PDBWRITEMODE {NONE = 0, SINGLE, MODEL, MULTI};

    Traj_PDBfile();
    static BaseIOtype* Alloc() { return (BaseIOtype*)new Traj_PDBfile(); }
    static void WriteHelp();
  private:
    typedef std::vector<int> Iarray;
    typedef std::vector<double> Darray;
    // Inherited functions
    bool ID_TrajFormat(CpptrajFile&);
    int setupTrajin(FileName const&, Topology*);
    int setupTrajout(FileName const&, Topology*, CoordinateInfo const&,int, bool);
    int openTrajin();
    void closeTraj();
    int readFrame(int,Frame&);
    int writeFrame(int,Frame const&);
    void Info();
    int processWriteArgs(ArgList&, DataSetList const&);
    int readVelocity(int, Frame&) { return 1; }
    int readForce(int, Frame&)    { return 1; }
    int processReadArgs(ArgList&) { return 0; }
#   ifdef MPI
    // Parallel functions
    int parallelOpenTrajout(Parallel::Comm const&);
    int parallelSetupTrajout(FileName const&, Topology*, CoordinateInfo const&,
                             int, bool, Parallel::Comm const&);
    int parallelWriteFrame(int, Frame const&);
    void parallelCloseTraj() {}
#   endif

    void WriteDisulfides(Frame const&);
    void WriteBonds();
    /// Used to set up B-factor/occupancy data from DataSets
    int AssignData(Darray&, DataSet*, Topology const&, bool, const char*) const;
    /// Used to scale Bfactor/occupancy data between set values
    void ScaleData(Darray&, double, double) const;

    typedef PDBfile::SSBOND SSBOND;
    enum TER_Mode { BY_MOL = 0, BY_RES, ORIGINAL_PDB, NO_TER };
    enum Radii_Mode { GB = 0, PARSE, VDW };
    enum CONECT_Mode { NO_CONECT = 0, HETATM_ONLY, ALL_BONDS };
    Radii_Mode radiiMode_;   ///< Radii to use if PQR.
    TER_Mode terMode_;       ///< TER card mode.
    CONECT_Mode conectMode_; ///< CONECT record mode.
    PDBWRITEMODE pdbWriteMode_;
    int pdbAtom_;
    int currentSet_;
    int ter_num_;       ///< Amount to increment atom number for TER
    bool dumpq_;        ///< If true print charges/radii in Occupancy column (PQR).
    bool pdbres_;       ///< If true convert Amber res names to PDBV3 style.
    bool pdbatom_;      ///< If true convert Amber atom names to PDBV3 style.
    bool write_cryst1_; ///< If false write CRYST1 for first frame only.
    bool include_ep_;   ///< If true include extra points.
    bool prependExt_;
    bool firstframe_;   ///< Set to false after first call to writeFrame
    bool bfacscale_;               ///< If specified scale B-factor data
    bool occscale_;                ///< If specified scale occupancy data
    bool bfacbyres_;               ///< If true do bfactor data by residue
    bool occbyres_;                ///< If true do occupancy data by residue
    std::string space_group_;
    Darray Bfactors_;              ///< Hold data for B-factor column.
    Darray Occupancy_;             ///< Hold data for occupancy column.
    Iarray TER_idxs_;              ///< TER card indices.
    Iarray atrec_;                 ///< Hold ATOM record #s for CONECT
    std::vector<bool> resIsHet_;   ///< True if residue needs HETATM records
    std::vector<SSBOND> ss_residues_;
    Iarray ss_atoms_;
    Topology *pdbTop_;
    PDBfile file_;
    std::vector<char> chainID_;      ///< Hold chainID for each residue.
    std::vector<NameType> resNames_; ///< Hold residue names.
    char chainchar_;
    DataSet* bfacdata_;
    DataSet* occdata_;
    double bfacmax_;
    double occmax_;
};
#endif
