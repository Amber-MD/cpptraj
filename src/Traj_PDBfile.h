#ifndef INC_TRAJ_PDBFILE_H
#define INC_TRAJ_PDBFILE_H
#include "TrajectoryIO.h"
#include "PDBfile.h"
// Class: Traj_PDBfile
/// TrajectoryIO class for reading coordinates from PDB files.
class Traj_PDBfile: public TrajectoryIO {
  public:
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
    enum TER_Mode { BY_MOL = 0, BY_RES = 1, NO_TER = 2 };
    enum Radii_Mode { GB = 0, PARSE, VDW };
    Radii_Mode radiiMode_; ///< Radii to use if PQR.
    TER_Mode terMode_;     ///< TER card mode.
    int pdbAtom_;
    int currentSet_;
    int ter_num_;      ///< Amount to increment atom number for TER
    PDBWRITEMODE pdbWriteMode_;
    bool dumpq_;   ///< If true print charges/radii in Occupancy column (PQR).
    bool pdbres_;  ///< If true convert Amber res names to PDBV3 style.
    bool pdbatom_; ///< If true convert Amber atom names to PDBV3 style.
    bool write_cryst1_; ///< If false write CRYST1 for first frame only.
    bool include_ep_;   ///< If true include extra points.
    bool writeConect_;  ///< If true write CONECT records for each bond.
    bool prependExt_;
    std::string space_group_;
    std::vector<double> radii_;  ///< Hold radii for PQR format.
    std::vector<int> TER_idxs_;  ///< TER card indices.
    std::vector<int> atrec_;     ///< Hold ATOM record #s for CONECT
    std::vector<bool> resIsHet_; ///< True if residue needs HETATM records
    Topology *pdbTop_;
    PDBfile file_;

    std::vector<char> chainID_;      ///< Hold chainID for each residue.
    std::vector<NameType> resNames_; ///< Hold residue names.
    char chainchar_;

    // Inherited functions
    bool ID_TrajFormat(CpptrajFile&);
    int setupTrajin(FileName const&, Topology*);
    int setupTrajout(FileName const&, Topology*, CoordinateInfo const&,int, bool);
    int openTrajin();
    void closeTraj();
    int readFrame(int,Frame&);
    int writeFrame(int,Frame const&);
    void Info();
    int processWriteArgs(ArgList&);
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
};
#endif
