#ifndef INC_ACTION_GIST_H
#define INC_ACTION_GIST_H
#include "Action.h"
#include "GridMover.h"
#include "ImageOption.h"
#include "Timer.h"
#include "EwaldOptions.h"
#include "CharMask.h"
#include "GridBin.h"
#include "PairList.h"
#include <map>
#ifdef CUDA
#include "cuda_kernels/GistCudaSetup.cuh"
#endif
#ifdef LIBPME
#include "GIST_PME.h"
#endif
class DataSet_3D;
class DataSet_MatrixFlt;
class DataSet_GridFlt;
class DataSet_GridDbl;

/// Class for applying Grid Inhomogenous Solvation Theory
/** \author Daniel R. Roe
  */
class Action_GIST : public Action {
  public:
    Action_GIST();
    #ifdef CUDA
    ~Action_GIST() {delete[] this->solvent_;}
    #endif
    //~Action_GIST(); // DEBUG MPI
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_GIST(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print();
#   ifdef MPI
    int SyncAction();
    int ParallelPostCalc();
#   endif

    typedef std::vector<float> Farray;
    typedef std::vector<int> Iarray;
    typedef std::vector<Farray> Xarray;
    typedef std::vector<double> Darray;

    struct SolventInfo
    {
        std::vector<int> i_element;
        std::vector<std::string> unique_elements;
        std::vector<int> element_count;
    };

    /**
     * @brief Print whitespace delimited data to a CpptrajFile
     */
    class DataFilePrinter
    {
      public:
        DataFilePrinter(CpptrajFile& df, const TextFormat& fltFmt, const TextFormat& intFmt)
          : df_(&df), fltFmt_(fltFmt), intFmt_(intFmt), is_new_line_(true) {}

        void print(double d)             { maybe_space(); df_->Printf(fltFmt_.Fmt().c_str(), d); }
        void print(int i)                { maybe_space(); df_->Printf(intFmt_.Fmt().c_str(), i); }
        void print(unsigned int i)       { print(static_cast<int>(i)); }
        void print(const char* s)        { maybe_space(); df_->Printf("%s", s); }
        void print(const std::string& s) { print(s.c_str()); }
        
        void newline() {
            df_->Printf(" \n");
            is_new_line_ = true;
        }

        template<typename T>
        DataFilePrinter& operator<<(T val) { print(val); return *this; }

      private:
        void maybe_space() {
          if (!is_new_line_) {
              df_->Printf(" ");
          }
          is_new_line_ = false;
        }

        CpptrajFile* df_;
        TextFormat fltFmt_;
        TextFormat intFmt_;
        bool is_new_line_;
    };
#   ifdef MPI
    void sync_Xarray(Xarray&) const;
#   endif
    static inline void Ecalc(double, double, double, NonbondType const&, double&, double&);
    void NonbondEnergy_pme(Frame const&);
    void NonbondEnergy(Frame const&, Topology const&);
    void Order_PL(Frame const&);
    void Order(Frame const&);
    // void SumEVV();
    void CollectEnergies();
    void CalcAvgVoxelEnergy_PME(double, DataSet_3D&, DataSet_3D&, Farray&) const;
    void CalcAvgVoxelEnergy(double, DataSet_3D&, DataSet_3D&, Farray&, Farray&,
                            DataSet_3D&, DataSet_3D&, Farray&);
    DataSet_3D* AddDatasetAndFile(const std::string& name, const std::string& filename, DataSet::DataType dtype);
    int setSolventProperties(const Molecule& mol, const Topology& top);
    int checkSolventProperties(const Molecule& mol, const Topology& top) const;
    void setSolventType(const Topology& top);          ///< Set the solventType_ for each atom based on solventNames_.
    void setSoluteSolvent(const Topology& top);        ///< Set atomIsSolute, U_idxs_, (and solvent_).
    int calcVoxelIndex(double x, double y, double z);
    void analyzeSolventElements(const Molecule& mol, const Topology& top);
    bool setRigidAtomIndices(const Molecule& mol, const Topology& top);
    bool createMoleculeDatasets();
    bool createAtomDensityDatasets();
    std::vector<DataSet_3D*> getDensityDataSets();
    Farray DataSetAsArray(const DataSet_3D& ds) const;
    double SumDataSet(const DataSet_3D& ds) const;
    double SumDataSet(const std::string& name) const;
    void ScaleDataSet(DataSet_3D& ds, double factor) const;
    void ScaleFarray(Farray& ds, double factor) const;
    Vec3 calcMolCenter(const ActionFrame& frm, int begin, int end) const;
    bool isMainSolvent(int atom) const;

    template<typename ARRAY_TYPE>
    void CopyArrayToDataSet(const ARRAY_TYPE& arr, DataSet_3D& ds) const;

    template<typename T, typename ARR>
    std::vector<T> NormalizeDataSet(const DataSet_3D& ds, const ARR& norm) const;

    template<typename T>
    std::vector<T> WeightDataSet(const DataSet_3D& ds, double factor) const;
    template<typename T>
    std::vector<T> DensityWeightDataSet(const DataSet_3D& ds) const;

    int CalcTranslationalEntropy(unsigned int, unsigned int) const;

    int debug_;      ///< Action debug level
    int numthreads_; ///< Number of OpenMP threads

#ifdef CUDA
    // Additional data for GPU calculation

    std::vector<float> lJParamsA_;
    std::vector<float> lJParamsB_;
    std::vector<float> charges_;
    std::vector<int> atomTypes_;
    std::vector<int> NBIndex_;
    std::vector<int> molecule_;

    unsigned int numberAtoms_;
    int numberAtomTypes_;
    int headAtomType_;
    bool *solvent_; // This is needed additionally to atomIsSolute_ because we need a bool pointer to copy to the GPU.

    // Arrays on GPU
    int *NBindex_c_;
    void *molecule_c_;
    void *paramsLJ_c_;
    float *max_c_;
    float *min_c_;
    float *result_eww_c_;
    float *result_esw_c_;
    int *result_O_c_;
    int *result_N_c_;

    std::vector<float> E_UV_f_;
    std::vector<float> E_VV_f_;

    // CUDA only functions
    void freeGPUMemory(void);
    void copyToGPU(void);
    void NonbondCuda(ActionFrame const&);

#endif

#   ifdef LIBPME
    GIST_PME gistPme_;     ///< Holds GIST PME functionality
#   endif
    bool usePme_;          ///< If true, try to use GIST PME
    EwaldOptions pmeOpts_; ///< Hold PME options for GIST PME

    static const Vec3 x_lab_;
    static const Vec3 y_lab_;
    static const Vec3 z_lab_;
    static const double maxD_;
    static const double QFAC_;
    static const int OFF_GRID_; ///< Value in atom_voxel_ that indicates atom is off the grid
    static const int UNKNOWN_MOLECULE_;

    double gridspacing_;
    Vec3 gridcntr_;
    int griddim_[3];
    DataSet_3D* masterGrid_; ///< Grid that will be used to determine voxels for all grids
    const GridBin* gridBin_; ///< Hold the GridBin class from masterGrid_
    GridBin borderGrid_;     ///< Hold dims for grid + 1.5 Ang buffer region.
    Matrix_3x3 borderGridUcell0_; ///< Hold border grid original unit cell (in case of rotation).

    Cpptraj::GridMover mover_; ///< Used to move the master grid if necessary
    AtomMask moveMask_;        ///< Select atoms used to move the grid if necessary

    PairList pairList_;        ///< Pair list for order calc
    bool use_PL_;              ///< If true, user wants to use pair list
    bool PL_active_;           ///< If true, pairlist can be used
    double PL_cut_;            ///< Pair list cutoff

    std::vector<std::string> rigidAtomNames_;
    int rigidAtomIndices_[3]; ///< the 3 atoms that define the orientation of a solvent molecule;

    // NOTE: '*' = Updated in DoAction(). '+' = Updated in Setup().
    ImageOption imageOpt_;  ///< Used to determine if imaging should be used.*
    // GIST float grid datasets
    DataSet_3D* gO_;        ///< Solvent oxygen density
    DataSet_3D* gH_;        ///< Solvent hydrogen density
    DataSet_3D* Esw_;       ///< Solute-water energy
    DataSet_3D* Eww_;       ///< Water-water energy
    DataSet_3D* dTStrans_;  ///< Solvent translation entropy
    DataSet_3D* dTSorient_; ///< Solvent orentational entropy
    DataSet_3D* dTSsix_;    ///< Solvent entropy estimate from combined trans/rotate
    DataSet_3D* neighbor_;  ///< Number of neighbors within 3.5 Angstrom*
    DataSet_3D* dipole_;    ///< Mean dipole moment
    // GIST double grid datasets
    DataSet_3D* order_;     ///< Average tetrahedral order parameter for solvent, qtet
    DataSet_3D* dipolex_;   ///< Water dipole (X)*
    DataSet_3D* dipoley_;   ///< Water dipole (Y)*
    DataSet_3D* dipolez_;   ///< Water dipole (Z)*
    // PME GIST double grid datasets
    DataSet_3D* PME_;       ///< The PME nonbond interaction( charge-charge + vdw) cal for water
    DataSet_3D* U_PME_;     ///< The PME nonbond energy for solute atoms

    SolventInfo solventInfo_;
    std::vector<DataSet_3D*> atomDensitySets_;
    std::vector<DataSet_3D*> molDensitySets_;
    std::vector<DataSet_3D*> molEswSets_;
    std::vector<DataSet_3D*> molEwwSets_;

    std::vector<std::string> solventNames_;
    std::vector<int> solventType_;
    DataSetList* DSL_;
    DataFileList* DFL_;
    std::string dsname_;
    // GIST matrix datasets
    DataSet_MatrixFlt* ww_Eij_; ///< Water-water interaction energy matrix.*

    std::string soluteMask_; ///< Solute mask supplied by the user using [solute]
    //Iarray mol_nums_;     ///< Absolute molecule number of each solvent molecule.+ //TODO needed?
    Iarray O_idxs_;         ///< First atom indices for each solvent molecule (where atomIsSolute is true).+
    Iarray OnGrid_idxs_;    ///< Indices for each non-solute atom where the molecule is on the grid.*
    Iarray atom_voxel_;     ///< Absolute grid voxel for each atom (OFF_GRID_ if atom not on grid).*
    std::vector<bool> atomIsSolute_; ///< True if atom is solute.+
    std::vector<bool> atomIsSolventO_; ///< True if atom is sovent O. Used to choose atoms for neighbor calc.+
    Iarray U_idxs_;         ///< Atom indices for solute atoms only.+
    Iarray U_onGrid_idxs_;  ///< Indices for each solute atom on the grid.*
    Iarray N_solvent_;       ///< Number of solvent centers in each voxel.*
    Iarray N_main_solvent_; ///< Number of main solvent (usually water) centers in each voxel.*
    Iarray N_solute_atoms_; ///< Number of solute atoms in each voxel.*
    Iarray N_hydrogens_;    ///< Number of hydrogen atoms in each voxel.*
#   ifdef _OPENMP
    std::vector<Iarray> EIJ_V1_; ///< Hold any interaction energy voxel 1 each frame.*
    std::vector<Iarray> EIJ_V2_; ///< Hold any interaction energy voxel 2 each frame.*
#   endif

    std::vector<Farray> neighborPerThread_; ///< Number of water neighbors within 3.5 Ang.*
#   ifdef _OPENMP
    std::vector<Farray> EIJ_EN_;   ///< Hold any interaction energies each frame.*
#   endif

    Xarray voxel_xyz_; ///< Coords for all waters in each voxel.*
    Xarray voxel_Q_;   ///< w4, x4, y4, z4 for all waters in each voxel.*

    Darray OnGrid_XYZ_;             ///< XYZ coordinates for on-grid waters.*
    std::vector<Darray> E_UV_;  ///< Solute-solvent van der Waals energy for each atom.*
    std::vector<Darray> E_VV_;  ///< Solvent-solvent van der Waals energy for each atom.*
    // PME energy terms
    Darray E_pme_;     ///< Total nonbond interaction energy(VDW + electrostatic) calculated by PME for water TODO grid?
    Darray U_E_pme_;   ///< Total nonbond interaction energy(VDW + Elec) calculated by PME for solute TODO grid?

    // Timing data
    Timer gist_init_;
    Timer gist_setup_;
    Timer gist_action_;
    Timer gist_print_;
    Timer gist_grid_;
    Timer gist_nonbond_;
    Timer gist_nonbond_dist_;
    Timer gist_nonbond_UV_;
    Timer gist_nonbond_VV_;
    Timer gist_nonbond_OV_;
    Timer gist_euler_;
    Timer gist_dipole_;
    Timer gist_order_;
    Timer gist_print_OE_; ///< Orientational entropy calc timer
    Timer gist_print_TE_; ///< Translational entropy calc timer
    Timer gist_print_write_; ///< GIST results file write timer

    Topology* CurrentParm_;    ///< Current topology, for energy calc.
    CpptrajFile* datafile_;    ///< GIST output
    CpptrajFile* eijfile_;     ///< Eij matrix output
    CpptrajFile* infofile_;    ///< GIST info
    std::string prefix_;       ///< Output file name prefix
    std::string ext_;
    TextFormat fltFmt_;        ///< Output file format for floating point values
    TextFormat intFmt_;        ///< Output file format for integer values.
    Darray Q_;                 ///< Solvent molecule charges (for dipole calc)
    double BULK_DENS_;         ///< Bulk water density
    double temperature_;       ///< Temperature
    double NeighborCut2_;      ///< Cutoff for determining water neighbors (squared).
//    double system_potential_energy_; ///< the emsemble average potential energy ( Eelec + Vdw ) for the frames (pme only)
//    double solute_potential_energy_; ///< the ensemble average potential energy on solute atoms (pme only)
    unsigned int MAX_GRID_PT_; ///< Max number of grid points (voxels).
    unsigned int NSOLVENT_;    ///< Number of solvent molecules.
    unsigned int N_ON_GRID_;   ///< Number of water atoms on the grid.*
    unsigned int nMolAtoms_;   ///< Number of atoms in a water molecule.+
    int NFRAME_;               ///< Total # frames analyzed
    int max_nwat_;             ///< Max number of waters in any voxel
    int nNnSearchLayers_;      ///< Number of layers of voxels to search for nearest neighbors in the entropy search.
    int n_linear_solvents_;    ///< Count how many near-linear solvents occur during the GIST calculation.*
#   ifdef DEBUG_GIST
    CpptrajFile* debugOut_; ///> DEBUG
#   endif
    bool doOrder_;             ///< If true do the order calc
    bool doEij_;               ///< If true do the i-j energy calc
    bool skipE_;               ///< If true skip the nonbond energy calc
    bool skipS_;               ///< If true does not calculate entropy
    bool exactNnVolume_;       ///< If true use the exact volume equation for the NN entropy
    bool useCom_;              ///< If true use the COM as the molecular center; If false, use the first atom according to rigidAtomIndices.
    bool setupSuccessful_;     ///< Used to skip Print() if setup failed.
#   ifdef MPI
    Parallel::Comm trajComm_;  ///< Communicator across trajectory
#   endif
    int watCountSubvol_;       ///< Hold nwts; if -1, trans. entropy calc has not been run
};
#endif
