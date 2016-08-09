#ifndef INC_ACTION_GIST2_H
#define INC_ACTION_GIST2_H
// TODO change protect eventually
#include "Action.h"
#include "ImagedAction.h"
#include "DataSet_3D.h"
#include "DataSet_MatrixFlt.h"
#include "Timer.h"
/// Class for applying Grid Inhomogenous Solvation Theory
class Action_GIST : public Action {
  public:
    Action_GIST();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_GIST(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print();

    inline void TransEntropy(float,float,float,float,float,float,float,int,double&,double&) const;
    static inline double Dist2(ImagingType, const double*, const double*, Box const&,
                               Matrix_3x3 const&, Matrix_3x3 const&);
    static inline void Ecalc(double, double, double, NonbondType const&, double&, double&);
    void NonbondEnergy(Frame const&, Topology const&);
    void Order(Frame const&);

    static const Vec3 x_lab_;
    static const Vec3 y_lab_;
    static const Vec3 z_lab_;
    static const double QFAC_;

    ImagedAction image_; ///< Imaging routines.
    // NOTE: '*' = Updated in DoAction(). '+' = Updated in Setup().
    // GIST float grid datasets
    DataSet_3D* gO_;        ///< Solvent oxygen density
    DataSet_3D* gH_;        ///< Solvent hydrogen density
    DataSet_3D* Esw_;       ///< Solute-water energy
    DataSet_3D* Eww_;       ///< Water-water energy
    DataSet_3D* dTStrans_;  ///< Solvent translation entropy
    DataSet_3D* dTSorient_; ///< Solvent orentational entropy
    DataSet_3D* dTSsix_;
    DataSet_3D* neighbor_norm_;
    DataSet_3D* dipole_; // pol
    // GIST double grid datasets
    DataSet_3D* order_norm_; // qtet
    DataSet_3D* dipolex_;    ///< Water dipole (X)*
    DataSet_3D* dipoley_;    ///< Water dipole (Y)*
    DataSet_3D* dipolez_;    ///< Water dipole (Z)*
    // GIST matrix datasets
    DataSet_MatrixFlt* ww_Eij_; ///< Water-water interaction energy matrix.*

    typedef std::vector<int> Iarray;
    //Iarray mol_nums_;    ///< Absolute molecule number of each solvent molecule.+ //TODO needed?
    Iarray O_idxs_;      ///< Oxygen atom indices for each solvent molecule.+
    Iarray water_voxel_; ///< Absolute grid voxel for each solvent molecule.* TODO long int?
    Iarray U_idxs_;      ///< Solute atom indices for each solute atom.+
    Iarray N_waters_;    ///< Number of waters (oxygen atoms) in each voxel.*
    Iarray N_hydrogens_; ///< Number of hydrogen atoms in each voxel.*

    typedef std::vector<float> Farray;
    Farray neighbor_; ///< Number of water neighbors within 3.5 Ang.*

    typedef std::vector<Farray> Xarray;
    Xarray voxel_xyz_; ///< Coords for all waters in each voxel.*
    Xarray voxel_Q_;   ///< w4, x4, y4, z4 for all waters in each voxel.*

    typedef std::vector<double> Darray;
    Darray E_UV_VDW_;  ///< Solute-solvent van der Waals energy for each voxel.*
    Darray E_VV_VDW_;  ///< Solvent-solvent van der Waals energy for each voxel.*
    Darray E_UV_Elec_; ///< Solute-solvent electrostatic energy for each voxel.*
    Darray E_VV_Elec_; ///< Solvent-solvent electrostatic energy for each voxel.*

    Vec3 G_max_; ///< Grid max + 1.5 Ang.

    // Timing data
    Timer gist_init_;
    Timer gist_setup_;
    Timer gist_action_;
    Timer gist_print_;
    Timer gist_grid_;
    Timer gist_nonbond_;
    Timer gist_nonbond_UV_;
    Timer gist_nonbond_VV_;
    Timer gist_euler_;
    Timer gist_dipole_;
    Timer gist_order_;

    Topology* CurrentParm_;    ///< Current topology, for energy calc.
    CpptrajFile* datafile_;    ///< GIST output
    std::string prefix_;       ///< Output file name prefix
    double BULK_DENS_;         ///< Bulk water density
    double temperature_;       ///< Temperature
    double q_O_;               ///< Charge on water oxygen
    double q_H1_;              ///< Charge on water H1
    double q_H2_;              ///< Charge on water H2 (sanity check)
    double NeighborCut2_;      ///< Cutoff for determining water neighbors (squared).
    unsigned int MAX_GRID_PT_; ///< Max number of grid points (voxels).
    unsigned int NSOLVENT_;    ///< Number of solvent molecules.
    int NFRAME_;               ///< Total # frames analyzed
    int max_nwat_;             ///< Max number of waters in any voxel
    bool doOrder_;             ///< If true do the order calc
    bool doEij_;               ///< If true do the i-j energy calc
    bool skipE_;               ///< If true skip the nonbond energy calc
};
#endif
