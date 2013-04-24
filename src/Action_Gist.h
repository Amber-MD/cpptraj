#ifndef INC_ACTION_GIST_H
#define INC_ACTION_GIST_H
#include "Action.h"
#include "Vec3.h"
#include "Matrix_3x3.h"
#include "ImagedAction.h"
//#include "Residue.h"
//#include "Topology.h"

// Class: Action_Gist
/// Calculate water energy and entropy
//class Action_Gist: public Action, ImagedAction  {
class Action_Gist: public Action  {
  public:
    Action_Gist();
    ~Action_Gist();

    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Gist(); }
    static void Help();

    void Print();
    void PrintDX(std::string const&);
  private:
    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);

    /** File containing all GIST energies for each voxel */
    std::string datafile_;
    /** Data set list for all data sets created here */
    DataSetList myDSL_;
    Topology* CurrentParm_;          ///< Set to the current topology file.
    bool watermodel_;
    bool useTIP3P_;
    bool useTIP4P_;
    bool useTIP4PEW_;
//    AtomMask Mask1_;                 ///< Calculate energy for atoms in mask
//    AtomMask Mask2_;                 ///< Calculate energy for atoms in mask
//    AtomMask Mask3_;                 ///< Calculate energy for atoms in mask
//    std::string solvname_;
//    std::vector<Residue> solvent_residues_;

    // other constants
    int NFRAME_;
    double BULK_DENS_;
    double BULK_E_;

    //Grid Stuff   
    Vec3 gridcntr_;    
    int* griddim_;
    int* gridindex_;
    //Vec3 griddim_; 
    Vec3 gridorig_; 
    double gridspacn_;
    Vec3 x_lab, y_lab, z_lab;
    Vec3 x_res, y_res, z_res;
    int resnum;
    int resnum2;
    int MAX_GRID_PT_;
    std::vector <int> gridwat_;
    double Vvox_;
    std::vector<double> grid_x_;
    std::vector<double> grid_y_;
    std::vector<double> grid_z_;
    /// Return X coordinate of bin center
    double Xcrd(int i) { return (double)i*gridspacn_ - gridcntr_[0] + 0.5*gridspacn_; }
    /// Return Y coordinate of bin center
    double Ycrd(int j) { return (double)j*gridspacn_ - gridcntr_[1] + 0.5*gridspacn_; }
    /// Return Z coordinate of bin center
    double Zcrd(int k) { return (double)k*gridspacn_ - gridcntr_[2] + 0.5*gridspacn_; }

    //non-bond energy stuff
    double ELJ_;                     ///< Total VDW energy over all selected atoms.
    double Eelec_;                   ///< Total elec. energy over all selected atoms.
    double kes_;                     ///< Electrostatic constant, 1.0 when using Amber units
    //std::vector<double> atom_evdw_;  ///< Cumulative Evdw on each atom
    //std::vector<double> atom_eelec_; ///< Cumulative Eelec on each atom
    std::vector<double> atom_charge_;
    std::vector<double> wh_evdw_;
    std::vector<double> wh_eelec_;
    std::vector<double> ww_evdw_;
    std::vector<double> ww_eelec_;
    std::vector < std::vector <double> > ww_Eij_;
    std::vector<double> dEwh_dw_;
    std::vector<double> dEww_dw_ref_;
    std::vector<double> dEwh_norm_;
    std::vector<double> dEww_norm_ref_;
//    DataSet_Matrix* ww_Eij_;	// upper left triangular matrix - not what we want
//    float ** ww_Eij_ = new float * [MAX_GRID_PT];	//lower left triangular matrix
//    for (a=1; a<MAX_GRID_PT; a++)
//	Eij[a] = new float [a];

    // entropy stuff
    std::vector<double> nwat_;
    std::vector<double> nw_angle_;
    std::vector<double> TSNN_;
    std::vector<double> TSwNN_;
    std::vector<double> TStrans_dw_;
    std::vector<double> TStrans_norm_;
    double TSNNtot_;
    double max_nwat_;
    double TStranstot_;
    std::vector < std::vector <double> > the_vox_;
    std::vector < std::vector <double> > phi_vox_;
    std::vector < std::vector <double> > psi_vox_;
    std::vector <double> g_;
    std::vector <double> dens_;

    void NonbondEnergy2(Frame *, Topology *);
    void Grid(Frame*, Topology *);
    void EulerAngle(Frame *, Topology *);
};
#endif
