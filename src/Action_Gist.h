#ifndef INC_ACTION_GIST_H
#define INC_ACTION_GIST_H
#include "Action.h"
#include "Vec3.h"
#include "Matrix_3x3.h"
#include "ImagedAction.h"

// Class: Action_Gist
/// Calculate water energy and entropy
class Action_Gist: public Action, ImagedAction  {
  public:
    Action_Gist();
    ~Action_Gist();

    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Gist(); }
    static void Help();

    void Print();
    void PrintDX(std::string const&, std::vector<float>&);
    void PrintDX_double(std::string const&, std::vector<double>&);
    void PrintOutput(std::string const&);

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
/*    bool watermodel_;
    bool useTIP3P_;
    bool useTIP4P_;
    bool useTIP4PEW_;
    bool useTIP5P_;
    bool useTIP3PFW_;
    bool useSPCE_;
    bool useSPCFW_;
*/    bool doOrder_;
    bool doEij_;

    //time 
    clock_t gist_t_begin,gist_t_end;
    double gist_t_diff;

    // other constants
    int NFRAME_;		// total number of frames analyzed
//    double BULK_E_;		// bulk water energy

    //Grid Stuff   
    Vec3 gridcntr_;    		// coordiantes of grid center
    std::vector <int> griddim_;		// grid dimension (bin number in each direction)
    Vec3 gridorig_;		// coordinates of grid origin
    double gridspacn_;		// grid spacing
    double BULK_DENS_;		// reference density from bulk water
    int MAX_GRID_PT_;
    std::vector <int> gridwat_;		// voxel index of each water
    std::vector <int> nwat_;		// total number of water found in each voxel
    std::vector <int> nH_;			// total number of hydrogen found in each voxel
    std::vector <int> nw_angle_;	// total nuber of Euler angles found in each voxel
    std::vector <float> g_;		// normalized water density
    std::vector <float> gH_;		// normalized H density
    std::vector <float> dens_;		// water density
    std::vector <float> grid_x_;	// voxel index in x
    std::vector <float> grid_y_;
    std::vector <float> grid_z_;
    std::vector <float> neighbor_;		// number of water neighbor within 3.5A
    std::vector <double> qtet_;		// tetahedral order parameter
    double Vvox_;			// voxel volume
    /// Return X coordinate of bin center
    double Xcrd(int i) { return (double)i*gridspacn_ + gridorig_[0] + 0.5*gridspacn_; }
    /// Return Y coordinate of bin center
    double Ycrd(int j) { return (double)j*gridspacn_ + gridorig_[1] + 0.5*gridspacn_; }
    /// Return Z coordinate of bin center
    double Zcrd(int k) { return (double)k*gridspacn_ + gridorig_[2] + 0.5*gridspacn_; }
    
    double Lx, Ly, Lz;		// box length
    double G_max_x, G_max_y, G_max_z;		// grid max length
    void pbc(Vec3& r) {
	if (r[0] < -Lx/2) r[0] += Lx;
        else if (r[0] > Lx/2) r[0] -= Lx;
        if (r[1] < -Ly/2) r[1] += Ly;
        else if (r[1] > Ly/2) r[1] -= Ly;
        if (r[2] < -Lz/2) r[2] += Lz;
        else if (r[2] > Lz/2) r[2] -= Lz;
    } 

    //general loop    
    Topology::mol_iterator solvmol, solvmol2;
    int voxel;
    int resnum,resnum2,resindex1,resindex2;
    double theta, phi, psi;
    
    //non-bond energy stuff
    Matrix_3x3 ucell, recip;
    /*    double x_0, x_1, x_2, x_3, x_4;
    double y_0, y_1, y_2, y_3, y_4;
    double z_0, z_1, z_2, z_3, z_4;*/
    
    /*    std::vector <double> x_;
    std::vector <double> y_;
    std::vector <double> z_;*/
    std::vector <double> wh_evdw_;
    std::vector <double> wh_eelec_;
    std::vector <double> ww_evdw_;
    std::vector <double> ww_eelec_;
    std::vector < std::vector <float> > ww_Eij_;
    std::vector <float> Ewh_dw_;
    std::vector <float> Ewh_norm_;
    std::vector <float> Eww_norm_;
    std::vector <float> Eww_dw_;
    std::vector <float> neighbor_dw_;
    std::vector <float> neighbor_norm_;
    std::vector <float> pol_;

//    DataSet_Matrix* ww_Eij_;	// upper left triangular matrix - not what we want
//    float ** ww_Eij_ = new float * [MAX_GRID_PT];	//lower left triangular matrix
//    for (a=1; a<MAX_GRID_PT; a++)
//	Eij[a] = new float [a];

    // entropy stuff
    std::vector <float> dTSorient_dw_;
    std::vector <double> dTSorient_norm_;
    std::vector <float> dTStrans_dw_;
    std::vector <float> dTStrans_norm_;
    double dTSorienttot_;
    int max_nwat_;
    double dTStranstot_;
    double Ewhtot_; 
    double Ewwtot_;
    std::vector < std::vector <float> > the_vox_;
    std::vector < std::vector <float> > phi_vox_;
    std::vector < std::vector <float> > psi_vox_;

    // dipole stuffs
    std::vector <double> dipolex_;
    std::vector <double> dipoley_;
    std::vector <double> dipolez_;

    void NonbondEnergy(Frame *);
    void Grid(Frame *);
    void EulerAngle(Frame *);
    void Dipole(Frame *);
    void Order(Frame *);
};
#endif
