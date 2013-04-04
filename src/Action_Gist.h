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

    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Gist(); }
    static void Help();

  private:
    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print() {}

    DataSet* gist_;  ///< Will hold DataSet of calculated gist values
    Topology* CurrentParm_;          ///< Set to the current topology file.
    bool watermodel_;
    bool useTIP3P_;
    bool useTIP4P_;
    bool useTIP4PEW_;
    AtomMask Mask1_;                 ///< Calculate energy for atoms in mask
    AtomMask Mask2_;                 ///< Calculate energy for atoms in mask
    AtomMask Mask3_;                 ///< Calculate energy for atoms in mask
    std::string solvname_;
    std::vector<Residue> solvent_residues_;

    //    Box gridbox_;    
    Vec3 gridcntr_;    
    Vec3 griddim_; 
    Vec3 gridorig_; 
    double gridspacn_;
    Vec3 x_lab, y_lab, z_lab;
    Vec3 x_res, y_res, z_res;

    //non-bond energy stuff
    double ELJ_;                     ///< Total VDW energy over all selected atoms.
    double Eelec_;                   ///< Total elec. energy over all selected atoms.
    double kes_;                     ///< Electrostatic constant, 1.0 when using Amber units
    std::vector<double> atom_evdw_;  ///< Cumulative Evdw on each atom
    std::vector<double> atom_eelec_; ///< Cumulative Eelec on each atom
    std::vector<double> atom_charge_;


    //grid stuff
    int resnum;
    double MAX_GRID_PT;
    std::vector<double> gridwat_;
    std::vector < std::vector <double> > the_vox_;
    std::vector < std::vector <double> > phi_vox_;
    std::vector < std::vector <double> > psi_vox_;

    void NonbondEnergy2(Frame *, Topology *, AtomMask &, AtomMask &);
    //    void Grid(Frame*, Topology *, AtomMask &);
    void Grid(Frame*, Topology *);
    void EulerAngle(Frame *, Topology *);

};
#endif
