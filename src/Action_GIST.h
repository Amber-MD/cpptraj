#ifndef INC_ACTION_GIST_H
#define INC_ACTION_GIST_H
#include "Action.h"
// Class: Action_GIST
/// Calculate water energy and entropy
//class Action_GIST: public Action, ImagedAction  {
class Action_GIST: public Action  {
  public:
    Action_GIST();

    static DispatchObject* Alloc() { return (DispatchObject*)new Action_GIST(); }
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

    //    Box gridbox_;    
    Vec3 gridcntr_;    
    Vec3 griddim_; 
    double gridspacn_;

    //non-bond energy stuff
    double ELJ_;                     ///< Total VDW energy over all selected atoms.
    double Eelec_;                   ///< Total elec. energy over all selected atoms.
    double kes_;                     ///< Electrostatic constant, 1.0 when using Amber units
    std::vector<double> atom_evdw_;  ///< Cumulative Evdw on each atom
    std::vector<double> atom_eelec_; ///< Cumulative Eelec on each atom
    std::vector<double> atom_charge_;

    void NonbondEnergy2(Frame *, Topology *, AtomMask &, AtomMask &);
    
};
#endif
