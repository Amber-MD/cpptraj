#ifndef INC_ACTION_LIE_H
#define INC_ACTION_LIE_H
#include "Action.h"
#include "ImagedAction.h"
// Class: Action_LIE
/** Action to calculate the Linear Interaction Energy (effectively the nonbonded 
  * energies between two different masks
  */
class Action_LIE: public Action, ImagedAction {
  public:
    Action_LIE();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_LIE(); }
    void Help() const;
  private:
    DataSet *elec_;  //< EEL data set
    DataSet *vdw_;   //< VDW data set
    bool dovdw_;     //< Calculate VDW contribution
    bool doelec_;    //< Calculate EEL contribution
    AtomMask Mask1_; //< Ligand mask
    AtomMask Mask2_; //< Surroundings mask
    double cut2vdw_; //< Square of the cutoff for VDW
    double dielc_;   //< dielectric constant
    double cut2elec_;//< Square of the cutoff for EEL
    double onecut2_;//< 1 / sqrt of electrostatic cutoff

    Topology* CurrentParm_; //< Topology to get params from

    // Hold atom charges * 18.2223
    std::vector<double> atom_charge_;

    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
    void Print() {}

    // Specific routines to this action
    int SetupParms(Topology const&);
    double Calculate_LJ(Frame const&, Topology const&) const;
    double Calculate_Elec(Frame const&) const;
};
#endif
