#ifndef INC_ACTION_GIBBSEOFHYDRATION_H
#define INC_ACTION_GIBBSEOFHYDRATION_H
#include "Action.h"
#include "Grid.h"
/** \author Mark J. Williamson
  * \author C++ adaptation by DRR
  *  For more on the theory, please see eq. 1 in http://dx.doi.org/10.1021/ci100462t 
  */
class Action_GridFreeEnergy : public Action {
  public:
    Action_GridFreeEnergy();

    static DispatchObject* Alloc() { return (DispatchObject*)new Action_GridFreeEnergy(); }
    static void Help();


    // Action member
    void Print();
  private:
    // Action members
    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);

    /// maximum expected voxel occupancy count
    // TODO Work out a smart way of calculating this since it will be a 
    //      function of the trajectory.Maybe an upper limit of this is:
    //          numberOfVoxels * numberOfFrames?
    int maxVoxelOccupancyCount_;
    // Temperature to calculate gfe at
    double tempInKevin_;
    /// Output filename
    std::string filename_;
    /// Atom mask
    AtomMask mask_;
    Grid grid_;

};
#endif
