#ifndef INC_ACTION_GIBBSEOFHYDRATION_H
#define INC_ACTION_GIBBSEOFHYDRATION_H
#include "Action.h"
#include "Grid.h"
/** \author Mark J. Williamson
  * \author C++ adaptation by DRR
  */
class Action_GibbsEofHydration : public Action {
  public:
    Action_GibbsEofHydration();

    // Action member
    void print();
  private:
    // Action members
    int init();
    int setup();
    int action();

    /// maximum expected voxel occupancy count
    // TODO Work out a smart way of calculating this since it will be a 
    //      function of the trajectory.Maybe an upper limit of this is:
    //          numberOfVoxels * numberOfFrames?
    int maxVoxelOccupancyCount_;
    /// Output filename
    std::string filename_;
    /// Atom mask
    AtomMask mask_;
    Grid grid_;
};
#endif
