#ifndef INC_ACTION_GRID_H
#define INC_ACTION_GRID_H
#include "Action.h"
#include "Grid.h"
class Action_Grid : public Action {
  public:
    Action_Grid();

    void print();
  private:
    double max_;
    double madura_;
    double smooth_;
    bool invert_;
    AtomMask mask_;
    std::string filename_;
    std::string pdbname_;
    Grid grid_;

    int init();
    int setup();
    int action();
};
#endif
