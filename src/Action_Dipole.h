#ifndef INC_ACTION_DIPOLE_H
#define INC_ACTION_DIPOLE_H
#include "Action.h"
#include "Grid.h"
class Action_Dipole : public Action {
  public:
    Action_Dipole();
    void print();
  private:
    int init();
    int setup();
    int action();

    Grid grid_;
    std::vector<double> dipolex_;
    std::vector<double> dipoley_;
    std::vector<double> dipolez_;
    std::string filename_;
    AtomMask mask_;
    double max_;
};
#endif
