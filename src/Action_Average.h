#ifndef INC_ACTION_AVERAGE_H
#define INC_ACTION_AVERAGE_H
#include "Action.h"
// Class: Action_Average
/// Sum up all coordinates and print the averaged coords in given format.
class Action_Average: public Action {
  public:
    Action_Average();

    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Average(); }
    static void Help();

    ~Action_Average();

    void print();
  private:
    int init();
    int setup();
    int action();

    AtomMask Mask1_;
    Frame* AvgFrame_;
    Topology* AvgParm_;
    ArgList trajArgs_;
    bool parmStripped_;
    int Natom_;
    int Nframes_;
    int start_;
    int stop_;
    int offset_;
    int targetFrame_;
    std::string avgfilename_;

};
#endif  
