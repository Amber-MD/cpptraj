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

    void Print();
  private:
    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);

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
