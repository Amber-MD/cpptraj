#ifndef INC_ACTION_RMSAVGCORR_H
#define INC_ACTION_RMSAVGCORR_H
#include "Action.h"
#include "CoordList.h"
// Class: Action_RmsAvgCorr
/// Calculate rmsd using running avg structures 
class Action_RmsAvgCorr: public Action {
  public:
    Action_RmsAvgCorr();

    static DispatchObject* Alloc() { return (DispatchObject*)new Action_RmsAvgCorr(); }
    static void Help();


    void Print();
  private:
    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);

    std::string separateName_;
    AtomMask Mask0_;
    CoordList ReferenceCoords_;
    Topology* ReferenceParm_;
    DataSet* Ct_;
    int parmNatom_;
    int maxwindow_;
    bool useMass_;
};
#endif  
