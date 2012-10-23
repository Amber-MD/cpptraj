#ifndef INC_ACTION_ATOMICFLUCT_H
#define INC_ACTION_ATOMICFLUCT_H
#include "Action.h"
class Action_AtomicFluct : public Action {
  public :
    Action_AtomicFluct();

    static DispatchObject* Alloc() { return (DispatchObject*)new Action_AtomicFluct(); }
    static void Help();

    void Print();
  private :
    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);

    enum outputType { BYATOM = 0, BYRES, BYMASK };

    Frame SumCoords_;
    Frame SumCoords2_;
    AtomMask Mask;
    int sets_;
    int start_;
    int stop_;
    int offset_;
    int targetSet_;
    bool bfactor_;
    std::string outfilename_;
    Topology *fluctParm_;
    outputType outtype_;
    std::string setname_;
    DataSet* dataout_;
    DataFile* outfile_;
};
#endif
