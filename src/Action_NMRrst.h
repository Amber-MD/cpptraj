#ifndef INC_ACTION_NMRRST_H
#define INC_ACTION_NMRRST_H
#include "Action.h"
#include "ImagedAction.h"
#include "BufferedLine.h"
// Class: Action_NMRrst
class Action_NMRrst: public Action {
  public:
    Action_NMRrst();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_NMRrst(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print() {}

    int ReadXplor( BufferedLine& );
    int ReadAmber( BufferedLine& );
    /// Type to hold NOE data.
    struct noeDataType {
      AtomMask dMask1_; ///< First mask for distance pair.
      AtomMask dMask2_; ///< Second mask for distance pair.
      double bound_;    ///< Lower bound.
      double boundh_;   ///< Upper bound.
      double rexp_;     ///< Expected distance
      DataSet* dist_;   ///< Distance DataSet.
    };
    typedef std::vector<noeDataType> noeArray;
    noeArray NOEs_;
    
    ImagedAction Image_;
    bool useMass_;
};
#endif
