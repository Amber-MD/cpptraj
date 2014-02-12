#ifndef INC_ACTION_VELOCITYAUTOCORR_H
#define INC_ACTION_VELOCITYAUTOCORR_H
#include "Action.h"
class Action_VelocityAutoCorr : public Action {
  public:
    Action_VelocityAutoCorr();
    static DispatchObject* Alloc() { return (DispatchObject*)new Action_VelocityAutoCorr(); }
    static void Help();
  private:
    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void Print() {}

    DataSet* velAC_;
    bool useVelInfo_; ///< If true use actual velocities in frame if present
    size_t frameIdx_;
    AtomMask mask_;
    Frame previousFrame_; ///< Hold previous frame coords (!useVelInfo only)
    typedef std::vector<Vec3> Varray;
    Varray Velocity0_; ///< Hold velocities at time 0.
    typedef std::vector<double> Darray;
    Darray Norm_; ///< Hold v0*v0 for normalization.
};
#endif
