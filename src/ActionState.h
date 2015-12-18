#ifndef INC_ACTIONSTATE_H
#define INC_ACTIONSTATE_H
#include "DataFileList.h"
/*! \file ActionState.h
    \brief Classes used to wrap arguments passed to Actions.
 */
/** The ActionInit class is used to pass in the master DataSetList and
  * DataFileList.
  */
class ActionInit {
  public:
    ActionInit() : dsl_(0), dfl_(0) {} // NOTE: For pytraj/cython
    ActionInit(DataSetList& dslIn, DataFileList& dflIn) :
      dsl_(&dslIn), dfl_(&dflIn) {}
    DataSetList& DSL()              { return *dsl_; }
    DataSetList const& DSL()  const { return *dsl_; }
    /// Can be used by Actions that want to set up DataSets in Setup/DoAction/Print.
    DataSetList* DslPtr()           { return dsl_;  }
    DataFileList& DFL()             { return *dfl_; }
    DataFileList const& DFL() const { return *dfl_; }
  private:
    DataSetList* dsl_;
    DataFileList* dfl_;
};
/* The ActionSetup class is used to pass in the current Topology, CoordinateInfo,
 * and expected number of frames associated with the current Topology.
 */
class ActionSetup {
  public:
    // NOTE: Blank constructor/Set() for use during ensemble
    ActionSetup() : top_(0), cInfo_(0), nFrames_(0) {}
    void Set(Topology* p, CoordinateInfo const& c, int n) {
      top_ = p;
      cInfo_ = (CoordinateInfo*)&c;
      nFrames_ = n;
    }
    ActionSetup(Topology* topIn, CoordinateInfo const& cInfoIn, int nIn) :
      top_(topIn), cInfo_((CoordinateInfo*)&cInfoIn), nFrames_(nIn) {}
    Topology const& Top()             const { return *top_;    }
    /// Current Topology pointer.
    /** Used to set up output trajectories and for certain Actions so that
      * current Topology can be used in other functions (e.g. DoAction).
      */
    Topology* TopAddress()                  { return top_;     }
    CoordinateInfo const& CoordInfo() const { return *cInfo_;  }
    int Nframes()                     const { return nFrames_; }
    void SetTopology( Topology* p )         { top_ = p;        }
    void SetCoordInfo( CoordinateInfo* c )  { cInfo_ = c;      }
  private:
    Topology* top_;
    CoordinateInfo* cInfo_;
    int nFrames_;
};
/** The ActionFrame class is used to pass in the current Frame. */
class ActionFrame {
  public:
    ActionFrame() : frm_(0), trajoutNum_(0) {}
    ActionFrame(Frame* fIn, int t) : frm_(fIn), trajoutNum_(t) {}
    Frame const& Frm()  const { return *frm_; }
    Frame& ModifyFrm()        { return *frm_; }
    int TrajoutNum()    const { return trajoutNum_; }
    void SetTrajoutNum(int t) { trajoutNum_ = t; }
    // NOTE: Used during ensemble.
    Frame* FramePtr()         { return frm_;  }
    void SetFrame( Frame* f ) { frm_ = f;     }
  private:
    Frame* frm_;
    int trajoutNum_; ///< Current output trajectory frame number.
};
#endif
