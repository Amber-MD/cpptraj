#ifndef INC_ACTIONSTATE_H
#define INC_ACTIONSTATE_H
#include "DataFileList.h"
class ActionInit {
  public:
    ActionInit(DataSetList& dslIn, DataFileList& dflIn) :
      dsl_(&dslIn), dfl_(&dflIn) {}
    DataSetList& DSL()              { return *dsl_; }
    DataSetList const& DSL()  const { return *dsl_; }
    DataSetList* DslPtr()           { return dsl_;  }
    DataFileList& DFL()             { return *dfl_; }
    DataFileList const& DFL() const { return *dfl_; }
    DataSetList* DSL_Ptr()          { return dsl_;  }
  private:
    DataSetList* dsl_;
    DataFileList* dfl_;
};

class ActionSetup {
  public:
    // NOTE: Blank constructor/Set for use during ensemble
    ActionSetup() : top_(0), cInfo_(0), nFrames_(0) {}
    void Set(Topology* p, CoordinateInfo const& c, int n) {
      top_ = p;
      cInfo_ = (CoordinateInfo*)&c;
      nFrames_ = n;
    }
    ActionSetup(Topology* topIn, CoordinateInfo const& cInfoIn, int nIn) :
      top_(topIn), cInfo_((CoordinateInfo*)&cInfoIn), nFrames_(nIn) {}
    Topology const& Top()             const { return *top_;    }
    // NOTE: Used to set up output trajectories.
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

class ActionFrame {
  public:
    ActionFrame() : frm_(0) {}
    ActionFrame(Frame* fIn) : frm_(fIn) {}
    Frame const& Frm()  const { return *frm_; }
    Frame* FramePtr()         { return frm_;  }
    void SetFrame( Frame* f ) { frm_ = f;     }
  private:
    Frame* frm_;
};
#endif
