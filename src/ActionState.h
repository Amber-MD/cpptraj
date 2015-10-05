#ifndef INC_ACTIONSTATE_H
#define INC_ACTIONSTATE_H
#include "DataFileList.h"
class ActionInit {
  public:
    ActionInit(DataSetList& dslIn, DataFileList& dflIn) :
      dsl_(&dslIn), dfl_(&dflIn) {}
    DataSetList& DSL()              { return *dsl_; }
    DataSetList const& DSL()  const { return *dsl_; }
    DataFileList& DFL()             { return *dfl_; }
    DataFileList const& DFL() const { return *dfl_; }
    DataSetList* DSL_Ptr()          { return dsl_;  }
  private:
    DataSetList* dsl_;
    DataFileList* dfl_;
};

class ActionSetup {
  public:
    ActionSetup(Topology* topIn, CoordinateInfo* cInfoIn, int nIn) :
      top_(topIn), cInfo_(cInfoIn), nFrames_(nIn) {}
    Topology const& Top()             const { return *top_;    }
    //Topology* TopPtr()                      { return top_;    }
    CoordinateInfo const& CoordInfo() const { return *cInfo_;  }
    int Nframes()                     const { return nFrames_; }
    //CoordinateInfo* CinfoPtr()              { return cInfo_;  }
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
    void SetFrame( Frame* f ) { frm_ = f;     }
  private:
    Frame* frm_;
};
#endif
