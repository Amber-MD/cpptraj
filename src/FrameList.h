#ifndef INC_FRAMELIST_H
#define INC_FRAMELIST_H
#include "TopologyList.h"
class ReferenceFrame {
  public:
    ReferenceFrame() : frame_(0), parm_(0), num_(0) {}
    ReferenceFrame(Frame* fIn, Topology* pIn, std::string const& nameIn, int nIn ) :
      frame_(fIn),
      parm_(pIn),
      name_(nameIn),
      num_(nIn)
    {}
    /*~ReferenceFrame() {
      if (frame_ != 0) delete frame_;
    }*/ // NOTE: Delete in FrameList
    //ReferenceFrame(const ReferenceFrame&);
    //ReferenceFrame& operator=(const ReferenceFrame&);
    bool operator==(const ReferenceFrame& rhs) const {
      if (frame_ != rhs.frame_) return false;
      if (parm_ != rhs.parm_) return false;
      return true;
    }
    void SetRef( Frame* newFrame, Topology* newParm ) {
      if (frame_ != 0) delete frame_;
      frame_ = newFrame;
      parm_ = newParm;
    }
    bool error()      const { return num_ == -1;    }
    bool empty()      const { return frame_ == 0;   }
    Frame* Coord()          { return frame_;        }
    Topology* Parm()        { return parm_;         }
    int Num()         const { return num_;          }
    std::string const& FrameName() const { return name_;   }
  private:
    Frame* frame_;     ///< Reference coords, allocated.
    Topology* parm_;   ///< Pointer to assiociated parm in TopologyList.
    std::string name_; ///< Unique name assigned to this ref structure (filename/tag).
    int num_;          ///< Frame number.
};

// Class: FrameList
/// Hold a series of Frame classes. 
/** Optionally hold a corresponding name and number (in the case of 
  * reference frames, this is the name of the reference coordinate 
  * file and the frame number). 
  */
// TODO: Eventually store a vector of Frames, not Frame*s
class FrameList : public FileList {
  public:
    FrameList();
    ~FrameList();
    void Clear();
    /// Return the current active reference frame
    Frame* ActiveReference();
    /// Set the active reference frame
    void SetActiveRef(int);
    /// Add a reference frame base on given args
    int AddReference(ArgList&, TopologyList&);
    /// Get reference frame based on given args
    ReferenceFrame GetFrameFromArgs(ArgList&) const;
    /// Get reference frame with given name.
    ReferenceFrame GetFrameByName(std::string const&) const;
    /// Replace the given reference frame with given Frame/Topology.
    int ReplaceFrame(ReferenceFrame const&, Frame*, Topology*);
    /// Print all reference frames.
    void List() const;
    /// \return the number of reference frames.
    int NumFrames() const { return (int)frames_.size(); }
  private:
    std::vector<ReferenceFrame> frames_;
    std::vector<Topology*> StrippedRefParms_;
    int refFrameNum_;
    static const ReferenceFrame ErrorFrame_;
};
#endif
