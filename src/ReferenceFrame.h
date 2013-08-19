#ifndef INC_REFERENCEFRAME_H
#define INC_REFERENCEFRAME_H
class ReferenceFrame {
  public:
    ReferenceFrame() : frame_(0), parm_(0), num_(0) {}
    ReferenceFrame(Frame* fIn, Topology* pIn, FileName const& nameIn,
                   std::string const& tagIn, int nIn ) :
      frame_(fIn),
      parm_(pIn),
      name_(nameIn),
      tag_(tagIn),
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
    Frame* Coord()                    { return frame_;      }
    Topology* Parm()                  { return parm_;       }
    bool error()                const { return num_ == -1;  }
    bool empty()                const { return frame_ == 0; }
    int Num()                   const { return num_;        }
    FileName const& FrameName() const { return name_;       }
    std::string const& Tag()    const { return tag_;        }
  private:
    Frame* frame_;     ///< Reference coords, allocated.
    Topology* parm_;   ///< Pointer to assiociated parm in TopologyList.
    FileName name_;    ///< Ref structure filename.
    std::string tag_;  ///< Ref structure optional tag.
    int num_;          ///< Frame number.
};
#endif
