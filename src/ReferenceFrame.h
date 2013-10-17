#ifndef INC_REFERENCEFRAME_H
#define INC_REFERENCEFRAME_H
#include "Topology.h"
#include "ArgList.h"
class ReferenceFrame {
  public:
    ReferenceFrame() : frame_(0), parm_(0), num_(-1), strippedParm_(false) {}
    ~ReferenceFrame();
    Frame* Coord()                    { return frame_;      }
    Topology* Parm()                  { return parm_;       }
    bool error()                const { return num_ == -1;  }
    bool empty()                const { return frame_ == 0; }
    FileName const& FrameName() const { return name_;       }
    std::string const& Tag()    const { return tag_;        }
    int LoadRef(std::string const&, Topology*, int);
    int LoadRef(std::string const&, ArgList&, Topology*, int);
    int StripRef( std::string const& );
    int StripRef( AtomMask const& );
    void RefInfo() const;
  private:
    Frame* frame_;      ///< Reference coords, allocated.
    Topology* parm_;    ///< Pointer to assiociated parm in TopologyList.
    FileName name_;     ///< Ref structure filename.
    std::string tag_;   ///< Ref structure optional tag.
    int num_;           ///< Frame number.
    bool strippedParm_; ///< True if parm was stripped and should be deleted.
};
#endif
