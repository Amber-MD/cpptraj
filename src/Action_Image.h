#ifndef INC_ACTION_IMAGE_H
#define INC_ACTION_IMAGE_H
// Class: Action_Image
/// Action to wrap coordinates back into primary box
#include "Action.h"
class Action_Image: public Action {
  public:
    Action_Image();

    static DispatchObject* Alloc() { return (DispatchObject*)new Action_Image(); }
    static void Help();

    ~Action_Image();
    void Print() {}
  private:
    enum ImageMode { BYMOL = 0, BYRES, BYATOM };
    static const char* ImageModeString[];
    ImageMode imageMode_;
    /// Only atoms in Mask1 will be imaged
    AtomMask Mask1_;
    /// If defined, image w.r.t. the center of atoms in ComMask.
    AtomMask *ComMask_;
    /// If true image w.r.t. coordinate origin, otherwise box center
    bool origin_;
    /// If true molecules will be imaged w.r.t. their center, otherwise first atom will be used
    bool center_;
    /// True if orthorhombic cell, false otherwise.
    bool ortho_;
    bool useMass_;
    bool truncoct_;
    enum TriclinicArg {OFF, FORCE, FAMILIAR};
    TriclinicArg triclinic_;
    int debug_;
    /// Vector containing atom ranges to be imaged (first to last)
    std::vector<int> imageList_; 

    Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
                          DataFileList*, int);
    Action::RetType Setup(Topology*, Topology**);
    Action::RetType DoAction(int, Frame*, Frame**);
    void CheckRange(int,int);
};
#endif
