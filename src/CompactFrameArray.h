#ifndef INC_COMPACTFRAMEARRAY_H
#define INC_COMPACTFRAMEARRAY_H
#include "CoordinateInfo.h"
/// Hold a compact version of the time-depedent information contained in an array of Frames
class CompactFrameArray {
  public:
    CompactFrameArray();

    /// Set up frame array to hold coordinate info with specified # atoms, optionally # frames.
    int SetupFrameArray(CoordinateInfo const&, unsigned int, int);
  private:
    void addComponent(long int&, CoordinateInfo::Component, long int);

    typedef std::vector<float> Farray;
    Farray compactFrames_;                              ///< Array storing all info.
    bool hasComponent_[CoordinateInfo::NCOMPONENTS];    ///< True if array contains indicated component
    std::vector<CoordinateInfo::Component> components_; ///< List info contained in this array.
    std::vector<long int> offsets_;                     ///< Offsets for each component present.
    unsigned int frameSize_;                            ///< Size of each individual frame.
};
#endif
