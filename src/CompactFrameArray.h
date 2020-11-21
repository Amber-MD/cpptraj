#ifndef INC_COMPACTFRAMEARRAY_H
#define INC_COMPACTFRAMEARRAY_H
#include "CoordinateInfo.h"
#include <algorithm>
/// Hold a compact version of the time-depedent information contained in an array of Frames
template <class T> class CompactFrameArray {
  public:
    CompactFrameArray() { std::fill(hasComponent_, hasComponent_+CoordinateInfo::NCOMPONENTS, false); }
  private:
    typedef std::vector<T> ArrayType;
    std::vector<ArrayType> compactFrames_;              ///< Array storing all info.
    bool hasComponent_[CoordinateInfo::NCOMPONENTS];    ///< True if array contains indicated component
    std::vector<CoordinateInfo::Component> components_; ///< List info contained in this array.
    unsigned int frameSize_;                            ///< Size of each individual frame.
    std::vector<long int> offsets_;                     ///< Offsets for each component present.
};
#endif
