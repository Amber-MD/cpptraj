#ifndef INC_COMPACTFRAMEARRAY_H
#define INC_COMPACTFRAMEARRAY_H
#include "CoordinateInfo.h"
/// Hold a compact version of the time-depedent information contained in an array of Frames
class CompactFrameArray {
  public:
    CompactFrameArray();

    /// \return True if set up for any components
    bool HasComponents() const { return !(components_.empty()); }

    /// Allocate for specified number of frames
    void Resize(int);
    /// Set up frame array to hold coordinate info with specified # atoms, optionally # frames.
    int SetupFrameArray(CoordinateInfo const&, unsigned int, int);

    /// Copy component from double array to specified frame
    int SetFromDblPtr(unsigned int, const double*, CoordinateInfo::Component);
    /// Copy component from specified frame to double array
    int GetToDblPtr(double*, unsigned int, CoordinateInfo::Component) const;
    /// Copy component from integer array to specified frame
    int SetFromIntPtr(unsigned int, const int*, CoordinateInfo::Component);
    /// Copy component from specified frame to integer array
    int GetToIntPtr(int*, unsigned int, CoordinateInfo::Component) const;

    /// Copy component from double value to specified frame
    int SetFromDblVal(unsigned int, double, CoordinateInfo::Component);
    /// \return component for specified frame (assumes size is 1)
    float GetVal(unsigned int, CoordinateInfo::Component) const;
  private:
    void addComponent(long int&, CoordinateInfo::Component, long int);

    typedef std::vector<float> Farray;
    Farray compactFrames_;                              ///< Array storing all info.
    int componentIdx_[CoordinateInfo::NCOMPONENTS];     ///< Index into arrays for each component
    std::vector<CoordinateInfo::Component> components_; ///< List info contained in this array.
    std::vector<long int> offsets_;                     ///< Offsets for each component present.
    //unsigned int frameSize_;                            ///< Size of each individual frame.
};
#endif
