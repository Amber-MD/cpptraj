#ifndef INC_IMAGEOPTION_H
#define INC_IMAGEOPTION_H
/// Used to determine if imaging should be used. Imaging disabled when no box present.
class ImageOption {
  public:
    /// CONSTRUCTOR
    ImageOption() : useImage_(false), imagingEnabled_(false) {}
    /// Determine if imaging is desired.
    void InitImaging(bool imagingDesired) { useImage_ = imagingDesired; }
    /// Determine if imaging possible if desired.
    void SetupImaging(bool hasBox) {
      if (useImage_) {
        // Imaging desired if possible.
        if (hasBox)
          imagingEnabled_ = true;
        else
          imagingEnabled_ = false;
      } else
        // No imaging desired.
        imagingEnabled_ = false;
    }
    /// \return True if imaging is desired.
    bool UseImage() const { return useImage_; }
    /// \return True if imaging desired and currently possible.
    bool ImagingEnabled() const { return imagingEnabled_; }
  private:
    bool useImage_;       ///< If true, imaging desired if possible.
    bool imagingEnabled_; ///< If true, imaging is possible and should be on.
};
#endif
