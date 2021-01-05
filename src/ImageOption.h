#ifndef INC_IMAGEOPTION_H
#define INC_IMAGEOPTION_H
/// Used to determine if imaging should be used. Imaging disabled when no box present.
class ImageOption {
  public:
    /// Different imaging types.
    enum Type { NO_IMAGE = 0, ORTHO, NONORTHO };
    /// CONSTRUCTOR
    ImageOption() : type_(NO_IMAGE), useImage_(false), imagingEnabled_(false), forceNonortho_(false) {}
    /// Determine if imaging is desired. Optionally force nonortho.
    void InitImaging(bool imagingDesired, bool forceNonOrthoIn) {
      useImage_ = imagingDesired;
      forceNonortho_ = forceNonOrthoIn;
      type_ = NO_IMAGE;
    }
    /// Determine if imaging is desired.
    void InitImaging(bool imagingDesired) { InitImaging(imagingDesired, false); }
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
    bool UseImage()       const { return useImage_; }
    /// \return True if imaging desired and currently possible.
    bool ImagingEnabled() const { return imagingEnabled_; }
    /// \return Current imaging type
    Type ImagingType()    const { return type_; }
    /// Determine imaging type given current options and if cell is x_aligned orthogonal
    void SetImageType(bool is_x_aligned_ortho) {
      if (imagingEnabled_) {
        if (forceNonortho_)
          type_ = NONORTHO;
        else if (is_x_aligned_ortho)
          type_ = ORTHO;
        else
          type_ = NONORTHO;
      } else
        type_ = NO_IMAGE;
    }
  private:
    Type type_;           ///< Current imaging type based on current cell.
    bool useImage_;       ///< If true, imaging desired if possible.
    bool imagingEnabled_; ///< If true, imaging is possible and should be on.
    bool forceNonortho_;  ///< If true, non-ortho imaging will be used if possible.
};
#endif
