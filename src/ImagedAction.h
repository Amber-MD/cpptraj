#ifndef INC_IMAGEDACTION_H
#define INC_IMAGEDACTION_H
#include "DistRoutines.h"
class ImagedAction {
  public:
    ImagedAction() : imageType_(NOIMAGE), useImage_(false) {}

    void InitImaging( bool imageIn ) { useImage_ = imageIn; }
  
    void SetupImaging( Box::BoxType parmboxtype ) {
      if (!useImage_)
        // Imaging disabled
        imageType_ = NOIMAGE;
      else {
        if (parmboxtype == Box::NOBOX) {
          imageType_ = NOIMAGE;
          //if (debug>0)
          //  mprintf("    Warning: No box info in %s, disabling imaging.\n",currentParm->c_str());
        } else if (parmboxtype == Box::ORTHO)
          imageType_ = ORTHO;
        else
          imageType_ = NONORTHO;
      }
    }

    /// Return true if imaging is currently enabled.
    bool ImagingEnabled()   { return (imageType_ != NOIMAGE); }
    bool UseImage()         { return useImage_;  } ///< True if imaging is desired.
    ImagingType ImageType() { return imageType_; } ///< Return type of imaging.
  private:
    ImagingType imageType_; ///< Type of imaging to be performed.
    bool useImage_;         ///< If true, use imaging.
};
#endif
