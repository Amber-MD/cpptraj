#ifndef INC_RMSACTION_H
#define INC_RMSACTION_H
#include "ArgList.h"
#include "Topology.h"
/// Common functionality for RMS action
class RmsAction {
  public:
    RmsAction();
    /// Get nofit, norotate, and mass keywords.
    void GetRmsKeywords(ArgList&);
    /// Get target and reference masks.
    std::string GetRmsMasks(ArgList&);
    /// Print fit/rotate/mass status.
    void PrintRmsStatus();
    /// Set up target mask and frame for selected target atoms.
    int SetupRmsMask(Topology const&, const char*);
    AtomMask const& TgtMask() const { return tgtMask_; }
    bool Fit()                const { return fit_;     }
    bool UseMass()            const { return useMass_; }
    /// Perform fit/nofit [non-]mass-weighted calc with/without rotation.
    double CalcRmsd(Frame& TGT, Frame const& REF, Vec3 const& refTrans) {
      double R;
      // Set selected frame atoms. Masses have already been set.
      tgtFrame_.SetCoordinates(TGT, tgtMask_);
      if (!fit_) {
        R = tgtFrame_.RMSD_NoFit(REF, useMass_);
      } else {
        R = tgtFrame_.RMSD_CenteredRef(REF, rot_, tgtTrans_, useMass_);
        if (rotate_)
          TGT.Trans_Rot_Trans(tgtTrans_, rot_, refTrans);
        else {
          tgtTrans_ += refTrans;
          TGT.Translate(tgtTrans_);
        }
      }
      return R;
    }
  private:
    AtomMask tgtMask_; ///< Mask of selected target atoms.
    bool fit_;         ///< If true, best-fit RMS.
    bool rotate_;      ///< If true, rotate coordinates according to best-fit.
    bool useMass_;     ///< If true, mass-weight calculation.
    Vec3 tgtTrans_;    ///< Hold translation to origin.
    Matrix_3x3 rot_;   ///< Hold best-fit rotation matrix.
    Frame tgtFrame_;   ///< Hold selected target atoms.
};
#endif
