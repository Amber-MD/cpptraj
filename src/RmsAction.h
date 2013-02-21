#ifndef INC_RMSACTION_H
#define INC_RMSACTION_H
#include "ArgList.h"
#include "Topology.h"
/// Common functionality for RMS action
class RmsAction {
  public:
    RmsAction();
    void GetRmsKeywords(ArgList&);
    std::string GetRmsMasks(ArgList&);
    void PrintRmsStatus();
    int SetupRmsMask(Topology const&, const char*);
    AtomMask const& TgtMask() const { return tgtMask_; }
    bool Fit()                const { return fit_;     }
    bool UseMass()            const { return useMass_; }
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
    AtomMask tgtMask_;
    bool fit_;
    bool rotate_;
    bool useMass_;
    Vec3 tgtTrans_;
    Matrix_3x3 rot_;
    Frame tgtFrame_;
};
#endif
