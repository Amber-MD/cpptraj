// -*- mode: c++; -*-
#ifndef INC_ACTION_PAIRDIST_H
#define INC_ACTION_PAIRDIST_H

#include "Action.h"
#include "ImageOption.h"
#include "OnlineVarT.h"
/** \author Hannes H. Loeffler
  * Updated by DRR to ensure results do not change depending on # processes
  * used and/or initial max size of the histogram.
  */
class Action_PairDist : public Action {
 public:
  Action_PairDist();

    DispatchObject* Alloc() const { return (DispatchObject*)new Action_PairDist(); }

    void Help() const;

 private:
  Action::RetType Init(ArgList&, ActionInit&, int);
  Action::RetType Setup(ActionSetup&);
  Action::RetType DoAction(int, ActionFrame&);
# ifdef MPI
  int SyncAction();
  Parallel::Comm trajComm_;
# endif
  void Print();

  inline void BinHist(std::vector<double>&, const double*, const double*, Box const&);
  void UpdateHistogramFrames();

  DataSet *Pr_;    ///< distance vs P(r)
  DataSet *std_;   ///< distance vs std
  AtomMask mask1_;
  AtomMask mask2_;
  double delta_;                          ///< histogram resolution
  std::vector<Stats<double> > histogram_; ///< The final histogram
  unsigned long maxbin_;                  ///< Index of current maximum histogram bin
  ImageOption imageOpt_;                  ///< Used to decide if imaging should be used
  bool single_mask_;                      ///< True if only 1 mask specified, false if 2
  unsigned int nframes_;                  ///< Used to update ndata in bins not present at start
  //int frame_, idx1_, idx2_; // DEBUG
};
#endif
