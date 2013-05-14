// -*- mode: c++; -*-
#ifndef INC_ACTION_PAIRDIST_H
#define INC_ACTION_PAIRDIST_H
#include "Action.h"
#include "ImagedAction.h"

/** \author Hannes H. Loeffler
  */
class Action_PairDist : public Action, ImagedAction {
 public:
  Action_PairDist();

  static DispatchObject* Alloc() {
    return (DispatchObject*)new Action_PairDist();
  }

  static void Help();

 private:
  Action::RetType Init(ArgList&, TopologyList*, FrameList*, DataSetList*,
		       DataFileList*, int);
  Action::RetType Setup(Topology*, Topology**);
  Action::RetType DoAction(int, Frame*, Frame**);
  void Print();

  template <class T>
  class OnlineVar {
    public:
      OnlineVar() : n_(0.0), mean_(0.0), M2_(0.0) {}
      void accumulate(const T x)
      {
        T delta;
        n_++;
        delta = x - mean_;
        mean_ += delta / n_;
        M2_ += delta * (x - mean_);
      }
      T mean() const { return mean_; };
      T variance() const { return M2_ / (n_ - 1); };
      T nData() const { return n_; };
    private:
      T n_;
      T mean_;
      T M2_;
  };

  CpptrajFile output_;

  AtomMask mask1_;
  AtomMask mask2_;

  double delta_;		///< resolution

  std::vector<OnlineVar<double> > histogram_;
  unsigned long maxbin_;

  bool same_mask_;
  unsigned long ub1_;
  unsigned long ub2_;
};
#endif
