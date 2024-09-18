#ifndef INC_PAIRLISTENGINE_H
#define INC_PAIRLISTENGINE_H
#include "PairList.h"
class Topology;
class AtomMask;
namespace Cpptraj {
class EngineOpts;
/// Abstract base class for Pairlist engines
class PairListEngine {
  public:
    PairListEngine() {}
    // Virtual since inherited
    virtual ~PairListEngine() {}

    virtual int InitEngine(Box const&, EngineOpts const&, int) = 0;
    virtual int SetupEngine(Topology const&, AtomMask const&) = 0;

    virtual void FrameBeginCalc() = 0;
    virtual void SetupAtom0(PairList::AtmType const&) = 0;
    virtual void SetupAtom1(PairList::AtmType const&) = 0;
    virtual void CutoffSatisfied(double, PairList::AtmType const&, PairList::AtmType const&) = 0;
    virtual void CutoffNotSatisfied(double, PairList::AtmType const&, PairList::AtmType const&) = 0;
};
}
#endif
