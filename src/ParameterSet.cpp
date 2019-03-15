#include "ParameterSet.h"
#include "CpptrajFile.h"
#include "CpptrajStdio.h"

size_t ParameterSet::DataSize() const {
  return (atomTypes_.DataSize() +
          bondParm_.DataSize() +
          angleParm_.DataSize() +
          ubParm_.DataSize() +
          dihParm_.DataSize() +
          impParm_.DataSize());
}

void ParameterSet::Debug(const char* fnameIn) const {
  CpptrajFile Out;
  Out.OpenWrite( fnameIn );
  Out.Printf("Atom Types:\n");
  Out.Printf("\t%6s %8s %12s %12s %12s\n", "Name", "TypeIdx", "Radius", "Depth", "Mass");
  for (ParmHolder<AtomType>::const_iterator at = atomTypes_.begin(); at != atomTypes_.end(); ++at) {
    Out.Printf("\t%6s %8i %12.4f %12.4f %12.4f\n", *(at->first[0]), at - atomTypes_.begin(), at->second.LJ().Radius(), at->second.LJ().Depth(), at->second.Mass());
  }
  if (!nbParm_.empty()) {
    Out.Printf("LJ parameters:\n");
    Out.Printf("\t%6s %6s : %12s %12s\n", "Type1", "Type2", "A", "B");
    for (ParmHolder<NonbondType>::const_iterator nb = nbParm_.begin(); nb != nbParm_.end(); ++nb)
      Out.Printf("\t%6s %6s : %12.4E %12.4E\n", (*nb->first[0]), (*nb->first[1]), nb->second.A(), nb->second.B());
  }
  if (!bondParm_.empty()) {
    Out.Printf("Bond parameters:\n");
    Out.Printf("\t%6s %6s : %12s %12s\n", "Type1", "Type2", "Rk", "Req");
    for (ParmHolder<BondParmType>::const_iterator bp = bondParm_.begin(); bp != bondParm_.end(); ++bp)
      Out.Printf("\t%6s %6s : %12.4f %12.4f\n", *(bp->first[0]), *(bp->first[1]), bp->second.Rk(), bp->second.Req());
  }
  if (!angleParm_.empty()) {
    Out.Printf("Angle parameters:\n");
    Out.Printf("\t%6s %6s %6s : %12s %12s\n", "Type1", "Type2", "Type3", "Tk", "Teq");
    for (ParmHolder<AngleParmType>::const_iterator bp = angleParm_.begin(); bp != angleParm_.end(); ++bp)
      Out.Printf("\t%6s %6s %6s : %12.4f %12.4f\n", *(bp->first[0]), *(bp->first[1]), *(bp->first[2]), bp->second.Tk(), bp->second.Teq());
  }
  if (!ubParm_.empty()) {
    Out.Printf("UB parameters:\n");
    Out.Printf("\t%6s %6s : %12s %12s\n", "Type1", "Type2", "Uk", "Ueq");
    for (ParmHolder<BondParmType>::const_iterator bp = ubParm_.begin(); bp != ubParm_.end(); ++bp)
      Out.Printf("\t%s %s : %f %f\n", *(bp->first[0]), *(bp->first[1]), bp->second.Rk(), bp->second.Req());
  }
  if (!dihParm_.empty()) {
    Out.Printf("Dihedral parameters:\n");
    Out.Printf("\t%6s %6s %6s %6s %12s %12s %12s\n", "Type1", "Type2", "Type3", "Type4", "Pk", "Pn", "Phase");
    for (DihedralParmHolder::const_iterator it0 = dihParm_.begin(); it0 != dihParm_.end(); ++it0)
      for (DihedralParmArray::const_iterator it1 = it0->second.begin();
                                             it1 != it0->second.end(); ++it1)
        Out.Printf("\t%6s %6s %6s %6s : %12.4f %12.4f %12.4f\n", *(it0->first[0]), *(it0->first[1]), *(it0->first[2]), *(it0->first[3]), it1->Pk(), it1->Pn(), it1->Phase());
  }
  if (!impParm_.empty()) {
    Out.Printf("Improper parameters:\n");
    Out.Printf("\t%6s %6s %6s %6s %12s %12s %12s\n", "Type1", "Type2", "Type3", "Type4", "Pk", "Pn", "Phase");
    for (ParmHolder<DihedralParmType>::const_iterator bp = impParm_.begin(); bp != impParm_.end(); ++bp)
      Out.Printf("\t%6s %6s %6s %6s : %12.4f %12.4f %12.4f\n", *(bp->first[0]), *(bp->first[1]), *(bp->first[2]), *(bp->first[3]), bp->second.Pk(), bp->second.Pn(), bp->second.Phase());
  }
}

static inline void PrintParmType(BondParmType const& bp) { mprintf(" %12.4f %12.4f\n", bp.Rk(), bp.Req()); }
static inline void PrintParmType(AngleParmType const& ap) { mprintf(" %12.4f %12.4f\n", ap.Tk(), ap.Teq()); }
static inline void PrintParmType(DihedralParmType const& dp) { mprintf(" %12.4f %12.4f %12.4f\n", dp.Pk(), dp.Pn(), dp.Phase()); }
static inline void PrintParmType(DihedralParmArray const& dpa) {
  mprintf("\n");
  for (DihedralParmArray::const_iterator it = dpa.begin(); it != dpa.end(); ++it)
    mprintf("\t\t%12.4f %12.4f %12.4f\n", it->Pk(), it->Pn(), it->Phase());
}
static inline void PrintParmType(AtomType const& at) { mprintf(" %12.4f %12.4f %12.4f\n", at.LJ().Radius(), at.LJ().Depth(), at.Mass()); }
static inline void PrintParmType(NonbondType const& nb) { mprintf(" %12.4E %12.4E\n", nb.A(), nb.B()); }

/** Add update parameters.
  * \param0 Parameters to add to/update.
  * \param1 New parameters.
  * \param desc Description of parameters.
  */
template <typename T> int UpdateParameters(T& param0, T const& param1, const char* desc)
{
  // DEBUG
//  mprintf("DEBUG: Current %s Parameters:\n", desc);
//  for (typename T::const_iterator p = param0.begin(); p != param0.end(); ++p)
//    PrintParmType( p->second );

  int updateCount = 0;
  for (typename T::const_iterator newp = param1.begin(); newp != param1.end(); ++newp)
  {
    ParameterHolders::RetType ret = param0.AddParm( newp->first, newp->second, true );
    if (ret != ParameterHolders::ERR) {
      if (ret == ParameterHolders::ADDED) {
        mprintf("\tAdded NEW %s parameter:", desc);
        updateCount++;
      } else if (ret == ParameterHolders::UPDATED) {
        mprintf("\tUpdated %s parameter:", desc);
        updateCount++;
      } else if (ret == ParameterHolders::SAME)
        mprintf("\tParameter for %s already present:", desc);
      mprintf(" %s", newp->first.TypeString().c_str());
      PrintParmType( newp->second );
      //mprintf(" %s %s %12.4f %12.4f\n", 
      //        *(newp->first[0]), *(newp->first[1]), newp->second.Rk(), newp->second.Req());
    }
  }
  // DEBUG
//  mprintf("DEBUG: New %s Parameters:\n", desc);
//  for (typename T::const_iterator p = param0.begin(); p != param0.end(); ++p)
//    PrintParmType( p->second );

  return updateCount;
}

/** Update/add to parameters in this topology with those from given set. */
int ParameterSet::UpdateParams(ParameterSet const& set1, UpdateCount& uc) {
  ParameterSet& set0 = *this;
  // Check
  set0.Debug("originalp.dat");
  set1.Debug("incomingp.dat");

  // Bond parameters
  uc.nBondsUpdated_ = UpdateParameters< ParmHolder<BondParmType> >(set0.BP(), set1.BP(), "bond");
  // Angle parameters
  uc.nAnglesUpdated_ = UpdateParameters< ParmHolder<AngleParmType> >(set0.AP(), set1.AP(), "angle");
  // Dihedral/improper parameters
  uc.nDihedralsUpdated_ = UpdateParameters< DihedralParmHolder >(set0.DP(), set1.DP(), "dihedral");
  // Improper parameters
  uc.nImpropersUpdated_ = UpdateParameters< ParmHolder<DihedralParmType> >(set0.IP(), set1.IP(), "improper");
  // Urey-Bradley parameters
  uc.nUreyBradleyUpdated_ = UpdateParameters< ParmHolder<BondParmType> >(set0.UB(), set1.UB(), "Urey-Bradley");
  // Atom types
  uc.nAtomTypeUpdated_ = UpdateParameters< ParmHolder<AtomType> >(set0.AT(), set1.AT(), "atom type");
  // LJ Pairs
  uc.nLJparamsUpdated_ = UpdateParameters< ParmHolder<NonbondType> >(set0.NB(), set1.NB(), "LJ A-B");

  set0.Debug("newp.dat");
  return 0;
}
