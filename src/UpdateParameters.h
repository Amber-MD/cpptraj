#ifndef INC_UPDATEPARAMETERS_H
#define INC_UPDATEPARAMETERS_H
// NOTE: By design, this is intended for inclusion into the body of other
//       classes, e.g. ParameterSet and/or Topology. Therefore there are
//       no 'include' statements or forward declares as it is assumed
//       those files have already been included prior to this one.
static inline void PrintParmType(BondParmType const& bp) { mprintf(" %12.4f %12.4f\n", bp.Rk(), bp.Req()); }
static inline void PrintParmType(AngleParmType const& ap) { mprintf(" %12.4f %12.4f\n", ap.Tk(), ap.Teq()*Constants::RADDEG); }
static inline void PrintParmType(DihedralParmType const& dp) { mprintf(" %12.4f %12.4f %12.4f\n", dp.Pk(), dp.Pn(), dp.Phase()*Constants::RADDEG); }
static inline void PrintParmType(DihedralParmArray const& dpa) {
  mprintf("\n");
  for (DihedralParmArray::const_iterator it = dpa.begin(); it != dpa.end(); ++it)
    mprintf("\t\t%12.4f %12.4f %12.4f\n", it->Pk(), it->Pn(), it->Phase()*Constants::RADDEG);
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
#endif
