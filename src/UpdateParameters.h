#ifndef INC_UPDATEPARAMETERS_H
#define INC_UPDATEPARAMETERS_H
// NOTE: By design, this is intended for inclusion into the body of other
//       classes, e.g. ParameterSet and/or Topology. Therefore there are
//       no 'include' statements or forward declares as it is assumed
//       those files have already been included prior to this one.
static inline void PrintParmType(BondParmType const& bp) { mprintf(" %12.4f %12.4f", bp.Rk(), bp.Req()); }
static inline void PrintParmType(AngleParmType const& ap) { mprintf(" %12.4f %12.4f", ap.Tk(), ap.Teq()*Constants::RADDEG); }
static inline void PrintParmType(DihedralParmType const& dp) { mprintf(" %12.4f %12.4f %12.4f", dp.Pk(), dp.Pn(), dp.Phase()*Constants::RADDEG); }
static inline void PrintParmType(DihedralParmArray const& dpa) {
  mprintf("\n");
  for (DihedralParmArray::const_iterator it = dpa.begin(); it != dpa.end(); ++it)
    mprintf("\t\t%12.4f %12.4f %12.4f\n", it->Pk(), it->Pn(), it->Phase()*Constants::RADDEG);
}
static inline void PrintParmType(AtomType const& at) { mprintf(" R=%g Depth=%g Mass=%g", at.LJ().Radius(), at.LJ().Depth(), at.Mass()); }
static inline void PrintParmType(NonbondType const& nb) { mprintf(" %12.4E %12.4E", nb.A(), nb.B()); }
static inline void PrintParmType(HB_ParmType const& hb) { mprintf(" %12.4E %12.4E %12.4E", hb.Asol(), hb.Bsol(), hb.HBcut()); }

/** Add update parameters.
  * \param0 Parameters to add to/update.
  * \param1 New parameters.
  * \param desc Description of parameters.
  * \param verbose Verbosity: 0 - silent,  1 - updated only, 2 - updated & same, 3 - all
  */
template <typename T> int UpdateParameters(T& param0, T const& param1, const char* desc, int verbose)
{
  // DEBUG
//  mprintf("DEBUG: Current %s Parameters:\n", desc);
//  for (typename T::const_iterator p = param0.begin(); p != param0.end(); ++p)
//    PrintParmType( p->second );

  int updateCount = 0;
  for (typename T::const_iterator newp = param1.begin(); newp != param1.end(); ++newp)
  {
    Cpptraj::Parm::RetType ret = param0.AddParm( newp->first, newp->second, true );
    if (ret != Cpptraj::Parm::ERR) {
      bool print = false;
      if (ret == Cpptraj::Parm::ADDED) {
        if (verbose > 2) { mprintf("\tAdded NEW %s parameter:", desc); print = true; }
        updateCount++;
      } else if (ret == Cpptraj::Parm::UPDATED) {
        if (verbose > 0) { mprintf("\tUpdated %s parameter:", desc); print = true; }
        updateCount++;
      } else if (ret == Cpptraj::Parm::SAME) {
        if (verbose > 1) { mprintf("\tParameter for %s already present:", desc); print = true; }
      }
      if (print) {
        mprintf(" %s", newp->first.TypeString().c_str());
        //PrintParmType( newp->second );
        // NOTE: For AtomType its possible only a partial update occurred.
        //       Get the new parameter explicitly.
        typename T::const_iterator newIt = param0.GetParam( newp->first );
        if (ret == Cpptraj::Parm::UPDATED) {
          mprintf(" From");
          PrintParmType( param0.PreviousParm() );
          mprintf(" to");
          PrintParmType( newIt->second );
        } else
          PrintParmType( newIt->second );
        mprintf("\n");
      }
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
