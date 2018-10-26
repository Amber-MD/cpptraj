#include "Exec_UpdateParameters.h"
#include "CpptrajStdio.h"
#include "DataSet_Parameters.h"
#include "DataSet_Topology.h"

// Exec_UpdateParameters::Help()
void Exec_UpdateParameters::Help() const
{
  mprintf("\t%s setname <parm set>\n", DataSetList::TopArgs);
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
  return updateCount;
}

/** Update parameter set of given Topology with parameters from new set.
  * \param top The Topology to update.
  * \param set1 Set with new/updated parameters.
  */
int Exec_UpdateParameters::UpdateParams(Topology& top, ParameterSet const& set1) const
{
  ParameterSet set0 = top.GetParameters();
  // Check
  if (set0.AT().size() < 1)
    mprintf("Warning: No atom type information in '%s'\n", top.c_str());
  set0.Debug("originalp.dat");

  unsigned int updateCount;
  // Bond parameters
  updateCount = UpdateParameters< ParmHolder<BondParmType> >(set0.BP(), set1.BP(), "bond");
  if (updateCount > 0) {
    mprintf("\tRegenerating bond parameters.\n");
    top.AssignBondParams( set0.BP() );
  }
  // Angle parameters
  updateCount = UpdateParameters< ParmHolder<AngleParmType> >(set0.AP(), set1.AP(), "angle");
  if (updateCount > 0) {
    mprintf("\tRegenerating angle parameters.\n");
    top.AssignAngleParams( set0.AP() );
  }
  // Dihedral parameters
  updateCount = UpdateParameters< DihedralParmHolder >(set0.DP(), set1.DP(), "dihedral");
  if (updateCount > 0) {
    mprintf("\tRegenerating dihedral parameters.\n");
    top.AssignDihedralParams( set0.DP() );
  }
  // Urey-Bradley
  updateCount = UpdateParameters< ParmHolder<BondParmType> >(set0.UB(), set1.UB(), "Urey-Bradley");
  if (updateCount > 0) {
    mprintf("\tRegenerating UB parameters.\n");
    top.AssignUBParams( set0.UB() );
  }
  // Improper parameters
  updateCount = UpdateParameters< ParmHolder<DihedralParmType> >(set0.IP(), set1.IP(), "improper");
  if (updateCount > 0) {
    mprintf("\tRegenerating improper parameters.\n");
    top.AssignImproperParams( set0.IP() );
  }
  // Atom types
  updateCount = UpdateParameters< ParmHolder<AtomType> >(set0.AT(), set1.AT(), "atom type");
  updateCount += UpdateParameters< ParmHolder<NonbondType> >(set0.NB(), set1.NB(), "LJ A-B");
  if (updateCount > 0) {
    mprintf("\tRegenerating nonbond parameters.\n");
    top.AssignNonbondParams( set0.AT(), set0.NB() );
  }

  set0.Debug("newp.dat");
  return 0;
}

// Exec_UpdateParameters::Execute()
Exec::RetType Exec_UpdateParameters::Execute(CpptrajState& State, ArgList& argIn)
{
  std::string dsname = argIn.GetStringKey("setname");
  if (dsname.empty()) {
    mprinterr("Error: Specify parameter set.\n");
    return CpptrajState::ERR;
  }
  DataSet* ds = State.DSL().GetDataSet( dsname );
  if (ds == 0) {
    mprinterr("Error: Parameter data set '%s' not found.\n", dsname.c_str());
    return CpptrajState::ERR;
  }
  if (ds->Type() != DataSet::PARAMETERS && ds->Type() != DataSet::TOPOLOGY) {
    mprinterr("Error: Set '%s' is not a parameter or topology data set.\n", ds->legend());
    return CpptrajState::ERR;
  }
  Topology* dstop = State.DSL().GetTopology( argIn );
  if (dstop == 0) {
    mprinterr("Error: No topology specified.\n");
    return CpptrajState::ERR;
  }
  Topology& top = static_cast<Topology&>( *dstop );

  mprintf("\tUpdating parameters in topology '%s' using those in set '%s'\n",
          top.c_str(), ds->legend());

  if (ds->Type() == DataSet::PARAMETERS)
    UpdateParams(top, static_cast<DataSet_Parameters const&>( *ds ));
  else if (ds->Type() == DataSet::TOPOLOGY) {
    DataSet_Topology const& topds = static_cast<DataSet_Topology const&>( *ds );
    UpdateParams(top, topds.Top().GetParameters());
  } else // Sanity check
    return CpptrajState::ERR;

  return CpptrajState::OK;
}
