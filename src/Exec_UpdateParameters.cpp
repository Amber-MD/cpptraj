#include "Exec_UpdateParameters.h"
#include "CpptrajStdio.h"
#include "DataSet_Parameters.h"

// Exec_UpdateParameters::Help()
void Exec_UpdateParameters::Help() const
{
  mprintf("\t%s setname <parm set>\n", DataSetList::TopArgs);
}

/** Add/update bond parameters.
  * \param bp0 Bond parameters to add to/update.
  * \param bp1 New bond parameters.
  * \return number of bond parameters added or updated.
  */
int Exec_UpdateParameters::UpdateBondParams(ParmHolder<BondParmType>& bp0,
                                            ParmHolder<BondParmType> const& bp1) const
{
  int updateCount = 0;
  for (ParmHolder<BondParmType>::const_iterator newp = bp1.begin();
                                                newp != bp1.end(); ++newp)
  {
    ParameterHolders::RetType ret = bp0.AddParm( newp->first, newp->second, true );
    if (ret != ParameterHolders::ERR) {
      if (ret == ParameterHolders::ADDED) {
        mprintf("\tAdded NEW bond parameter:");
        updateCount++;
      } else if (ret == ParameterHolders::UPDATED) {
        mprintf("\tUpdated bond parameter:");
        updateCount++;
      } else if (ret == ParameterHolders::SAME)
        mprintf("\tBond parameter already present:");
      mprintf(" %s %s %12.4f %12.4f\n", 
              *(newp->first[0]), *(newp->first[1]), newp->second.Rk(), newp->second.Req());
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
  set0.Debug("originalp.dat");

  //StatArray status;
  //status = UpdateParamX<ParmHolder<BondParmType>>(set0.BP(), set1.BP());
  unsigned int updateCount;
  // Bond parameters
  updateCount = UpdateBondParams(set0.BP(), set1.BP());
  if (updateCount > 0) {
    mprintf("\tRegenerating bond parameters.\n");
    top.AssignBondParams( set0.BP() );
  }
  // Angle parameters
  updateCount = 0;
  for (ParmHolder<AngleParmType>::const_iterator newp = set1.AP().begin();
                                                 newp != set1.AP().end(); ++newp)
  {
    ParameterHolders::RetType ret = set0.AP().AddParm( newp->first, newp->second, true );
    if (ret != ParameterHolders::ERR) {
      if (ret == ParameterHolders::ADDED) {
        mprintf("\tAdded NEW angle parameter:");
        updateCount++;
      } else if (ret == ParameterHolders::UPDATED) {
        mprintf("\tUpdated angle parameter:");
        updateCount++;
      } else if (ret == ParameterHolders::SAME)
        mprintf("\tAngle parameter already present:");
      mprintf(" %s %s %s %12.4f %12.4f\n", 
              *(newp->first[0]), *(newp->first[1]), *(newp->first[2]), newp->second.Tk(), newp->second.Teq());
    }
  }
  if (updateCount > 0) {
    mprintf("\tRegenerating angle parameters.\n");
    top.AssignAngleParams( set0.AP() );
  }
  // Dihedral parameters
  updateCount = 0;
  for (DihedralParmHolder::const_iterator newp = set1.DP().begin();
                                          newp != set1.DP().end(); ++newp)
  {
    ParameterHolders::RetType ret = set0.DP().AddParm( newp->first, newp->second, true );
    if (ret != ParameterHolders::ERR) {
      if (ret == ParameterHolders::ADDED) {
        mprintf("\tAdded NEW dihedral parameters:\n");
        updateCount++;
      } else if (ret == ParameterHolders::UPDATED) {
        mprintf("\tUpdated dihedral parameters:\n");
        updateCount++;
      } else if (ret == ParameterHolders::SAME)
        mprintf("\tDihedral parameters already present:\n");
      for (DihedralParmArray::const_iterator it = newp->second.begin();
                                             it != newp->second.end(); ++it)
        mprintf("\t\t%s %s %s %s %12.4f %12.4f %12.4f\n", 
                *(newp->first[0]), *(newp->first[1]), *(newp->first[2]), *(newp->first[3]), it->Pk(), it->Pn(), it->Phase());
    }
  }
  if (updateCount > 0) {
    mprintf("\tRegenerating dihedral parameters.\n");
    top.AssignDihedralParams( set0.DP() );
  }
  // Urey-Bradley
  updateCount = UpdateBondParams( set0.UB(), set1.UB() );
  if (updateCount > 0) {
    mprintf("\tRegenerating UB parameters.\n");
    top.AssignUBParams( set0.UB() );
  }
  // Improper parameters
  updateCount = 0;
  for (ParmHolder<DihedralParmType>::const_iterator newp = set1.IP().begin();
                                                    newp != set1.IP().end(); ++newp)
  {
    ParameterHolders::RetType ret = set0.IP().AddParm( newp->first, newp->second, true );
    if (ret != ParameterHolders::ERR) {
      if (ret == ParameterHolders::ADDED) {
        mprintf("\tAdded NEW improper parameter:\n");
        updateCount++;
      } else if (ret == ParameterHolders::UPDATED) {
        mprintf("\tUpdated improper parameter:\n");
        updateCount++;
      } else if (ret == ParameterHolders::SAME)
        mprintf("\tImproper parameter already present:\n");
        mprintf("\t\t%s %s %s %s %12.4f %12.4f\n", 
                *(newp->first[0]), *(newp->first[1]), *(newp->first[2]), *(newp->first[3]), newp->second.Pk(), newp->second.Phase());
    }
  }
  if (updateCount > 0) {
    mprintf("\tRegenerating dihedral parameters.\n");
    top.AssignImproperParams( set0.IP() );
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
  if (ds->Type() != DataSet::PARAMETERS) {
    mprinterr("Error: Set '%s' is not a parameter data set.\n", ds->legend());
    return CpptrajState::ERR;
  }
  DataSet_Parameters const& prm = static_cast<DataSet_Parameters const&>( *ds );
  Topology* dstop = State.DSL().GetTopology( argIn );
  if (dstop == 0) {
    mprinterr("Error: No topology specified.\n");
    return CpptrajState::ERR;
  }
  Topology& top = static_cast<Topology&>( *dstop );


  mprintf("\tUpdating parameters in topology '%s' using those in set '%s'\n",
          top.c_str(), prm.legend());

  UpdateParams(top, prm);

  return CpptrajState::OK;
}
