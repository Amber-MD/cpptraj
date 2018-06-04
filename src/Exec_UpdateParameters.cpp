#include "Exec_UpdateParameters.h"
#include "CpptrajStdio.h"
#include "DataSet_Parameters.h"

// Exec_UpdateParameters::Help()
void Exec_UpdateParameters::Help() const
{
  mprintf("\t%s setname <parm set>\n", DataSetList::TopArgs);
}

static inline void BondTypes(ParmHolder<int>& ParmIndices, Topology const& top, BondArray const& bonds, int& skipCount) {
  for (BondArray::const_iterator b = bonds.begin(); b != bonds.end(); ++b)
  {
    if (b->Idx() == -1)
      skipCount++;
    else {
      AtomTypeHolder types(2);
      types.AddName( top[b->A1()].Type() );
      types.AddName( top[b->A2()].Type() );
      ParmIndices.AddParm( types, b->Idx(), false );
    }
  }
}

static inline void AngleTypes(ParmHolder<int>& ParmIndices, Topology const& top, AngleArray const& angles, int& skipCount) {
  for (AngleArray::const_iterator b = angles.begin(); b != angles.end(); ++b)
  {
    if (b->Idx() == -1)
      skipCount++;
    else {
      AtomTypeHolder types(3);
      types.AddName( top[b->A1()].Type() );
      types.AddName( top[b->A2()].Type() );
      types.AddName( top[b->A3()].Type() );
      ParmIndices.AddParm( types, b->Idx(), false );
    }
  }
}

static inline void DihedralTypes(ParmHolder<int>& ParmIndices, Topology const& top, DihedralArray const& dih, int& skipCount) {
  for (DihedralArray::const_iterator b = dih.begin(); b != dih.end(); ++b)
  {
    if (b->Idx() == -1)
      skipCount++;
    else {
      AtomTypeHolder types(3);
      types.AddName( top[b->A1()].Type() );
      types.AddName( top[b->A2()].Type() );
      types.AddName( top[b->A3()].Type() );
      types.AddName( top[b->A4()].Type() );
      ParmIndices.AddParm( types, b->Idx(), false );
    }
  }
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
  // We assume a parameter in topology is never repeated, so as soon
  // as it is found we can move to the next one.
  // ParmIndices will map type names for existing bonds, angles, etc. to
  // existing type indices.
  ParmHolder<int> ParmIndices;
  // Bond parameters
  int skipCount = 0;
  BondTypes(ParmIndices, top, top.Bonds(), skipCount);
  BondTypes(ParmIndices, top, top.BondsH(), skipCount);
  if (skipCount > 0)
    mprintf("Warning: %i bonds were missing parameters and were skipped.\n", skipCount);
  for (ParmHolder<BondParmType>::const_iterator it1 = prm.BP().begin();
                                                it1 != prm.BP().end(); ++it1)
  {
    for (ParmHolder<int>::const_iterator it0 = ParmIndices.begin();
                                         it0 != ParmIndices.end(); ++it0)
    {
      if (it1->first == it0->first)
      {
        BondParmType& bp = top.SetBondParm(it0->second);
        mprintf("\tUpdating bond parameter %s - %s from %f %f to %f %f\n",
                *(it0->first[0]), *(it0->first[1]), bp.Rk(), bp.Req(),
                it1->second.Rk(), it1->second.Req());
        bp = it1->second;
        break;
      }
    }
  }
  ParmIndices.clear();
  // Angle parameters
  skipCount = 0;
  AngleTypes(ParmIndices, top, top.Angles(), skipCount);
  AngleTypes(ParmIndices, top, top.AnglesH(), skipCount);
  if (skipCount > 0)
    mprintf("Warning: %i angles were missing parameters and were skipped.\n", skipCount);
  for (ParmHolder<AngleParmType>::const_iterator it1 = prm.AP().begin();
                                                 it1 != prm.AP().end(); ++it1)
  {
    for (ParmHolder<int>::const_iterator it0 = ParmIndices.begin();
                                         it0 != ParmIndices.end(); ++it0)
    {
      if (it1->first == it0->first)
      {
        AngleParmType& bp = top.SetAngleParm(it0->second);
        mprintf("\tUpdating angle parameter %s - %s - %s from %f %f to %f %f\n",
                *(it0->first[0]), *(it0->first[1]), *(it0->first[2]), bp.Tk(), bp.Teq(),
                it1->second.Tk(), it1->second.Teq());
        bp = it1->second;
        break;
      }
    }
  }
  ParmIndices.clear();
  // Dihedral parameters
  skipCount = 0;
  DihedralTypes(ParmIndices, top, top.Dihedrals(), skipCount);
  DihedralTypes(ParmIndices, top, top.DihedralsH(), skipCount);
  if (skipCount > 0)
    mprintf("Warning: %i dihedrals were missing parameters and were skipped.\n", skipCount);
  for (ParmHolder<DihedralParmType>::const_iterator it1 = prm.DP().begin();
                                                    it1 != prm.DP().end(); ++it1)
  {
    for (ParmHolder<int>::const_iterator it0 = ParmIndices.begin();
                                         it0 != ParmIndices.end(); ++it0)
    {
      if (it1->first == it0->first)
      {
        DihedralParmType& bp = top.SetDihedralParm(it0->second);
        mprintf("\tUpdating dihedral parameter %s - %s - %s  %s from %f %f %f to %f %f %f\n",
                *(it0->first[0]), *(it0->first[1]), *(it0->first[2]), *(it0->first[3]),
                bp.Pk(), bp.Pn(), bp.Phase(),
                it1->second.Pk(), it1->second.Pn(), it1->second.Phase());
        bp = it1->second;
        break;
      }
    }
  }
  ParmIndices.clear();
  if (top.Chamber().HasChamber()) {
    ChamberParmType& chm = top.SetChamber();
    // UB parameters
    skipCount = 0;
    BondTypes(ParmIndices, top, chm.UB(), skipCount);
    if (skipCount > 0)
      mprintf("Warning: %i Urey-Bradley terms were missing parameters and were skipped.\n", skipCount);
    for (ParmHolder<BondParmType>::const_iterator it1 = prm.UB().begin();
                                                  it1 != prm.UB().end(); ++it1)
    {
      for (ParmHolder<int>::const_iterator it0 = ParmIndices.begin();
                                           it0 != ParmIndices.end(); ++it0)
      {
        if (it1->first == it0->first)
        {
          BondParmType& bp = chm.SetUBparm(it0->second);
          mprintf("\tUpdating UB parameter %s - %s from %f %f to %f %f\n",
                  *(it0->first[0]), *(it0->first[1]), bp.Rk(), bp.Req(),
                  it1->second.Rk(), it1->second.Req());
          bp = it1->second;
          break;
        }
      }
    }
    ParmIndices.clear();
    // Improper parameters
    skipCount = 0;
    DihedralTypes(ParmIndices, top, chm.Impropers(), skipCount);
    if (skipCount > 0)
      mprintf("Warning: %i impropers were missing parameters and were skipped.\n", skipCount);
    for (ParmHolder<DihedralParmType>::const_iterator it1 = prm.IP().begin();
                                                      it1 != prm.IP().end(); ++it1)
    {
      for (ParmHolder<int>::const_iterator it0 = ParmIndices.begin();
                                           it0 != ParmIndices.end(); ++it0)
      {
        if (it1->first == it0->first)
        {
          DihedralParmType& bp = chm.SetImproperParm(it0->second);
          mprintf("\tUpdating improper parameter %s - %s - %s  %s from %f %f %f to %f %f %f\n",
                  *(it0->first[0]), *(it0->first[1]), *(it0->first[2]), *(it0->first[3]),
                  bp.Pk(), bp.Pn(), bp.Phase(),
                  it1->second.Pk(), it1->second.Pn(), it1->second.Phase());
          bp = it1->second;
          break;
        }
      }
    }
    ParmIndices.clear();
  }


  //for (ParmHolder<int>::const_iterator it = ParmIndices.begin(); it != ParmIndices.end(); ++it)
  //  mprintf("\t%s %s %i\n", *(it->first[0]), *(it->first[1]), it->second);
  return CpptrajState::OK;
}
