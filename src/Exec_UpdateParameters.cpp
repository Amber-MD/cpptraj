#include "Exec_UpdateParameters.h"
#include "CpptrajStdio.h"
#include "DataSet_Parameters.h"

// Exec_UpdateParameters::Help()
void Exec_UpdateParameters::Help() const
{
  mprintf("\t%s\n", DataSetList::TopArgs);
}

static inline void BondTypes(ParmHolder<int>& ParmIndices, Topology const& top, BondArray const& bonds) {
  for (BondArray::const_iterator b = bonds.begin(); b != bonds.end(); ++b)
  {
    AtomTypeHolder types(2);
    types.AddName( top[b->A1()].Type() );
    types.AddName( top[b->A2()].Type() );
    ParmIndices.AddParm( types, b->Idx(), false );
  }
}

static inline void AngleTypes(ParmHolder<int>& ParmIndices, Topology const& top, AngleArray const& angles) {
  for (AngleArray::const_iterator b = angles.begin(); b != angles.end(); ++b)
  {
    AtomTypeHolder types(3);
    types.AddName( top[b->A1()].Type() );
    types.AddName( top[b->A2()].Type() );
    types.AddName( top[b->A3()].Type() );
    ParmIndices.AddParm( types, b->Idx(), false );
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

  // Get list of existing parameters 
  ParmHolder<int> ParmIndices;
  // Bond parameters
  BondTypes(ParmIndices, top, top.Bonds());
  BondTypes(ParmIndices, top, top.BondsH());
  // Update bond parameters. We assume a parameter in top is never repeated.
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
  AngleTypes(ParmIndices, top, top.Angles());
  AngleTypes(ParmIndices, top, top.AnglesH());
  // Update bond parameters. We assume a parameter in top is never repeated.
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


  //for (ParmHolder<int>::const_iterator it = ParmIndices.begin(); it != ParmIndices.end(); ++it)
  //  mprintf("\t%s %s %i\n", *(it->first[0]), *(it->first[1]), it->second);
  return CpptrajState::OK;
}
