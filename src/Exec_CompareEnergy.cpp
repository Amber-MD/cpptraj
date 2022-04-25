#include <cmath>
#include "Exec_CompareEnergy.h"
#include "CpptrajStdio.h"
#include "EnergyKernel_HarmonicBond.h"

/** CONSRUCTOR */
Exec_CompareEnergy::Exec_CompareEnergy() :
  Exec(GENERAL)
{
  SetHidden(true);
}

// Exec_CompareEnergy::Help()
void Exec_CompareEnergy::Help() const
{
  mprintf("\tcrd0 <set0> crd1 <set1>\n");
}

/** \return COORDINATES set corresponding to setname. */
DataSet_Coords* Exec_CompareEnergy::GetCoordsSet(DataSetList const& DSL,
                                                 std::string const& setname)
{
  if (setname.empty()) {
    mprinterr("Error: crdout: Specify COORDS dataset name.\n");
    return 0;
  }
  DataSet_Coords* CRD = (DataSet_Coords*)DSL.FindSetOfGroup( setname, DataSet::COORDINATES );
  if (CRD == 0) {
    mprinterr("Error: crdout: No COORDS set with name %s found.\n", setname.c_str());
    return 0;
  }
  mprintf("\tUsing set '%s'\n", CRD->legend());
  return CRD;
}

/// For calculating non-imaged distance^2 between two xyz points
template <typename T> T Distance2Kernel(T const* a1, T const* a2) {
  T x = a1[0] - a2[0];
  T y = a1[1] - a2[1];
  T z = a1[2] - a2[2];
  //double D = x*x + y*y + z*z;
  //fprintf(stdout,"Mask1=%8.3f %8.3f %8.3f Mask2=%8.3f %8.3f %8.3f D=%8.3f\n",
  //        a1[0],a1[1],a1[2],a2[0],a2[1],a2[2],D);
  return (x*x + y*y + z*z);
}
  
/// \return Energy for given bond
static inline double EBOND(Frame const& frame0,
                           BondType const& b0,
                           BondParmArray const& bpa0)
{
  if (b0.Idx() < 0) {
    mprintf("Warning: Bond %i -- %i has no parameters.\n", b0.A1()+1, b0.A2()+1);
    return 0;
  }
  BondParmType const& bp0 = bpa0[b0.Idx()];
  double r20 = Distance2Kernel<double>( frame0.XYZ(b0.A1()),
                                        frame0.XYZ(b0.A2()) );
  double ene = EnergyKernel_HarmonicBond<double>( sqrt(r20),
                                                  bp0.Rk(), bp0.Req() );
  return ene;
}

/** Compare bond energies between two frames. */
void Exec_CompareEnergy::CalcBondEnergy(Frame const& frame0,
                                        BondArray const& bonds0,
                                        BondParmArray const& bpa0,
                                        Frame const& frame1,
                                        BondArray const& bonds1,
                                        BondParmArray const& bpa1,
                                        double& E0, double& E1,
                                        Stats<double>& avgDelta2)
const
{
  if (bonds0.size() != bonds1.size()) {
    mprintf("Warning: Different # of bonds (%zu vs %zu)\n", bonds0.size(), bonds1.size());
    return;
  }
  for (unsigned int bidx = 0; bidx != bonds0.size(); bidx++) {
    if (bonds0[bidx].A1() != bonds1[bidx].A1() ||
        bonds0[bidx].A2() != bonds1[bidx].A2())
    {
      mprintf("Warning: Bond atom # mismatch (%i-%i vs %i-%i)\n",
              bonds0[bidx].A1()+1, bonds0[bidx].A2()+1, bonds1[bidx].A1()+1, bonds1[bidx].A2()+1);
      continue;
    }
    if ( (mask1_.AtomInCharMask(bonds0[bidx].A1()) && mask2_.AtomInCharMask(bonds0[bidx].A2())) &&
         (mask1_.AtomInCharMask(bonds0[bidx].A2()) && mask2_.AtomInCharMask(bonds0[bidx].A1())) )
    {

      double ene0 = EBOND(frame0, bonds0[bidx], bpa0);
      E0 += ene0;
      double ene1 = EBOND(frame1, bonds1[bidx], bpa1);
      E1 += ene1;
      double delta = ene0 - ene1;
      bondout_->Printf("\t%8i %8i %12.4f %12.4f %12.4f\n",
                       bonds0[bidx].A1()+1, bonds0[bidx].A2()+1,
                       ene0, ene1, delta);
      avgDelta2.accumulate( delta*delta );
    }
  }

}

/** Do bond energy comparison. */
void Exec_CompareEnergy::BondEnergy(Frame const& frame0, Topology const& top0,
                                    Frame const& frame1, Topology const& top1)
const
{
  Stats<double> avgDelta2;
  double E0 = 0;
  double E1 = 0;
  CalcBondEnergy(frame0, top0.Bonds(), top0.BondParm(),
                 frame1, top1.Bonds(), top1.BondParm(), E0, E1, avgDelta2);
  CalcBondEnergy(frame0, top0.BondsH(), top0.BondParm(),
                 frame1, top1.BondsH(), top1.BondParm(), E0, E1, avgDelta2);
  double rmse = sqrt( avgDelta2.mean() );
  bondout_->Printf("Bond E0   = %f\n", E0);
  bondout_->Printf("Bond E1   = %f\n", E1);
  bondout_->Printf("Bond RMSE = %f\n", rmse);
}
  
/** Compare energies between two coords sets. */
int Exec_CompareEnergy::GetEnergies(DataSet_Coords* crd0, DataSet_Coords* crd1)
const
{
  if (crd0->Size() < 1) {
    mprinterr("Error: '%s' has no frames.\n", crd0->legend());
    return 1;
  }
  if (crd1->Size() < 1) {
    mprinterr("Error: '%s' has no frames.\n", crd1->legend());
    return 1;
  }

  Frame frame0 = crd0->AllocateFrame();
  Frame frame1 = crd1->AllocateFrame();
  Topology const& top0 = crd0->Top();
  Topology const& top1 = crd1->Top();

  unsigned int maxframes = std::max( crd0->Size(), crd1->Size() );
  unsigned int idx0 = 0;
  unsigned int idx1 = 0;
  for (unsigned int idx = 0; idx != maxframes; idx++) {
    crd0->GetFrame( idx0++, frame0 );
    crd1->GetFrame( idx1++, frame1 );

    BondEnergy(frame0, top0, frame1, top1);

    // Reset counters if needed
    if (idx0 == crd0->Size())
      idx0 = 0;
    if (idx1 == crd1->Size())
      idx1 = 0;
  }

  return 0;
}

// Exec_CompareEnergy::Execute()
Exec::RetType Exec_CompareEnergy::Execute(CpptrajState& State, ArgList& argIn)
{
  DataSet_Coords* crd0 = GetCoordsSet(State.DSL(), argIn.GetStringKey("crd0"));
  if (crd0 == 0) return CpptrajState::ERR;
  DataSet_Coords* crd1 = GetCoordsSet(State.DSL(), argIn.GetStringKey("crd1"));
  if (crd1 == 0) return CpptrajState::ERR;

  bondout_ = State.DFL().AddCpptrajFile( argIn.GetStringKey("bondout"),
                                                     "bond comparison",
                                                     DataFileList::TEXT,
                                                     true );
  if (bondout_ == 0) {
    mprinterr("Internal Error: Could not allocate bond comparison file.\n");
    return CpptrajState::ERR;
  }

  mask1_.SetMaskString( argIn.GetStringKey("mask1") );
  mask2_.SetMaskString( argIn.GetStringKey("mask2") );

  mprintf("\tMask 1: %s\n", mask1_.MaskString());
  mprintf("\tMask 2: %s\n", mask2_.MaskString());

  if (crd0->Top().SetupCharMask( mask1_ )) {
    mprinterr("Error: Setting up mask '%s' failed.\n", mask1_.MaskString());
    return CpptrajState::ERR;
  }
  if (mask1_.None()) {
    mprinterr("Error: No atoms selected by '%s'\n", mask1_.MaskString());
    return CpptrajState::ERR;
  }
  if (crd0->Top().SetupCharMask( mask2_ )) {
    mprinterr("Error: Setting up mask '%s' failed.\n", mask2_.MaskString());
    return CpptrajState::ERR;
  }
  if (mask2_.None()) {
    mprinterr("Error: No atoms selected by '%s'\n", mask2_.MaskString());
    return CpptrajState::ERR;
  }

  mask1_.MaskInfo();
  mask2_.MaskInfo();

  if (GetEnergies(crd0, crd1)) return CpptrajState::ERR;

  return CpptrajState::OK;
}
