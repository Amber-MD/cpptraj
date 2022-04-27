#include <cmath>
#include "Exec_CompareEnergy.h"
#include "CpptrajStdio.h"
#include "DataSet_double.h"
#include "EnergyKernel_Harmonic.h"

/** CONSRUCTOR */
Exec_CompareEnergy::Exec_CompareEnergy() :
  Exec(GENERAL),
  bondout_(0),
  bondDeltaE_(0),
  bondDeltaR_(0)
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
                           BondParmArray const& bpa0,
                           double& r0)
{
  if (b0.Idx() < 0) {
    mprintf("Warning: Bond %i -- %i has no parameters.\n", b0.A1()+1, b0.A2()+1);
    return 0;
  }
  BondParmType const& bp0 = bpa0[b0.Idx()];
  double r20 = Distance2Kernel<double>( frame0.XYZ(b0.A1()),
                                        frame0.XYZ(b0.A2()) );
  r0 = sqrt(r20);
  double ene = EnergyKernel_Harmonic<double>( r0, bp0.Rk(), bp0.Req() );
  return ene;
}

/// \return Energy for given bond
static inline double EBONDFXN(Frame const& frame0,
                              BondType const& b0,
                              BondParmType const& bp0,
                              double& r0)
{
  if (b0.Idx() < 0) {
    mprintf("Warning: Bond %i -- %i has no parameters.\n", b0.A1()+1, b0.A2()+1);
    return 0;
  }
  double r20 = Distance2Kernel<double>( frame0.XYZ(b0.A1()),
                                        frame0.XYZ(b0.A2()) );
  r0 = sqrt(r20);
  double ene = EnergyKernel_Harmonic<double>( r0, bp0.Rk(), bp0.Req() );
  return ene;
}

class Eresults {
  public:
    Eresults(DataSet_double* deltaE, DataSet_double* deltaR) : E0_(0), E1_(0), deltaE_(deltaE), deltaR_(deltaR) {}

    double E0_;
    double E1_;
    Stats<double> avgEDelta_;
    Stats<double> avgEDelta2_;
    Stats<double> avgRDelta_;
    Stats<double> avgRDelta2_;
    DataSet_double* deltaE_;
    DataSet_double* deltaR_;
};


template <typename T, typename P> void CalcEnergy( Eresults& result,
                                                   Topology const& top0,
                                                   Frame const& frame0,
                                                   std::vector<T> const& tarray0,
                                                   std::vector<P> const& parray0,
                                                   Topology const& top1,
                                                   Frame const& frame1,
                                                   std::vector<T> const& tarray1,
                                                   std::vector<P> const& parray1,
                                                   CharMask const& mask1, CharMask const& mask2,
                                                   double (*fxn)(Frame const&, T const&, P const&, double&) )
{
  if (tarray0.size() != tarray1.size()) {
    mprintf("Warning: Different # of bonds (%zu vs %zu)\n", tarray0.size(), tarray1.size());
    return;
  }
  for (unsigned int bidx = 0; bidx != tarray0.size(); bidx++) {
    T const& t0 = tarray0[bidx];
    T const& t1 = tarray1[bidx];
    if ( (t0.A1() != t1.A1()) || (t0.A2() != t1.A2()) ) { // FIXME
      mprintf("Warning: Atom # mismatch.\n");
      //mprintf("Warning: Bond atom # mismatch (%i-%i vs %i-%i)\n",
      //        tarray0[bidx].A1()+1, tarray0[bidx].A2()+1, tarray1[bidx].A1()+1, tarray1[bidx].A2()+1);
      continue;
    }
    if (t0.Idx() < 0) {
      mprintf("Warning: No parameters for 0\n");
      continue;
    }
    if (t1.Idx() < 0) {
      mprintf("Warning: No paramters for 1\n");
      continue;
    }
    if ( (mask1.AtomInCharMask(t0.A1()) && mask2.AtomInCharMask(t0.A2())) || 
         (mask1.AtomInCharMask(t0.A2()) && mask2.AtomInCharMask(t0.A1())) )
    {
      double r0 = 0;
      double r1 = 0;
      double ene0 = fxn(frame0, tarray0[bidx], parray0[t0.Idx()], r0);
      result.E0_ += ene0;
      double ene1 = fxn(frame1, tarray1[bidx], parray1[t1.Idx()], r1);
      result.E1_ += ene1;
      double edelta = ene1 - ene0;
      double rdelta = r1 - r0;
      result.deltaE_->AddElement( edelta );
      result.deltaR_->AddElement( rdelta );
      //bondout_->Printf("%-12s %-12s %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f\n",
      //                 top0.TruncResAtomName(bonds0[bidx].A1()).c_str(),
      //                 top0.TruncResAtomName(bonds0[bidx].A2()).c_str(),
      //                 ene0, ene1, edelta, r0, r1, rdelta);
      //avgEDelta.accumulate( edelta );
      //avgEDelta2.accumulate( edelta*edelta );
      //avgRDelta.accumulate( rdelta );
      //avgRDelta2.accumulate( rdelta*rdelta );
    }
  }
}


/** Compare bond energies between two frames. */
void Exec_CompareEnergy::CalcBondEnergy(Topology const& top0,
                                        Frame const& frame0,
                                        BondArray const& bonds0,
                                        BondParmArray const& bpa0,
                                        Topology const& top1,
                                        Frame const& frame1,
                                        BondArray const& bonds1,
                                        BondParmArray const& bpa1,
                                        double& E0, double& E1,
                                        Stats<double>& avgEDelta,
                                        Stats<double>& avgEDelta2,
                                        Stats<double>& avgRDelta,
                                        Stats<double>& avgRDelta2)
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
    if ( (mask1_.AtomInCharMask(bonds0[bidx].A1()) && mask2_.AtomInCharMask(bonds0[bidx].A2())) || 
         (mask1_.AtomInCharMask(bonds0[bidx].A2()) && mask2_.AtomInCharMask(bonds0[bidx].A1())) )
    {
      double r0 = 0;
      double r1 = 0;
      double ene0 = EBOND(frame0, bonds0[bidx], bpa0, r0);
      E0 += ene0;
      double ene1 = EBOND(frame1, bonds1[bidx], bpa1, r1);
      E1 += ene1;
      double edelta = ene1 - ene0;
      double rdelta = r1 - r0;
      bondDeltaE_->AddElement( edelta );
      bondDeltaR_->AddElement( rdelta );
      bondout_->Printf("%-12s %-12s %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f\n",
                       top0.TruncResAtomName(bonds0[bidx].A1()).c_str(),
                       top0.TruncResAtomName(bonds0[bidx].A2()).c_str(),
                       ene0, ene1, edelta, r0, r1, rdelta);
      avgEDelta.accumulate( edelta );
      avgEDelta2.accumulate( edelta*edelta );
      avgRDelta.accumulate( rdelta );
      avgRDelta2.accumulate( rdelta*rdelta );
    }
  }

}

/** Do bond energy comparison. */
void Exec_CompareEnergy::BondEnergy(Frame const& frame0, Topology const& top0,
                                    Frame const& frame1, Topology const& top1)
const
{
  bondout_->Printf("%-12s %-12s %12s %12s %12s %12s %12s %12s\n",
                   "#Name0", "Name1", "Ene0", "Ene1", "Edelta", "R0", "R1", "Rdelta");

  // DEBUG
  Eresults Ebond(bondDeltaE_, bondDeltaR_);
  CalcEnergy<BondType, BondParmType>(Ebond,
                                     top0, frame0, top0.Bonds(), top0.BondParm(),
                                     top1, frame1, top1.Bonds(), top1.BondParm(),
                                     mask1_, mask2_, EBONDFXN);
  // DEBUG

  Stats<double> avgEDelta, avgEDelta2, avgRDelta, avgRDelta2;
  double E0 = 0;
  double E1 = 0;
  CalcBondEnergy(top0, frame0, top0.Bonds(), top0.BondParm(),
                 top1, frame1, top1.Bonds(), top1.BondParm(), E0, E1,
                 avgEDelta, avgEDelta2, avgRDelta, avgRDelta2);
  CalcBondEnergy(top0, frame0, top0.BondsH(), top0.BondParm(),
                 top1, frame1, top1.BondsH(), top1.BondParm(), E0, E1,
                 avgEDelta, avgEDelta2, avgRDelta, avgRDelta2);
  double ermse = sqrt( avgEDelta2.mean() );
  double rrmse = sqrt( avgRDelta2.mean() );
  bondout_->Printf("#Bond E0       = %f\n", E0);
  bondout_->Printf("#Bond E1       = %f\n", E1);
  bondout_->Printf("#Bond <edelta> = %f\n", avgEDelta.mean());
  bondout_->Printf("#Bond ene RMSE = %f\n", ermse);
  bondout_->Printf("#Bond <rdelta> = %f\n", avgRDelta.mean());
  bondout_->Printf("#Bond len RMSE = %f\n", rrmse);
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

  std::string dsname = argIn.GetStringKey("name");
  if (dsname.empty())
    dsname = State.DSL().GenerateDefaultName("ECOMPARE");
  bondDeltaE_ = (DataSet_double*)State.DSL().AddSet(DataSet::DOUBLE, MetaData(dsname, "bondedelta"));
  if (bondDeltaE_ == 0) return CpptrajState::ERR;
  mprintf("\tBond energy delta set: %s\n", bondDeltaE_->legend());
  bondDeltaR_ = (DataSet_double*)State.DSL().AddSet(DataSet::DOUBLE, MetaData(dsname, "bondrdelta"));
  if (bondDeltaR_ == 0) return CpptrajState::ERR;
  mprintf("\tBond length delta set: %s\n", bondDeltaR_->legend());

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
