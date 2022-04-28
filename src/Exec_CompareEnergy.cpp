#include <cmath>
#include <cstdio> // sprintf
#include "Exec_CompareEnergy.h"
#include "CpptrajStdio.h"
#include "DataSet_double.h"
#include "EnergyKernel_Fourier.h"
#include "EnergyKernel_Harmonic.h"
#include "TorsionRoutines.h"

/** CONSRUCTOR */
Exec_CompareEnergy::Exec_CompareEnergy() :
  Exec(GENERAL),
  bondout_(0),
  bondDeltaE_(0),
  bondDeltaR_(0),
  angleout_(0),
  angleDeltaE_(0),
  angleDeltaR_(0)
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

// -----------------------------------------------------------------------------
/** Class to hold results from CalcEnergy */
class Eresults {
  public:
    Eresults(DataSet_double* deltaE, DataSet_double* deltaR) : E0_(0), E1_(0), deltaE_(deltaE), deltaR_(deltaR) {}

    void Print(CpptrajFile* outfile, const char* header) {
      double ermse = sqrt( avgEDelta2_.mean() );
      double rrmse = sqrt( avgRDelta2_.mean() );
      outfile->Printf("#%s E0       = %f\n", header, E0_);
      outfile->Printf("#%s E1       = %f\n", header, E1_);
      outfile->Printf("#%s <edelta> = %f\n", header, avgEDelta_.mean());
      outfile->Printf("#%s ene RMSE = %f\n", header, ermse);
      outfile->Printf("#%s <rdelta> = %f\n", header, avgRDelta_.mean());
      outfile->Printf("#%s len RMSE = %f\n", header, rrmse);
    }

    double E0_;
    double E1_;
    Stats<double> avgEDelta_;
    Stats<double> avgEDelta2_;
    Stats<double> avgRDelta_;
    Stats<double> avgRDelta2_;
    DataSet_double* deltaE_;
    DataSet_double* deltaR_;
};

/** Compare energies between two frames. */
template <typename T, typename P> void CalcEnergy( Eresults& result,
                                                   CpptrajFile* outfile,
                                                   Frame const& frame0,
                                                   std::vector<T> const& tarray0,
                                                   std::vector<P> const& parray0,
                                                   Frame const& frame1,
                                                   std::vector<T> const& tarray1,
                                                   std::vector<P> const& parray1,
                                                   std::vector<std::string> const& names,
                                                   double (*fxn)(Frame const&, T const&, P const&, double&) )
{
  for (unsigned int bidx = 0; bidx != tarray0.size(); bidx++) {
    T const& t0 = tarray0[bidx];
    T const& t1 = tarray1[bidx];

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
    outfile->Printf("%-s %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f\n",
                    names[bidx].c_str(),
                    ene0, ene1, edelta, r0, r1, rdelta);
    result.avgEDelta_.accumulate( edelta );
    result.avgEDelta2_.accumulate( edelta*edelta );
    result.avgRDelta_.accumulate( rdelta );
    result.avgRDelta2_.accumulate( rdelta*rdelta );
  }
}

// -----------------------------------------------------------------------------
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

/** Do bond energy comparison. */
void Exec_CompareEnergy::BondEnergy(Frame const& frame0, Topology const& top0,
                                    Frame const& frame1, Topology const& top1)
const
{
  bondout_->Printf("%-12s %-12s %12s %12s %12s %12s %12s %12s\n",
                   "#Name0", "Name1", "Ene0", "Ene1", "Edelta", "R0", "R1", "Rdelta");

  Eresults Ebond(bondDeltaE_, bondDeltaR_);
  CalcEnergy<BondType, BondParmType>(Ebond, bondout_,
                                     frame0, commonBonds0_, top0.BondParm(),
                                     frame1, commonBonds1_, top1.BondParm(),
                                     bondNames_, EBONDFXN);
  Ebond.Print( bondout_, "Bond" );
}
  
/** Set up array of selected bonds that top0 and top1 have in common. */
int Exec_CompareEnergy::SetupBondArray(Topology const& top0, BondArray const& bonds0, BondArray const& bonds1)
{
  char buffer[64];
  if (bonds0.size() != bonds1.size()) {
    mprintf("Warning: Different # of bonds (%zu vs %zu)\n", bonds0.size(), bonds1.size());
    return 1;
  }
  for (unsigned int bidx = 0; bidx != bonds0.size(); bidx++) {
    BondType const& b0 = bonds0[bidx];
    BondType const& b1 = bonds1[bidx];
    if ( (b0.A1() != b1.A1()) || (b0.A2() != b1.A2()) ) {
      mprintf("Warning: Bond atom # mismatch (%i-%i vs %i-%i)\n",
              b0.A1()+1, b0.A2()+1, b1.A1()+1, b1.A2()+1);
      continue;
    }
    if (b0.Idx() < 0) {
      mprintf("Warning: No parameters for bond 0\n");
      continue;
    }
    if (b1.Idx() < 0) {
      mprintf("Warning: No paramters for bond 1\n");
      continue;
    }

    if ( (mask1_.AtomInCharMask(b0.A1()) && mask2_.AtomInCharMask(b0.A2())) || 
         (mask1_.AtomInCharMask(b0.A2()) && mask2_.AtomInCharMask(b0.A1())) )
    {
      commonBonds0_.push_back( b0 );
      commonBonds1_.push_back( b1 );
      sprintf(buffer, "%-12s %-12s",
              top0.TruncResAtomName(b0.A1()).c_str(),
              top0.TruncResAtomName(b0.A2()).c_str());
      bondNames_.push_back( std::string(buffer) );
    }
  }

  return 0;
}

/** Set up arrays of selected bonds that top0 and top1 have in common. */
int Exec_CompareEnergy::SetupBondArrays(Topology const& top0, Topology const& top1) {
  commonBonds0_.clear();
  commonBonds1_.clear();
  bondNames_.clear();
  SetupBondArray(top0, top0.Bonds(), top1.Bonds());
  SetupBondArray(top0, top0.BondsH(), top1.BondsH());
  if (commonBonds0_.empty()) {
    mprinterr("Error: No bonds in common.\n");
    return 1;
  }
  return 0;
}

// -----------------------------------------------------------------------------
/// \return Energy for given angle 
static inline double EANGFXN(Frame const& frame0,
                             AngleType const& b0,
                             AngleParmType const& bp0,
                             double& r0)
{
  if (b0.Idx() < 0) {
    mprintf("Warning: Angle %i -- %i -- %i has no parameters.\n", b0.A1()+1, b0.A2()+1, b0.A3()+1);
    return 0;
  }
  double theta = CalcAngle( frame0.XYZ(b0.A1()),
                            frame0.XYZ(b0.A2()),
                            frame0.XYZ(b0.A3()) );
  double ene = EnergyKernel_Harmonic<double>( theta, bp0.Tk(), bp0.Teq() );
  return ene;
}

/** Do angle energy comparison. */
void Exec_CompareEnergy::AngleEnergy(Frame const& frame0, Topology const& top0,
                                    Frame const& frame1, Topology const& top1)
const
{
  angleout_->Printf("%-12s %-12s %-12s %12s %12s %12s %12s %12s %12s\n",
                   "#Name0", "Name1", "Name2", "Ene0", "Ene1", "Edelta", "R0", "R1", "Rdelta");

  Eresults Eangle(angleDeltaE_, angleDeltaR_);
  CalcEnergy<AngleType, AngleParmType>(Eangle, angleout_,
                                     frame0, commonAngles0_, top0.AngleParm(),
                                     frame1, commonAngles1_, top1.AngleParm(),
                                     angleNames_, EANGFXN);
  Eangle.Print( angleout_, "Angle" );
}
  
/** Set up array of selected angles that top0 and top1 have in common. */
int Exec_CompareEnergy::SetupAngleArray(Topology const& top0, AngleArray const& angles0, AngleArray const& angles1)
{
  char buffer[64];
  if (angles0.size() != angles1.size()) {
    mprintf("Warning: Different # of angles (%zu vs %zu)\n", angles0.size(), angles1.size());
    return 1;
  }
  for (unsigned int bidx = 0; bidx != angles0.size(); bidx++) {
    AngleType const& b0 = angles0[bidx];
    AngleType const& b1 = angles1[bidx];
    if ( (b0.A1() != b1.A1()) || (b0.A2() != b1.A2()) || (b0.A3() != b1.A3()) ) {
      mprintf("Warning: Angle atom # mismatch (%i-%i-%i vs %i-%i-%i)\n",
              b0.A1()+1, b0.A2()+1, b0.A3()+1, b1.A1()+1, b1.A2()+1, b1.A3()+1);
      continue;
    }
    if (b0.Idx() < 0) {
      mprintf("Warning: No parameters for angle 0\n");
      continue;
    }
    if (b1.Idx() < 0) {
      mprintf("Warning: No paramters for angle 1\n");
      continue;
    }

    if ( (mask1_.AtomInCharMask(b0.A1()) && mask2_.AtomInCharMask(b0.A2()) && mask3_.AtomInCharMask(b0.A3())) || 
         (mask1_.AtomInCharMask(b0.A3()) && mask2_.AtomInCharMask(b0.A2()) && mask3_.AtomInCharMask(b0.A1())) )
    {
      commonAngles0_.push_back( b0 );
      commonAngles1_.push_back( b1 );
      sprintf(buffer, "%-12s %-12s %-12s",
              top0.TruncResAtomName(b0.A1()).c_str(),
              top0.TruncResAtomName(b0.A2()).c_str(),
              top0.TruncResAtomName(b0.A3()).c_str());
      angleNames_.push_back( std::string(buffer) );
    }
  }

  return 0;
}

/** Set up arrays of selected angles that top0 and top1 have in common. */
int Exec_CompareEnergy::SetupAngleArrays(Topology const& top0, Topology const& top1) {
  commonAngles0_.clear();
  commonAngles1_.clear();
  angleNames_.clear();
  SetupAngleArray(top0, top0.Angles(), top1.Angles());
  SetupAngleArray(top0, top0.AnglesH(), top1.AnglesH());
  if (commonAngles0_.empty()) {
    mprinterr("Error: No angles in common.\n");
    return 1;
  }
  return 0;
}

// -----------------------------------------------------------------------------
/// \return Energy for given dihedral 
static inline double EDIHFXN(Frame const& frame0,
                             DihedralType const& b0,
                             DihedralParmType const& bp0,
                             double& r0)
{
  if (b0.Idx() < 0) {
    mprintf("Warning: Dihedral %i -- %i -- %i -- %i has no parameters.\n", b0.A1()+1, b0.A2()+1, b0.A3()+1, b0.A4()+1);
    return 0;
  }
  double theta = Torsion( frame0.XYZ(b0.A1()),
                          frame0.XYZ(b0.A2()),
                          frame0.XYZ(b0.A3()),
                          frame0.XYZ(b0.A4()) );
  double ene = EnergyKernel_Fourier<double>( theta, bp0.Pk(), bp0.Pn(), bp0.Phase() );
  return ene;
}

/** Do dihedral energy comparison. */
void Exec_CompareEnergy::DihedralEnergy(Frame const& frame0, Topology const& top0,
                                        Frame const& frame1, Topology const& top1)
const
{
  angleout_->Printf("%-12s %-12s %-12s %-12s %12s %12s %12s %12s %12s %12s\n",
                   "#Name0", "Name1", "Name2", "Name3", "Ene0", "Ene1", "Edelta", "R0", "R1", "Rdelta");

  Eresults Edihedral(dihedralDeltaE_, dihedralDeltaR_);
  CalcEnergy<DihedralType, DihedralParmType>(Edihedral, dihedralout_,
                                             frame0, commonDihedrals0_, top0.DihedralParm(),
                                             frame1, commonDihedrals1_, top1.DihedralParm(),
                                             dihedralNames_, EDIHFXN);
  Edihedral.Print( dihedralout_, "Dihedral" );
}
  
/** Set up array of selected dihedrals that top0 and top1 have in common. */
int Exec_CompareEnergy::SetupDihedralArray(Topology const& top0, DihedralArray const& dihedrals0, DihedralArray const& dihedrals1)
{
  char buffer[64];
  if (dihedrals0.size() != dihedrals1.size()) {
    mprintf("Warning: Different # of dihedrals (%zu vs %zu)\n", dihedrals0.size(), dihedrals1.size());
    return 1;
  }
  for (unsigned int bidx = 0; bidx != dihedrals0.size(); bidx++) {
    DihedralType const& b0 = dihedrals0[bidx];
    DihedralType const& b1 = dihedrals1[bidx];
    if ( (b0.A1() != b1.A1()) || (b0.A2() != b1.A2()) || (b0.A3() != b1.A3()) || (b0.A4() != b1.A4()) ) {
      mprintf("Warning: Dihedral atom # mismatch (%i-%i-%i-%i vs %i-%i-%i-%i)\n",
              b0.A1()+1, b0.A2()+1, b0.A3()+1, b0.A4()+1,
              b1.A1()+1, b1.A2()+1, b1.A3()+1, b1.A4()+1);
      continue;
    }
    if (b0.Idx() < 0) {
      mprintf("Warning: No parameters for dihedral 0\n");
      continue;
    }
    if (b1.Idx() < 0) {
      mprintf("Warning: No paramters for dihedral 1\n");
      continue;
    }

    if ( (mask1_.AtomInCharMask(b0.A1()) && mask2_.AtomInCharMask(b0.A2()) && mask3_.AtomInCharMask(b0.A3()) && mask4_.AtomInCharMask(b0.A4())) || 
         (mask1_.AtomInCharMask(b0.A4()) && mask2_.AtomInCharMask(b0.A3()) && mask3_.AtomInCharMask(b0.A2()) && mask4_.AtomInCharMask(b0.A1()))    )
    {
      commonDihedrals0_.push_back( b0 );
      commonDihedrals1_.push_back( b1 );
      sprintf(buffer, "%-12s %-12s %-12s %-12s",
              top0.TruncResAtomName(b0.A1()).c_str(),
              top0.TruncResAtomName(b0.A2()).c_str(),
              top0.TruncResAtomName(b0.A3()).c_str(),
              top0.TruncResAtomName(b0.A4()).c_str());
      dihedralNames_.push_back( std::string(buffer) );
    }
  }

  return 0;
}

/** Set up arrays of selected dihedrals that top0 and top1 have in common. */
int Exec_CompareEnergy::SetupDihedralArrays(Topology const& top0, Topology const& top1) {
  commonDihedrals0_.clear();
  commonDihedrals1_.clear();
  dihedralNames_.clear();
  SetupDihedralArray(top0, top0.Dihedrals(), top1.Dihedrals());
  SetupDihedralArray(top0, top0.DihedralsH(), top1.DihedralsH());
  if (commonDihedrals0_.empty()) {
    mprinterr("Error: No dihedrals in common.\n");
    return 1;
  }
  return 0;
}



// -----------------------------------------------------------------------------
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
    AngleEnergy(frame0, top0, frame1, top1);
    DihedralEnergy(frame0, top0, frame1, top1);

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
  angleout_ = State.DFL().AddCpptrajFile( argIn.GetStringKey("angleout"),
                                                     "angle comparison",
                                                     DataFileList::TEXT,
                                                     true );
  if (angleout_ == 0) {
    mprinterr("Internal Error: Could not allocate angle comparison file.\n");
    return CpptrajState::ERR;
  }
  dihedralout_ = State.DFL().AddCpptrajFile( argIn.GetStringKey("dihedralout"),
                                                     "dihedral comparison",
                                                     DataFileList::TEXT,
                                                     true );
  if (dihedralout_ == 0) {
    mprinterr("Internal Error: Could not allocate dihedral comparison file.\n");
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
  angleDeltaE_ = (DataSet_double*)State.DSL().AddSet(DataSet::DOUBLE, MetaData(dsname, "angleedelta"));
  if (angleDeltaE_ == 0) return CpptrajState::ERR;
  mprintf("\tAngle energy delta set: %s\n", angleDeltaE_->legend());
  angleDeltaR_ = (DataSet_double*)State.DSL().AddSet(DataSet::DOUBLE, MetaData(dsname, "anglerdelta"));
  if (angleDeltaR_ == 0) return CpptrajState::ERR;
  mprintf("\tAngle length delta set: %s\n", angleDeltaR_->legend());
  dihedralDeltaE_ = (DataSet_double*)State.DSL().AddSet(DataSet::DOUBLE, MetaData(dsname, "dihedraledelta"));
  if (dihedralDeltaE_ == 0) return CpptrajState::ERR;
  mprintf("\tDihedral energy delta set: %s\n", dihedralDeltaE_->legend());
  dihedralDeltaR_ = (DataSet_double*)State.DSL().AddSet(DataSet::DOUBLE, MetaData(dsname, "dihedralrdelta"));
  if (dihedralDeltaR_ == 0) return CpptrajState::ERR;
  mprintf("\tDihedral length delta set: %s\n", dihedralDeltaR_->legend());

  mask1_.SetMaskString( argIn.GetStringKey("mask1") );
  mask2_.SetMaskString( argIn.GetStringKey("mask2") );
  mask3_.SetMaskString( argIn.GetStringKey("mask3") );
  mask4_.SetMaskString( argIn.GetStringKey("mask4") );

  mprintf("\tMask 1: %s\n", mask1_.MaskString());
  mprintf("\tMask 2: %s\n", mask2_.MaskString());
  mprintf("\tMask 3: %s\n", mask3_.MaskString());
  mprintf("\tMask 4: %s\n", mask4_.MaskString());

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
  if (crd0->Top().SetupCharMask( mask3_ )) {
    mprinterr("Error: Setting up mask '%s' failed.\n", mask3_.MaskString());
    return CpptrajState::ERR;
  }
  if (mask3_.None()) {
    mprinterr("Error: No atoms selected by '%s'\n", mask3_.MaskString());
    return CpptrajState::ERR;
  }
  if (crd0->Top().SetupCharMask( mask4_ )) {
    mprinterr("Error: Setting up mask '%s' failed.\n", mask4_.MaskString());
    return CpptrajState::ERR;
  }
  if (mask4_.None()) {
    mprinterr("Error: No atoms selected by '%s'\n", mask4_.MaskString());
    return CpptrajState::ERR;
  }

  mask1_.MaskInfo();
  mask2_.MaskInfo();
  mask3_.MaskInfo();
  mask4_.MaskInfo();

  if (SetupBondArrays(crd0->Top(), crd1->Top())) return CpptrajState::ERR;
  if (SetupAngleArrays(crd0->Top(), crd1->Top())) return CpptrajState::ERR;
  if (SetupDihedralArrays(crd0->Top(), crd1->Top())) return CpptrajState::ERR;

  if (GetEnergies(crd0, crd1)) return CpptrajState::ERR;

  return CpptrajState::OK;
}
