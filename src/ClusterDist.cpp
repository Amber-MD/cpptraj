#include <cmath>
#include "ClusterDist.h"
#include "Constants.h" // RADDEG, DEGRAD
#include "ProgressBar.h"
#ifdef _OPENMP
#  include "omp.h"
#endif
// TODO: All DataSet stuff const&
/// Calculate smallest difference between two angles (in degrees).
static double DistCalc_Dih(double d1, double d2) {
  double diff = fabs(d1 - d2);
  if (diff > 180.0)
    return (360.0 - diff);
  else
    return diff;
}

/// Calculate basic difference.
static double DistCalc_Std(double d1, double d2) {
  return fabs(d1 - d2);
}

static bool IsTorsionArray( DataSet const& dsIn ) {
  if (dsIn.ScalarMode() == DataSet::M_TORSION ||
      dsIn.ScalarMode() == DataSet::M_PUCKER  ||
      dsIn.ScalarMode() == DataSet::M_ANGLE     )
    return true;
  return false;
}

// -----------------------------------------------------------------------------
/** Calculate unambiguous average dihedral angle (in degrees) by converting to 
  * cartesian coords using x = cos(theta), y = sin(theta), and:
  *   tan(avgtheta) = avgy / avgx = SUM[sin(theta)] / SUM[cos(theta)]
  * See Eq. 2 from Altis et al., J. Chem. Phys., 126 p. 244111 (2007).
  */
static double AvgCalc_Dih( DataSet* dsIn, ClusterDist::Cframes const& cframesIn ) {
  double sumy = 0.0;
  double sumx = 0.0;
  // TODO: Convert angles to radians prior to this call?
  for (ClusterDist::Cframes_it frm = cframesIn.begin(); frm != cframesIn.end(); ++frm) {
    double theta = dsIn->Dval( *frm ) * DEGRAD;
    sumy += sin( theta );
    sumx += cos( theta );
  }
  return atan2(sumy, sumx) * RADDEG; 
}

static double AvgCalc_Std( DataSet* dsIn, ClusterDist::Cframes const& cframesIn ) {
  double val = 0.0;
  for (ClusterDist::Cframes_it frm = cframesIn.begin(); frm != cframesIn.end(); ++frm)
    val += dsIn->Dval( *frm );
  return (val / (double)cframesIn.size());
}

// ---------- Distance calc routines for single DataSet ------------------------
ClusterDist_Num::ClusterDist_Num( DataSet* dsIn ) :
  data_(dsIn)
{
  if (IsTorsionArray( *dsIn ))
    dcalc_ = DistCalc_Dih;
  else
    dcalc_ = DistCalc_Std;
}

ClusterMatrix ClusterDist_Num::PairwiseDist(int sieve) {
  int f1, f2;
  int f2end = data_->Size();
  ClusterMatrix frameDistances( f2end, sieve );
  int f1end = f2end - sieve;
#ifdef _OPENMP
#pragma omp parallel private(f1, f2)
{
#pragma omp for schedule(dynamic)
#endif
  for (f1 = 0; f1 < f1end; f1 += sieve) {
    for (f2 = f1 + sieve; f2 < f2end; f2 += sieve)
      frameDistances.SetElement( f1, f2, dcalc_(data_->Dval(f1), data_->Dval(f2)) );
  }
#ifdef _OPENMP
}
#endif
  return frameDistances;
}

double ClusterDist_Num::CentroidDist(Centroid* c1, Centroid* c2) {
  return dcalc_(((Centroid_Num*)c1)->cval_, ((Centroid_Num*)c2)->cval_);
}

double ClusterDist_Num::FrameCentroidDist(int f1, Centroid* c1) {
  return dcalc_(data_->Dval(f1), ((Centroid_Num*)c1)->cval_);
}

/** Calculate avg value of given frames. */
void ClusterDist_Num::CalculateCentroid(Centroid* centIn, Cframes const& cframesIn) {
  Centroid_Num* cent = (Centroid_Num*)centIn;
  if (IsTorsionArray(*data_))
    cent->cval_ = AvgCalc_Dih(data_, cframesIn);
  else
    cent->cval_ = AvgCalc_Std(data_, cframesIn);
}

/** \return A new centroid of the given frames. */
Centroid* ClusterDist_Num::NewCentroid( Cframes const& cframesIn ) {
  Centroid_Num* cent = new Centroid_Num();
  CalculateCentroid( cent, cframesIn );
  return cent;
}

// ---------- Distance calc routines for multiple DataSets (Euclid) ------------
ClusterDist_Euclid::ClusterDist_Euclid(DsArray const& dsIn) :
  dsets_(dsIn)
{
  for (DsArray::iterator ds = dsets_.begin(); ds != dsets_.end(); ++ds) {
    if ( IsTorsionArray( *(*ds) ) )
      dcalcs_.push_back( DistCalc_Dih );
    else
      dcalcs_.push_back( DistCalc_Std );
  }
}

ClusterMatrix ClusterDist_Euclid::PairwiseDist(int sieve) {
  int f1, f2;
  double dist, diff;
  DcArray::iterator dcalc;
  DsArray::iterator ds;
  int f2end = dsets_[0]->Size();
  ClusterMatrix frameDistances( f2end, sieve );
  int f1end = f2end - sieve;
#ifdef _OPENMP
#pragma omp parallel private(f1, f2, dist, diff, dcalc, ds)
{
#pragma omp for schedule(dynamic)
#endif
  for (f1 = 0; f1 < f1end; f1 += sieve) {
    for (f2 = f1 + sieve; f2 < f2end; f2 += sieve) {
      dist = 0.0;
      dcalc = dcalcs_.begin();
      for (ds = dsets_.begin(); ds != dsets_.end(); ++ds, ++dcalc) {
        diff = (*dcalc)((*ds)->Dval(f1), (*ds)->Dval(f2));
        dist += (diff * diff);
      }
      frameDistances.SetElement( f1, f2, sqrt(dist) );
    }
  }
#ifdef _OPENMP
}
#endif
  return frameDistances;
}

double ClusterDist_Euclid::CentroidDist(Centroid* c1, Centroid* c2) {
  double dist = 0.0;
  std::vector<double>::iterator c2val = ((Centroid_Multi*)c2)->cvals_.begin();
  DcArray::iterator dcalc = dcalcs_.begin();
  for (std::vector<double>::iterator c1val = ((Centroid_Multi*)c1)->cvals_.begin();
                                     c1val != ((Centroid_Multi*)c1)->cvals_.end();
                                     ++c1val, ++dcalc)
  {
    double diff = (*dcalc)(*c1val, *(c2val++));
    dist += (diff * diff);
  }
  return sqrt(dist);
}

double ClusterDist_Euclid::FrameCentroidDist(int f1, Centroid* c1) {
  double dist = 0.0;
  std::vector<double>::iterator c1val = ((Centroid_Multi*)c1)->cvals_.begin();
  DcArray::iterator dcalc = dcalcs_.begin();
  for (DsArray::iterator ds = dsets_.begin(); ds != dsets_.end(); ++ds) {
    double diff = (*dcalc)((*ds)->Dval(f1), *(c1val++));
    dist += (diff * diff);
  }
  return sqrt(dist);
}

void ClusterDist_Euclid::CalculateCentroid(Centroid* centIn, Cframes const& cframesIn) {
  Centroid_Multi* cent = (Centroid_Multi*)centIn;
  cent->cvals_.clear();
  for (DsArray::iterator ds = dsets_.begin(); ds != dsets_.end(); ++ds) {
    if (IsTorsionArray(*(*ds)))
      cent->cvals_.push_back( AvgCalc_Dih(*ds, cframesIn) );
    else
      cent->cvals_.push_back( AvgCalc_Std(*ds, cframesIn) );
  }
}

Centroid* ClusterDist_Euclid::NewCentroid(Cframes const& cframesIn) {
  Centroid_Multi* cent = new Centroid_Multi();
  CalculateCentroid(cent, cframesIn);
  return cent;
}

// ---------- Distance calc routines for COORDS DataSet using DME --------------
ClusterDist_DME::ClusterDist_DME(DataSet* dIn, std::string const& maskexpr) :
  coords_((DataSet_Coords*)dIn),
  mask_(maskexpr)
{
  coords_->Top().SetupIntegerMask( mask_ );
  frm1_.SetupFrameFromMask(mask_, coords_->Top().Atoms());
}

ClusterMatrix ClusterDist_DME::PairwiseDist(int sieve) {
  int f1, f2;
  Frame frm2 = frm1_;
  int f2end = coords_->Size();
  ClusterMatrix frameDistances( f2end, sieve );
  int f1end = f2end - sieve;
#ifdef _OPENMP
  Frame frm1 = frm1_;
# define frm1_ frm1
#pragma omp parallel private(f1, f2) firstprivate(frm1, frm2)
{
#pragma omp for schedule(dynamic)
#endif
  for (f1 = 0; f1 < f1end; f1 += sieve) { 
    coords_->GetFrame( f1, frm1_, mask_ );
    for (f2 = f1 + sieve; f2 < f2end; f2 += sieve) {
      coords_->GetFrame( f2, frm2,  mask_ );
      frameDistances.SetElement( f1, f2, frm1_.DISTRMSD( frm2 ) );
    }
  }
#ifdef _OPENMP
# undef frm1_
} // END pragma omp parallel
#endif
  return frameDistances;
}

double ClusterDist_DME::CentroidDist(Centroid* c1, Centroid* c2) {
  return ((Centroid_Coord*)c1)->cframe_.DISTRMSD( ((Centroid_Coord*)c2)->cframe_ );
}

double ClusterDist_DME::FrameCentroidDist(int f1, Centroid* c1) {
  coords_->GetFrame( f1, frm1_, mask_ );
  return frm1_.DISTRMSD( ((Centroid_Coord*)c1)->cframe_ );
}

/** Compute the centroid (avg) coords for each atom from all frames in this
  * cluster.
  */
void ClusterDist_DME::CalculateCentroid(Centroid* centIn, Cframes const& cframesIn) {
  Centroid_Coord* cent = (Centroid_Coord*)centIn;
  // Reset atom count for centroid.
  cent->cframe_.ClearAtoms();
  for (Cframes_it frm = cframesIn.begin(); frm != cframesIn.end(); ++frm)
  {
    coords_->GetFrame( *frm, frm1_, mask_ );
    if (cent->cframe_.empty())
      cent->cframe_ = frm1_;
    else
      cent->cframe_ += frm1_;
  }
  cent->cframe_.Divide( (double)cframesIn.size() );
  //mprintf("\t\tFirst 3 centroid coords (of %i): %f %f %f\n", cent->cframe_.Natom(), 
  //        cent->cent->cframe_[0], cent->cframe_[1],cent->cframe_[2]);
}

Centroid* ClusterDist_DME::NewCentroid( Cframes const& cframesIn ) {
  // TODO: Incorporate mass?
  Centroid_Coord* cent = new Centroid_Coord( mask_.Nselected() );
  CalculateCentroid( cent, cframesIn );
  return cent;
}

// ---------- Distance calc routines for COORDS DataSets using RMSD ------------
ClusterDist_RMS::ClusterDist_RMS(DataSet* dIn, std::string const& maskexpr, 
                                 bool nofit, bool useMass) :
  coords_((DataSet_Coords*)dIn),
  mask_(maskexpr),
  nofit_(nofit),
  useMass_(useMass)
{
  coords_->Top().SetupIntegerMask( mask_ );
  frm1_.SetupFrameFromMask(mask_, coords_->Top().Atoms());
}

ClusterMatrix ClusterDist_RMS::PairwiseDist(int sieve) {
  double rmsd;
  int f1, f2;
  Frame frm2 = frm1_;
  int f2end = coords_->Size();
  ClusterMatrix frameDistances( f2end, sieve );
  int f1end = f2end - sieve;
  ParallelProgress progress(f1end);
#ifdef _OPENMP
  Frame frm1 = frm1_;
# define frm1_ frm1
#pragma omp parallel private(f1, f2, rmsd) firstprivate(frm1, frm2, progress)
{
  progress.SetThread(omp_get_thread_num());
#pragma omp for schedule(dynamic)
#endif
  for (f1 = 0; f1 < f1end; f1 += sieve) {
    progress.Update(f1);
    coords_->GetFrame( f1, frm1_, mask_ );
    for (f2 = f1 + sieve; f2 < f2end; f2 += sieve) {
      coords_->GetFrame( f2, frm2,  mask_ );
      if (nofit_) 
        rmsd = frm1_.RMSD_NoFit( frm2, useMass_ );
      else
        rmsd = frm1_.RMSD( frm2, useMass_ );
      frameDistances.SetElement( f1, f2, rmsd );
    }
  }
#ifdef _OPENMP
# undef frm1_
} // END pragma omp parallel
#endif
  progress.Finish();
  return frameDistances;
}

double ClusterDist_RMS::CentroidDist(Centroid* c1, Centroid* c2) {
  if (nofit_)
    return ((Centroid_Coord*)c1)->cframe_.RMSD_NoFit( ((Centroid_Coord*)c2)->cframe_, useMass_ );
  else // Centroid is already at origin.
    return ((Centroid_Coord*)c1)->cframe_.RMSD_CenteredRef( ((Centroid_Coord*)c2)->cframe_, 
                                                            useMass_ );
}

double ClusterDist_RMS::FrameCentroidDist(int f1, Centroid* c1) {
  coords_->GetFrame( f1, frm1_, mask_ );
  if (nofit_)
    return frm1_.RMSD_NoFit( ((Centroid_Coord*)c1)->cframe_, useMass_ );
  else // Centroid is already at origin.
    return frm1_.RMSD_CenteredRef( ((Centroid_Coord*)c1)->cframe_, useMass_ );
}

/** Compute the centroid (avg) coords for each atom from all frames in this
  * cluster. If fitting,  RMS fit to centroid as it is being built.
  */
void ClusterDist_RMS::CalculateCentroid(Centroid* centIn,  Cframes const& cframesIn) {
  Matrix_3x3 Rot;
  Vec3 Trans;
  Centroid_Coord* cent = (Centroid_Coord*)centIn;
  // Reset atom count for centroid.
  cent->cframe_.ClearAtoms();
  for (Cframes_it frm = cframesIn.begin(); frm != cframesIn.end(); ++frm)
  {
    coords_->GetFrame( *frm, frm1_, mask_ );
    if (cent->cframe_.empty()) {
      cent->cframe_ = frm1_;
      if (!nofit_)
        cent->cframe_.CenterOnOrigin(useMass_);
    } else {
      if (!nofit_) {
        frm1_.RMSD_CenteredRef( cent->cframe_, Rot, Trans, useMass_ );
        frm1_.Rotate( Rot );
      }
      cent->cframe_ += frm1_;
    }
  }
  cent->cframe_.Divide( (double)cframesIn.size() );
  //mprintf("\t\tFirst 3 centroid coords (of %i): %f %f %f\n", cent->cframe_.Natom(), 
  //        cent->cent->cframe_[0], cent->cframe_[1],cent->cframe_[2]);
}

Centroid* ClusterDist_RMS::NewCentroid( Cframes const& cframesIn ) {
  // TODO: Incorporate mass?
  Centroid_Coord* cent = new Centroid_Coord( mask_.Nselected() );
  CalculateCentroid( cent, cframesIn );
  return cent;
}
