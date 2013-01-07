#include <cmath>
#include "ClusterDist.h"

ClusterMatrix ClusterDist_Num::PairwiseDist(int sieve) {
  int f2end = data_->Size();
  ClusterMatrix frameDistances( f2end );
  int f1end = f2end - sieve;
  for (int f1 = 0; f1 < f1end; f1 += sieve) {
    for (int f2 = f1 + sieve; f2 < f2end; f2 += sieve)
      frameDistances.SetElement( f1, f2, fabs(data_->Dval(f1) - data_->Dval(f2)) );
  }
  return frameDistances;
}

double ClusterDist_Num::CentroidDist(Centroid* c1, Centroid* c2) {
  return fabs(((Centroid_Num*)c1)->cval_ - ((Centroid_Num*)c2)->cval_);
}

double ClusterDist_Num::FrameCentroidDist(int f1, Centroid* c1) {
  return fabs(data_->Dval(f1) - ((Centroid_Num*)c1)->cval_);
}

void ClusterDist_Num::CalculateCentroid(Centroid* centIn, Cframes const& cframesIn) {
  Centroid_Num* cent = (Centroid_Num*)centIn;
  cent->cval_ = 0.0;
  for (Cframes_it frm = cframesIn.begin(); frm != cframesIn.end(); ++frm)
    cent->cval_ += data_->Dval( *frm );
  cent->cval_ /= (double)cframesIn.size();
}

Centroid* ClusterDist_Num::NewCentroid( Cframes const& cframesIn ) {
  Centroid_Num* cent = new Centroid_Num();
  CalculateCentroid( cent, cframesIn );
  return cent;
}

// -----------------------------------------------------------------------------
ClusterDist_DME::ClusterDist_DME(DataSet* dIn, std::string const& maskexpr) :
  coords_((DataSet_Coords*)dIn),
  mask_(maskexpr)
{
  coords_->Top().SetupIntegerMask( mask_ );
  frm1_.SetupFrameFromMask(mask_, coords_->Top().Atoms());
}

ClusterMatrix ClusterDist_DME::PairwiseDist(int sieve) {
  Frame frm2 = frm1_;
  int f2end = coords_->Size();
  ClusterMatrix frameDistances( f2end );
  int f1end = f2end - sieve;
  for (int f1 = 0; f1 < f1end; f1 += sieve) { 
    coords_->GetFrame( f1, frm1_, mask_ );
    for (int f2 = f1 + sieve; f2 < f2end; f2 += sieve) {
      coords_->GetFrame( f2, frm2,  mask_ );
      frameDistances.SetElement( f1, f2, frm1_.DISTRMSD( frm2 ) );
    }
  }
  return frameDistances;
}

double ClusterDist_DME::CentroidDist(Centroid* c1, Centroid* c2) {
  return ((Centroid_Coord*)c1)->cframe_.DISTRMSD( ((Centroid_Coord*)c2)->cframe_ );
}

double ClusterDist_DME::FrameCentroidDist(int f1, Centroid* c1) {
  coords_->GetFrame( f1, frm1_, mask_ );
  return frm1_.DISTRMSD( ((Centroid_Coord*)c1)->cframe_ );
}

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

// -----------------------------------------------------------------------------
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
  Frame frm2 = frm1_;
  int f2end = coords_->Size();
  ClusterMatrix frameDistances( f2end );
  int f1end = f2end - sieve;
  for (int f1 = 0; f1 < f1end; f1 += sieve) {
    coords_->GetFrame( f1, frm1_, mask_ );
    for (int f2 = f1 + sieve; f2 < f2end; f2 += sieve) {
      coords_->GetFrame( f2, frm2,  mask_ );
      if (nofit_) 
        rmsd = frm1_.RMSD_NoFit( frm2, useMass_ );
      else
        rmsd = frm1_.RMSD( frm2, useMass_ );
      frameDistances.SetElement( f1, f2, rmsd );
    }
  }
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
