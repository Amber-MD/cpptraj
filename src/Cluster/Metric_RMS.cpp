#include "Metric_RMS.h"
#include "Centroid_Coord.h"

/// CONSTRUCTOR
Cpptraj::Cluster::Metric_RMS::Metric_RMS(DataSet_Coords* dIn, AtomMask const& maskIn, 
                                 bool nofit, bool useMass) :
  coords_(dIn),
  mask_(maskIn),
  nofit_(nofit),
  useMass_(useMass)
{
  frm1_.SetupFrameFromMask(mask_, coords_->Top().Atoms());
  frm2_ = frm1_;
}

/** \return RMSD between two given frames. */
double Cpptraj::Cluster::Metric_RMS::FrameDist(int f1, int f2) {
  coords_->GetFrame( f1, frm1_, mask_ );
  coords_->GetFrame( f2, frm2_, mask_ );
  if (nofit_)
    return frm1_.RMSD_NoFit( frm2_, useMass_ );
  else
    return frm1_.RMSD( frm2_, useMass_ );
}

/** \return RMSD between two given centroids. */
double Cpptraj::Cluster::Metric_RMS::CentroidDist(Centroid* c1, Centroid* c2) {
  if (nofit_)
    return ((Centroid_Coord*)c1)->Cframe().RMSD_NoFit( ((Centroid_Coord*)c2)->Cframe(), useMass_ );
  else // Centroid is already at origin.
    return ((Centroid_Coord*)c1)->Cframe().RMSD_CenteredRef( ((Centroid_Coord*)c2)->Cframe(), 
                                                            useMass_ );
}

/** \return RMSD between given frame and centroid. */
double Cpptraj::Cluster::Metric_RMS::FrameCentroidDist(int f1, Centroid* c1) {
  coords_->GetFrame( f1, frm1_, mask_ );
  if (nofit_)
    return frm1_.RMSD_NoFit( ((Centroid_Coord*)c1)->Cframe(), useMass_ );
  else // Centroid is already at origin.
    return frm1_.RMSD_CenteredRef( ((Centroid_Coord*)c1)->Cframe(), useMass_ );
}

/** Compute the centroid (avg) coords for each atom from all frames in this
  * cluster. If fitting,  RMS fit to centroid as it is being built.
  */
void Cpptraj::Cluster::Metric_RMS::CalculateCentroid(Centroid* centIn,  Cframes const& cframesIn) {
  Matrix_3x3 Rot;
  Vec3 Trans;
  Centroid_Coord* cent = (Centroid_Coord*)centIn;
  // Reset atom count for centroid.
  cent->Cframe().ClearAtoms();
  for (Cframes_it frm = cframesIn.begin(); frm != cframesIn.end(); ++frm)
  {
    coords_->GetFrame( *frm, frm1_, mask_ );
    if (cent->Cframe().empty()) {
      cent->Cframe() = frm1_;
      if (!nofit_)
        cent->Cframe().CenterOnOrigin(useMass_);
    } else {
      if (!nofit_) {
        frm1_.RMSD_CenteredRef( cent->Cframe(), Rot, Trans, useMass_ );
        frm1_.Rotate( Rot );
      }
      cent->Cframe() += frm1_;
    }
  }
  cent->Cframe().Divide( (double)cframesIn.size() );
  //mprintf("\t\tFirst 3 centroid coords (of %i): %f %f %f\n", cent->Cframe().Natom(), 
  //        cent->cent->Cframe()[0], cent->Cframe()[1],cent->Cframe()[2]);
}

/** \return Average structure of given frames. */
Cpptraj::Cluster::Centroid* Cpptraj::Cluster::Metric_RMS::NewCentroid( Cframes const& cframesIn ) {
  // TODO: Incorporate mass?
  Centroid_Coord* cent = new Centroid_Coord( mask_.Nselected() );
  CalculateCentroid( cent, cframesIn );
  return cent;
}

// Subtract Notes
// FIXME: Handle single frame
// FIXME: Check if frame is in cluster?
void Cpptraj::Cluster::Metric_RMS::FrameOpCentroid(int frame, Centroid* centIn, double oldSize,
                                      CentOpType OP)
{
  Matrix_3x3 Rot;
  Vec3 Trans;
  Centroid_Coord* cent = (Centroid_Coord*)centIn;
  coords_->GetFrame( frame, frm1_, mask_ );
  if (!nofit_) {
    frm1_.RMSD_CenteredRef( cent->Cframe(), Rot, Trans, useMass_ );
    frm1_.Rotate( Rot );
  }
  cent->Cframe().Multiply( oldSize );
  if (OP == ADDFRAME) {
    cent->Cframe() += frm1_;
    cent->Cframe().Divide( oldSize + 1 );
  } else { // SUBTRACTFRAME
    cent->Cframe() -= frm1_;
    cent->Cframe().Divide( oldSize - 1 );
  }
}

/** \return Description of RMS calc. */
std::string Cpptraj::Cluster::Metric_RMS::Description() const {
  std::string description("rms " + mask_.MaskExpression());
  if (nofit_) description.append(" nofit");
  if (useMass_) description.append(" mass");
  return description;
}
