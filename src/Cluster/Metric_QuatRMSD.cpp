#include "Metric_QuatRMSD.h"
#include "Centroid_Coord.h"
#include "Cframes.h"
#include "../CpptrajStdio.h"
#include "../QuaternionRMSD.h"

/** Initialize the metric. */
int Cpptraj::Cluster::Metric_QuatRMSD::Init(DataSet_Coords* dIn, AtomMask const& maskIn, 
                                            bool useMass)
{
  // TODO better error handles
  if (dIn == 0) {
    mprinterr("Internal Error: Metric_QuatRMSD::Init() called with null data set.\n");
    return 1;
  }
# ifdef DEBUG_CLUSTER
  mprintf("DEBUG: Init QRMSD metric for '%s', mask '%s', nofit=%i, usemass=%i\n",
          dIn->legend(), maskIn.MaskString(), 0, (int)useMass);
# endif
  coords_ = dIn;
  mask_ = maskIn;
//  nofit_ = nofit;
  useMass_ = useMass;

  return 0;
}

/** Set up the metric. */
int Cpptraj::Cluster::Metric_QuatRMSD::Setup() {
  if (coords_->Top().SetupIntegerMask( mask_ )) return 1;
  mask_.MaskInfo();
# ifdef DEBUG_CLUSTER
  mprintf("DEBUG: QRMSD metric topology: %s %s %i\n", coords_->legend(),
          coords_->Top().c_str(), coords_->Top().Natom());
# endif
  if (frm1_.SetupFrameFromMask(mask_, coords_->Top().Atoms())) return 1;
  frm2_ = frm1_;
# ifdef DEBUG_CLUSTER
  mprintf("DEBUG: Setup QRMSD metric for %i atoms, %zu frames.\n", frm1_.Natom(), coords_->Size());
# endif
  return 0;
}

/** \return RMSD between two given frames. */
double Cpptraj::Cluster::Metric_QuatRMSD::FrameDist(int f1, int f2) {
  coords_->GetFrame( f1, frm1_, mask_ );
  coords_->GetFrame( f2, frm2_, mask_ );
# ifdef DEBUG_CLUSTER
  double rms;
  frm1_.printAtomCoord(0);
  frm2_.printAtomCoord(0);
  //if (nofit_)
  //  rms = frm1_.RMSD_NoFit( frm2_, useMass_ );
  //else
    frm2_.CenterOnOrigin(useMass_);
    rms = QuaternionRMSD_CenteredRef( frm2_, frm1_, useMass_ );
  mprintf("\tMetric_QuatRMSD::FrameDist(%i, %i)= %g\n", f1, f2, rms);
  return rms;
# else
  //if (nofit_)
  //  return frm1_.RMSD_NoFit( frm2_, useMass_ );
  //else
    frm2_.CenterOnOrigin(useMass_);
    return QuaternionRMSD_CenteredRef( frm2_, frm1_, useMass_ );
# endif
}

/** \return RMSD between two given centroids. */
double Cpptraj::Cluster::Metric_QuatRMSD::CentroidDist(Centroid* c1, Centroid* c2) {
  //if (nofit_)
  //  return ((Centroid_Coord*)c1)->Cframe().RMSD_NoFit( ((Centroid_Coord*)c2)->Cframe(), useMass_ );
  //else // Centroid is already at origin.
    return QuaternionRMSD_CenteredRef( ((Centroid_Coord*)c2)->Cframe(),
                                       ((Centroid_Coord*)c1)->Cframe(), useMass_ );
}

/** \return RMSD between given frame and centroid. */
double Cpptraj::Cluster::Metric_QuatRMSD::FrameCentroidDist(int f1, Centroid* c1) {
  coords_->GetFrame( f1, frm1_, mask_ );
  //if (nofit_)
  //  return frm1_.RMSD_NoFit( ((Centroid_Coord*)c1)->Cframe(), useMass_ );
  //else // Centroid is already at origin.
    return QuaternionRMSD_CenteredRef( ((Centroid_Coord*)c1)->Cframe(), frm1_, useMass_ );
}

/** Compute the centroid (avg) coords for each atom from all frames in this
  * cluster. If fitting,  RMS fit to centroid as it is being built.
  */
void Cpptraj::Cluster::Metric_QuatRMSD::CalculateCentroid(Centroid* centIn,  Cframes const& cframesIn) {
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
      //if (!nofit_)
        cent->Cframe().CenterOnOrigin(useMass_);
    } else {
      //if (!nofit_) {
        QuaternionRMSD_CenteredRef( cent->Cframe(), frm1_, Rot, Trans, useMass_ );
        //Rot.Print("CalculateCentroid"); // DEBUG
        frm1_.InverseRotate( Rot );
      //}
      cent->Cframe() += frm1_;
    }
  }
  //mprintf("DEBUG: Metric_QuatRMSD::CalculateCentroid divide by %zu\n", cframesIn.size()); 
  cent->Cframe().Divide( (double)cframesIn.size() );
  //mprintf("\t\tFirst 3 centroid coords (of %i): %f %f %f\n", cent->Cframe().Natom(), 
  //        cent->cent->Cframe()[0], cent->Cframe()[1],cent->Cframe()[2]);
}

/** \return Average structure of given frames. */
Cpptraj::Cluster::Centroid* Cpptraj::Cluster::Metric_QuatRMSD::NewCentroid( Cframes const& cframesIn ) {
  // TODO: Incorporate mass?
  Centroid_Coord* cent = new Centroid_Coord( mask_.Nselected() );
  CalculateCentroid( cent, cframesIn );
  return cent;
}

// Subtract Notes
// TODO: Handle single frame
// TODO: Check if frame is in cluster?
void Cpptraj::Cluster::Metric_QuatRMSD::FrameOpCentroid(int frame, Centroid* centIn, double oldSize,
                                      CentOpType OP)
{
  Matrix_3x3 Rot;
  Vec3 Trans;
  Centroid_Coord* cent = (Centroid_Coord*)centIn;
  coords_->GetFrame( frame, frm1_, mask_ );
  //if (!nofit_) {
    QuaternionRMSD_CenteredRef( cent->Cframe(), frm1_, Rot, Trans, useMass_ );
    frm1_.InverseRotate( Rot );
  //}
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
std::string Cpptraj::Cluster::Metric_QuatRMSD::Description() const {
  std::string description("qrmsd " + mask_.MaskExpression());
  //if (nofit_) description.append(" nofit");
  if (useMass_) description.append(" mass");
  return description;
}

void Cpptraj::Cluster::Metric_QuatRMSD::Info() const {
  mprintf("\tMetric: Quaternion RMSD");
  if (mask_.MaskExpression() == "*")
    mprintf(" (all atoms)");
  else
    mprintf(" (mask '%s')", mask_.MaskString());
  if (useMass_)
    mprintf(", mass-weighted");
  //if (nofit_)
  //  mprintf(", no fitting");
  //else
    mprintf(" best-fit");
  mprintf("\n");
}
