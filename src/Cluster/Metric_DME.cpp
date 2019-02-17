#include "Metric_DME.h"
#include "Centroid_Coord.h"
#include "../CpptrajStdio.h"

/** Initialize the metric. */
int Cpptraj::Cluster::Metric_DME::Init(DataSet_Coords* dIn, AtomMask const& maskIn)
{
  // TODO better error handles
  if (dIn == 0) {
    mprinterr("Internal Error: Metric_DME::Init() called with null data set.\n");
    return 1;
  }
  mprintf("DEBUG: Init DME metric for '%s', mask '%s'\n",
          dIn->legend(), maskIn.MaskString());
  coords_ = dIn;
  mask_ = maskIn;

  return 0;
}

/** Set up the metric. */
int Cpptraj::Cluster::Metric_DME::Setup() {
  if (coords_->Top().SetupIntegerMask( mask_ )) return 1;
  mprintf("DEBUG: DME metric topology: %s %s %i\n", coords_->legend(),
          coords_->Top().c_str(), coords_->Top().Natom());
  if (frm1_.SetupFrameFromMask(mask_, coords_->Top().Atoms())) return 1;
  frm2_ = frm1_;
  mprintf("DEBUG: Setup DME metric for %i atoms, %zu frames.\n", frm1_.Natom(), coords_->Size());
  return 0;
}

/** \return DME between two given frames. */
double Cpptraj::Cluster::Metric_DME::FrameDist(int f1, int f2) {
  coords_->GetFrame( f1, frm1_, mask_ );
  coords_->GetFrame( f2, frm2_, mask_ );
  return frm1_.DISTRMSD( frm2_ );
}

/** \return DME between two given centroids. */
double Cpptraj::Cluster::Metric_DME::CentroidDist(Centroid* c1, Centroid* c2) {
  return ((Centroid_Coord*)c1)->Cframe().DISTRMSD( ((Centroid_Coord*)c2)->Cframe() );
}

/** \return RMSD between given frame and centroid. */
double Cpptraj::Cluster::Metric_DME::FrameCentroidDist(int f1, Centroid* c1) {
  coords_->GetFrame( f1, frm1_, mask_ );
  return frm1_.DISTRMSD( ((Centroid_Coord*)c1)->Cframe() );
}

/** Compute the centroid (avg) coords for each atom from all frames in this
  * cluster. NOTE: For DME the centroid should probably be calculated via
  * internal coordinates; use RMS best-fit as a cheat.
  */
void Cpptraj::Cluster::Metric_DME::CalculateCentroid(Centroid* centIn,  Cframes const& cframesIn) {
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
      cent->Cframe().CenterOnOrigin(false);
    } else {
      frm1_.RMSD_CenteredRef( cent->Cframe(), Rot, Trans, false );
      frm1_.Rotate( Rot );
      cent->Cframe() += frm1_;
    }
  }
  cent->Cframe().Divide( (double)cframesIn.size() );
  //mprintf("\t\tFirst 3 centroid coords (of %i): %f %f %f\n", cent->Cframe().Natom(), 
  //        cent->cent->Cframe()[0], cent->Cframe()[1],cent->Cframe()[2]);
}

/** \return Average structure of given frames. */
Cpptraj::Cluster::Centroid* Cpptraj::Cluster::Metric_DME::NewCentroid( Cframes const& cframesIn ) {
  // TODO: Incorporate mass?
  Centroid_Coord* cent = new Centroid_Coord( mask_.Nselected() );
  CalculateCentroid( cent, cframesIn );
  return cent;
}

// Subtract Notes
// FIXME: Handle single frame
// FIXME: Check if frame is in cluster?
void Cpptraj::Cluster::Metric_DME::FrameOpCentroid(int frame, Centroid* centIn, double oldSize,
                                      CentOpType OP)
{
  Matrix_3x3 Rot;
  Vec3 Trans;
  Centroid_Coord* cent = (Centroid_Coord*)centIn;
  coords_->GetFrame( frame, frm1_, mask_ );
  frm1_.RMSD_CenteredRef( cent->Cframe(), Rot, Trans, false );
  frm1_.Rotate( Rot );
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
std::string Cpptraj::Cluster::Metric_DME::Description() const {
  return "dme " + mask_.MaskExpression();
}

void Cpptraj::Cluster::Metric_DME::Info() const {
  mprintf("\tMetric: DME");
  if (mask_.MaskExpression() == "*")
    mprintf(" (all atoms)\n");
  else
    mprintf(" (mask '%s')\n", mask_.MaskString());
}
