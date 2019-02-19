#include "Metric_SRMSD.h"
#include "Centroid_Coord.h"
#include "../CpptrajStdio.h"

int Cpptraj::Cluster::Metric_SRMSD::Init(DataSet_Coords* dIn, AtomMask const& maskIn, 
                                     bool nofit, bool useMass, int debugIn)
{
  // TODO better error handles
  if (dIn == 0) {
    mprinterr("Internal Error: Metric_SRMSD::Init() called with null data set.\n");
    return 1;
  }
  coords_ = dIn;
  mask_ = maskIn;
  SRMSD_.InitSymmRMSD(!nofit, useMass, debugIn);
  return 0;
}

int Cpptraj::Cluster::Metric_SRMSD::Setup() {
  if (coords_->Top().SetupIntegerMask( mask_ )) return 1;
  mprintf("DEBUG: SRMSD metric topology: %s %s %i\n", coords_->legend(),
          coords_->Top().c_str(), coords_->Top().Natom());
  // false = no remap warning
  if (SRMSD_.SetupSymmRMSD(coords_->Top(), mask_, false)) return 1;
  frm1_.SetupFrameFromMask(mask_, coords_->Top().Atoms());
  frm2_ = frm1_;
  return 0;
}

double Cpptraj::Cluster::Metric_SRMSD::FrameDist(int f1, int f2) {
  coords_->GetFrame( f1, frm1_, mask_ );
  coords_->GetFrame( f2, frm2_, mask_ );
  return SRMSD_.SymmRMSD(frm1_, frm2_);
}

double Cpptraj::Cluster::Metric_SRMSD::CentroidDist(Centroid* c1, Centroid* c2) {
  // Centroid is already at origin.
  return SRMSD_.SymmRMSD_CenteredRef( ((Centroid_Coord*)c1)->Cframe(),
                                      ((Centroid_Coord*)c2)->Cframe() );
}

double Cpptraj::Cluster::Metric_SRMSD::FrameCentroidDist(int f1, Centroid* c1) {
  coords_->GetFrame( f1, frm1_, mask_ );
  // Centroid is already at origin.
  return SRMSD_.SymmRMSD_CenteredRef( frm1_, ((Centroid_Coord*)c1)->Cframe() );
}

/** Compute the centroid (avg) coords for each atom from all frames in this
  * cluster. If fitting,  RMS fit to centroid as it is being built.
  */
void Cpptraj::Cluster::Metric_SRMSD::CalculateCentroid(Centroid* centIn,  Cframes const& cframesIn) {
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
      if (SRMSD_.Fit())
        cent->Cframe().CenterOnOrigin(SRMSD_.UseMass());
    } else {
      SRMSD_.SymmRMSD_CenteredRef( frm1_, cent->Cframe() );
      // Remap atoms
      frm2_.SetCoordinatesByMap( frm1_, SRMSD_.AMap() );
      if (SRMSD_.Fit()) {
        frm2_.Translate( SRMSD_.TgtTrans() );
        frm2_.Rotate( SRMSD_.RotMatrix() );
      }
      cent->Cframe() += frm2_;
    }
  }
  cent->Cframe().Divide( (double)cframesIn.size() );
  //mprintf("\t\tFirst 3 centroid coords (of %i): %f %f %f\n", cent->cframe_.Natom(), 
  //        cent->cent->cframe_[0], cent->cframe_[1],cent->cframe_[2]);
}

Cpptraj::Cluster::Centroid* Cpptraj::Cluster::Metric_SRMSD::NewCentroid( Cframes const& cframesIn ) {
  // TODO: Incorporate mass?
  Centroid_Coord* cent = new Centroid_Coord( mask_.Nselected() );
  CalculateCentroid( cent, cframesIn );
  return cent;
}

void Cpptraj::Cluster::Metric_SRMSD::FrameOpCentroid(int frame, Centroid* centIn, double oldSize,
                                        CentOpType OP)
{
  Matrix_3x3 Rot;
  Vec3 Trans;
  Centroid_Coord* cent = (Centroid_Coord*)centIn;
  coords_->GetFrame( frame, frm1_, mask_ );
  SRMSD_.SymmRMSD_CenteredRef( frm1_, cent->Cframe() );
  // Remap atoms
  frm2_.SetCoordinatesByMap( frm1_, SRMSD_.AMap() );
  if (SRMSD_.Fit()) {
    frm2_.Translate( SRMSD_.TgtTrans() );
    frm2_.Rotate( SRMSD_.RotMatrix() );
  }
  cent->Cframe().Multiply( oldSize );
  if (OP == ADDFRAME) {
    cent->Cframe() += frm2_;
    cent->Cframe().Divide( oldSize + 1 );
  } else { // SUBTRACTFRAME
    cent->Cframe() -= frm2_;
    cent->Cframe().Divide( oldSize - 1 );
  }
}

/** \return Description of symmetric RMS calc. */
std::string Cpptraj::Cluster::Metric_SRMSD::Description() const {
  std::string description("srmsd " + mask_.MaskExpression());
  if (!SRMSD_.Fit()) description.append(" nofit");
  if (SRMSD_.UseMass()) description.append(" mass");
  return description;
}

void Cpptraj::Cluster::Metric_SRMSD::Info() const {
  mprintf("\tMetric: Symmetric RMSD");
  if (mask_.MaskExpression() == "*")
    mprintf(" (all atoms)");
  else
    mprintf(" (mask '%s')", mask_.MaskString());
  if (SRMSD_.UseMass())
    mprintf(", mass-weighted");
  if (!SRMSD_.Fit())
    mprintf(", no fitting");
  else
    mprintf(" best-fit");
  mprintf("\n");
}
