#include "Metric_QuatRMSD.h"
#include "Centroid_Coord.h"
#include "Cframes.h"
#include "../CpptrajStdio.h"
#include "../QuaternionRMSD.h"

int Cpptraj::Cluster::Metric_QuatRMSD::Init(DataSet_Coords* dIn, AtomMask const& maskIn, 
                                     bool nofit, bool useMass, int debugIn)
{
  // TODO better error handles
  if (dIn == 0) {
    mprinterr("Internal Error: Metric_QuatRMSD::Init() called with null data set.\n");
    return 1;
  }
  coords_ = dIn;
  mask_ = maskIn;
  useMass_ = useMass;

  return 0;
}

int Cpptraj::Cluster::Metric_QuatRMSD::Setup() {
  if (coords_->Top().SetupIntegerMask( mask_ )) return 1;
  mask_.MaskInfo();
# ifdef DEBUG_CLUSTER
  mprintf("DEBUG: QuatRMSD metric topology: %s %s %i\n", coords_->legend(),
          coords_->Top().c_str(), coords_->Top().Natom());
# endif
  frm1_.SetupFrameFromMask(mask_, coords_->Top().Atoms());
  frm2_ = frm1_;
  return 0;
}

double Cpptraj::Cluster::Metric_QuatRMSD::FrameDist(int f1, int f2) {
  coords_->GetFrame( f1, frm1_, mask_ );
  coords_->GetFrame( f2, frm2_, mask_ );
  return QuaternionRMSD_CenteredRef(frm2_, frm1_,  useMass_); 
}

double Cpptraj::Cluster::Metric_QuatRMSD::CentroidDist(Centroid* c1, Centroid* c2) {
  // Centroid is already at origin.
  return QuaternionRMSD_CenteredRef( ((Centroid_Coord*)c2)->Cframe(),
                                     ((Centroid_Coord*)c1)->Cframe(), useMass_ );
}

double Cpptraj::Cluster::Metric_QuatRMSD::FrameCentroidDist(int f1, Centroid* c1) {
  coords_->GetFrame( f1, frm1_, mask_ );
  // Centroid is already at origin.
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
      //if (SRMSD_.Fit())
        cent->Cframe().CenterOnOrigin(useMass_);
    } else {
      Matrix_3x3 Rot;
      Vec3 TgtTrans;
      QuaternionRMSD_CenteredRef( cent->Cframe(), frm1_, Rot, TgtTrans, useMass_ );
      //if (SRMSD_.Fit()) {
        //frm1_.Translate( TgtTrans );
        frm1_.Rotate( Rot );
      //}
      cent->Cframe() += frm1_;
    }
  }
  cent->Cframe().Divide( (double)cframesIn.size() );
  //mprintf("\t\tFirst 3 centroid coords (of %i): %f %f %f\n", cent->cframe_.Natom(), 
  //        cent->cent->cframe_[0], cent->cframe_[1],cent->cframe_[2]);
}

Cpptraj::Cluster::Centroid* Cpptraj::Cluster::Metric_QuatRMSD::NewCentroid( Cframes const& cframesIn ) {
  // TODO: Incorporate mass?
  Centroid_Coord* cent = new Centroid_Coord( mask_.Nselected() );
  CalculateCentroid( cent, cframesIn );
  return cent;
}

void Cpptraj::Cluster::Metric_QuatRMSD::FrameOpCentroid(int frame, Centroid* centIn, double oldSize,
                                        CentOpType OP)
{
  Matrix_3x3 Rot;
  Vec3 Trans;
  Centroid_Coord* cent = (Centroid_Coord*)centIn;
  coords_->GetFrame( frame, frm1_, mask_ );
  QuaternionRMSD_CenteredRef( cent->Cframe(), frm1_, Rot, Trans, useMass_ );
  //if (SRMSD_.Fit()) {
    //frm2_.Translate( SRMSD_.TgtTrans() );
    frm1_.Rotate( Rot );
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

/** \return Description of quaternion RMS calc. */
std::string Cpptraj::Cluster::Metric_QuatRMSD::Description() const {
  std::string description("qrmsd " + mask_.MaskExpression());
  //if (!SRMSD_.Fit()) description.append(" nofit");
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
  //if (!SRMSD_.Fit())
  //  mprintf(", no fitting");
  //else
    mprintf(" best-fit");
  mprintf("\n");
}
