#include <cmath> // fabs, atan2, sin, cos
#include "Metric_Data.h"
#include "Centroid_Multi.h"
#include "../Constants.h" // RADDEG, DEGRAD
#include "../CpptrajStdio.h"

int Cpptraj::Cluster::Metric_Data::Init(DsArray const& dsIn)
{
  ntotal_ = 0;
  for (DsArray::const_iterator ds = dsIn.begin(); ds != dsIn.end(); ++ds) {
    if ( (*ds)->Group() != DataSet::SCALAR_1D) {
      mprinterr("Error: Set '%s' is not 1D scalar - cannot use for clustering.\n",
                (*ds)->legend());
      return 1;
    }
    dsets_.push_back( (DataSet_1D*)*ds );
    if ( dsets_.back()->Meta().IsTorsionArray() )
      dcalcs_.push_back( DistCalc_Dih );
    else
      dcalcs_.push_back( DistCalc_Std );
  }
  return 0;
}

/** Total number of points is number in smallest set. */
int Cpptraj::Cluster::Metric_Data::Setup() {
  ntotal_ = 0;
  for (D1Array::const_iterator ds = dsets_.begin(); ds != dsets_.end(); ++ds)
  {
    if (ntotal_ == 0)
      ntotal_ = (*ds)->Size();
    else {
      if ( (*ds)->Size() < ntotal_ ) {
        mprintf("Warning: Set '%s' size %zu is smaller than current size %u.\n",
                (*ds)->legend(), (*ds)->Size(), ntotal_);
        ntotal_ = (*ds)->Size();
      }
    }
  }
  return 0;
}

/// Calculate smallest difference between two angles (in degrees).
double Cpptraj::Cluster::Metric_Data::DistCalc_Dih(double d1, double d2) {
  double diff = fabs(d1 - d2);
  if (diff > 180.0)
    return (360.0 - diff);
  else
    return diff;
}

/// Calculate basic difference.
double Cpptraj::Cluster::Metric_Data::DistCalc_Std(double d1, double d2) {
  return fabs(d1 - d2);
}

/** Calculate unambiguous average dihedral angle (in degrees) by converting to 
  * cartesian coords using x = cos(theta), y = sin(theta), and:
  *   tan(avgtheta) = avgy / avgx = SUM[sin(theta)] / SUM[cos(theta)]
  * See Eq. 2 from Altis et al., J. Chem. Phys., 126 p. 244111 (2007).
  */
double Cpptraj::Cluster::Metric_Data::AvgCalc_Dih(DataSet_1D const& dsIn, 
                                                  Cframes const& cframesIn,
                                                  double& sumx, double& sumy)
{
  sumy = 0.0;
  sumx = 0.0;
  // TODO: Convert angles to radians prior to this call?
  for (Cframes::const_iterator frm = cframesIn.begin(); frm != cframesIn.end(); ++frm) {
    double theta = dsIn.Dval( *frm ) * Constants::DEGRAD;
    sumy += sin( theta );
    sumx += cos( theta );
  }
  return atan2(sumy, sumx) * Constants::RADDEG; 
}

double Cpptraj::Cluster::Metric_Data::AvgCalc_Std(DataSet_1D const& dsIn,
                                                  Cframes const& cframesIn)
{
  double val = 0.0;
  for (Cframes::const_iterator frm = cframesIn.begin(); frm != cframesIn.end(); ++frm)
    val += dsIn.Dval( *frm );
  return (val / (double)cframesIn.size());
}

/* Update centroid value for adding/removing a frame.
 * \param fval value of frame being added/removed.
 * \param cval current centroid value.
 * \param isTorsion data is periodic.
 * \param oldSize Previous size of the centroid.
 * \param OP Operation being performed.
 */
double Cpptraj::Cluster::Metric_Data::DistCalc_FrameCentroid(double fval, double cval,
                                     bool isTorsion,
                                     double oldSize, CentOpType OP,
                                     double& sumx, double& sumy)
{
  double newcval;
  if (isTorsion) {
    double ftheta = fval * Constants::DEGRAD;
    if (OP == ADDFRAME) {
      sumy += sin( ftheta );
      sumx += cos( ftheta );
    } else { // SUBTRACTFRAME
      sumy -= sin( ftheta );
      sumx -= cos( ftheta );
    }
    newcval = atan2(sumy, sumx) * Constants::RADDEG;
  } else {
    newcval = cval * oldSize;
    if (OP == ADDFRAME) {
      newcval += fval;
      newcval /= ( oldSize + 1 );
    } else { // SUBTRACTFRAME
      newcval -= fval;
      newcval /= ( oldSize - 1 );
    }
  }
  return newcval;
}

void Cpptraj::Cluster::Metric_Data::CalculateCentroid(Centroid* centIn, Cframes const& cframesIn) {
  Centroid_Multi* cent = (Centroid_Multi*)centIn;
  cent->Cvals().resize( dsets_.size(), 0.0 );
  cent->SumX().resize( dsets_.size(), 0.0 );
  cent->SumY().resize( dsets_.size(), 0.0 );
  for (unsigned int idx = 0; idx != dsets_.size(); ++idx) {
    if (dsets_[idx]->Meta().IsTorsionArray())
      cent->Cvals()[idx] = AvgCalc_Dih(*dsets_[idx], cframesIn,
                                       cent->SumX()[idx], cent->SumY()[idx]);
    else
      cent->Cvals()[idx] = AvgCalc_Std(*dsets_[idx], cframesIn);
  }
//  mprintf("DEBUG: Centroids:");
//  for (unsigned int i = 0; i != cent->cvals_.size(); i++)
//    mprintf("   %f (sumy=%f sumx=%f)", cent->cvals_[i], cent->Sumy_[i], cent->Sumx_[i]);
//  mprintf("\n");
}

Cpptraj::Cluster::Centroid* Cpptraj::Cluster::Metric_Data::NewCentroid(Cframes const& cframesIn) {
  Centroid_Multi* cent = new Centroid_Multi();
  CalculateCentroid(cent, cframesIn);
  return cent;
}

//static const char* OPSTRING[] = {"ADD", "SUBTRACT"}; // DEBUG

void Cpptraj::Cluster::Metric_Data::FrameOpCentroid(int frame, Centroid* centIn,
                                                           double oldSize, CentOpType OP)
{
  Centroid_Multi* cent = (Centroid_Multi*)centIn;
//  mprintf("DEBUG: Old Centroids:");
//  for (unsigned int i = 0; i != cent->cvals_.size(); i++)
//    mprintf("   sumy=%f sumx=%f", cent->Sumy_[i], cent->Sumx_[i]);
//    //mprintf(" %f", cent->cvals_[i]);
//  mprintf("\n");
  for (unsigned int i = 0; i != dsets_.size(); ++i)
    cent->Cvals()[i] = DistCalc_FrameCentroid(dsets_[i]->Dval(frame), 
                          cent->Cvals()[i], dsets_[i]->Meta().IsTorsionArray(), oldSize, OP,
                          cent->SumX()[i], cent->SumY()[i]);
//  mprintf("DEBUG: New Centroids after %s frame %i:", OPSTRING[OP], frame);
//  for (unsigned int i = 0; i != cent->cvals_.size(); i++)
//    mprintf(" %f", cent->cvals_[i]);
//  mprintf("\n");
}

std::string Cpptraj::Cluster::Metric_Data::SetNames(std::string const& descrip) const {
  std::string description(descrip);
  for (D1Array::const_iterator ds = dsets_.begin(); ds != dsets_.end(); ++ds)
    if (ds == dsets_.begin())
      description.append( (*ds)->Meta().PrintName() );
    else
      description.append( "," + (*ds)->Meta().PrintName() );
  return description;
}
