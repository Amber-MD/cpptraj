#include <cmath> // sqrt
#include "Metric_Data_Euclid.h"
#include "Centroid_Multi.h"

int Cpptraj::Cluster::Metric_Data_Euclid::Init(DsArray const& dsIn)
{
  for (DsArray::const_iterator ds = dsIn.begin(); ds != dsIn.end(); ++ds) {
    dsets_.push_back( (DataSet_1D*)*ds );
    if ( dsets_.back()->Meta().IsTorsionArray() )
      dcalcs_.push_back( DistCalc_Dih );
    else
      dcalcs_.push_back( DistCalc_Std );
  }
  return 0;
}

double Cpptraj::Cluster::Metric_Data_Euclid::FrameDist(int f1, int f2) {
  double dist = 0.0;
  DcArray::const_iterator dcalc = dcalcs_.begin();
  for (D1Array::iterator ds = dsets_.begin(); ds != dsets_.end(); ++ds, ++dcalc) {
    double diff = (*dcalc)((*ds)->Dval(f1), (*ds)->Dval(f2));
    dist += (diff * diff);
  }
  return sqrt(dist);
}

double Cpptraj::Cluster::Metric_Data_Euclid::CentroidDist(Centroid* c1, Centroid* c2) {
  double dist = 0.0;
  Centroid_Multi::Darray::const_iterator c2val = ((Centroid_Multi*)c2)->Cvals().begin();
  DcArray::const_iterator dcalc = dcalcs_.begin();
  for (Centroid_Multi::Darray::const_iterator c1val = ((Centroid_Multi*)c1)->Cvals().begin();
                                              c1val != ((Centroid_Multi*)c1)->Cvals().end();
                                            ++c1val, ++dcalc)
  {
    double diff = (*dcalc)(*c1val, *(c2val++));
    dist += (diff * diff);
  }
  return sqrt(dist);
}

double Cpptraj::Cluster::Metric_Data_Euclid::FrameCentroidDist(int f1, Centroid* c1) {
  double dist = 0.0;
  Centroid_Multi::Darray::const_iterator c1val = ((Centroid_Multi*)c1)->Cvals().begin();
  DcArray::const_iterator dcalc = dcalcs_.begin();
  for (D1Array::const_iterator ds = dsets_.begin(); ds != dsets_.end(); ++ds) {
    double diff = (*dcalc)((*ds)->Dval(f1), *(c1val++));
    dist += (diff * diff);
  }
  return sqrt(dist);
}

void Cpptraj::Cluster::Metric_Data_Euclid::CalculateCentroid(Centroid* centIn, Cframes const& cframesIn) {
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

Cpptraj::Cluster::Centroid* Cpptraj::Cluster::Metric_Data_Euclid::NewCentroid(Cframes const& cframesIn) {
  Centroid_Multi* cent = new Centroid_Multi();
  CalculateCentroid(cent, cframesIn);
  return cent;
}

//static const char* OPSTRING[] = {"ADD", "SUBTRACT"}; // DEBUG

void Cpptraj::Cluster::Metric_Data_Euclid::FrameOpCentroid(int frame, Centroid* centIn,
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

std::string Cpptraj::Cluster::Metric_Data_Euclid::Description() const {
  std::string description("data ");
  for (D1Array::const_iterator ds = dsets_.begin(); ds != dsets_.end(); ++ds)
    if (ds == dsets_.begin())
      description.append( (*ds)->Meta().PrintName() );
    else
      description.append( "," + (*ds)->Meta().PrintName() );
  return description;
}
