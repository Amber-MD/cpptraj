#include "ExtendedSimilarity.h"
#include "DataSet_Coords.h"

/** CONSTRUCTOR */
ExtendedSimilarity::ExtendedSimilarity()
{}

/** Strings corresponding to MetricType */
const char* ExtendedSimilarity::MetricStr_[] = {
  "Mean-squared deviation",
  "Bhattacharyya's U coefficient",
  "Faiman's coefficient",
  "Gleason's coefficient",
  "Jaccard's coefficient",
  "Jaccard-Tanimoto coefficient",
  "Rogers-Tanimoto coefficient",
  "Russell-Rao coefficient",
  "Simpson's coefficient",
  "Sokal-Sneath 1 coefficient",
  "Sokal-Sneath 2 coefficient",
  "No metric"
};

/** \return Extended comparison value for COORDS set. */
double ExtendedSimilarity::Comparison(DataSet_Coords& crdIn, MetricType metricIn)
const
{
  unsigned int Ncoords = crdIn.Top().Natom() * 3;
  //unsigned int Nelements = crdIn.Size() * Ncoords;
  Darray c_sum( Ncoords, 0.0 );
  Darray sq_sum_total( Ncoords, 0.0 );
  Frame frmIn = crdIn.AllocateFrame();
  // Get sum and sum squares for each coordinate
  for (unsigned int idx = 0; idx < crdIn.Size(); idx++) {
    crdIn.GetFrame(idx, frmIn);
    for (unsigned int icrd = 0; icrd < Ncoords; icrd++) {
      c_sum[icrd] = frmIn[icrd];
      sq_sum_total[icrd] = frmIn[icrd] * frmIn[icrd];
    }
  }
  return ExtendedSimilarity::Comparison(c_sum, sq_sum_total, metricIn);
}

/** \return Extended comparison value. */
double ExtendedSimilarity::Comparison(Darray const& c_sum, Darray const& sq_sum, MetricType metricIn)
const
{
  return 0;
}
