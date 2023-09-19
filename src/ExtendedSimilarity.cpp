#include "ExtendedSimilarity.h"
#include "DataSet_Coords.h"
#include "CpptrajStdio.h"

/** CONSTRUCTOR for MSD - sum of squares array, number of atoms. */
ExtendedSimilarity::Opts::Opts(Darray const& sq_sum, unsigned int natoms) :
  metric_(MSD),
  sq_sum_ptr_(&sq_sum),
  natoms_(natoms)
{ }

// -----------------------------------------------------------------------------

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

/** \return Character string corresponding to given metric type. */
const char* ExtendedSimilarity::metricStr(MetricType m) {
  return MetricStr_[m];
}

/** \return Extended comparison value for COORDS set. */
/*double ExtendedSimilarity::Comparison(DataSet_Coords& crdIn, MetricType metricIn)
const
{
  unsigned int Ncoords = crdIn.Top().Natom() * 3;
  unsigned int Nelements = crdIn.Size() * Ncoords;
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
  return ExtendedSimilarity::Comparison(c_sum, sq_sum_total, metricIn,
                                        Nelements-1, crdIn.Top().Natom());
}*/

/** \return Extended comparison value.
  * \param c_sum Column sum of the data
  * \param sq_sum Column sum of the squared data (MSD only)
  * \param metricIn Similarity metric to use
  * \param Nframes Number of samples (frames)
  * \param Natoms Number of atoms (MSD only)
  */
double ExtendedSimilarity::Comparison(Darray const& c_sum, unsigned int Nframes, 
                                      Opts const& opts)
const
{

  double val = 0;
  switch (opts.Metric()) {
    case MSD : val = msd_condensed(c_sum, opts.Sq_sum(), Nframes, opts.Natoms()); break;
    default:
      mprinterr("Internal Error: ExtendedSimilarity::Comparison(): Metric '%s' is unhandled.\n",
                MetricStr_[opts.Metric()]);
  }

  return val;
}

/** Mean-squared deviation
  * \param c_sum (Natoms*3) Column sum of the data
  * \param sq_sum (Natoms*3) Column sum of the squared data
  * \param Nframes Number of samples (frames) 
  * \param Natoms Number of atoms in the system
  */
double ExtendedSimilarity::msd_condensed(Darray const& c_sum, Darray const& sq_sum, unsigned int Nframes, unsigned int Natoms)
const
{
  if (c_sum.size() != sq_sum.size()) {
    mprinterr("Internal Error: ExtendedSimilarity::msd_condensed(): Array sizes are not equal.\n");
    return 0;
  }
  //msd = np.sum(2 * (N * sq_sum - c_sum ** 2)) / (N * (N - 1))
  double msd = 0;
  for (unsigned int idx = 0; idx != c_sum.size(); idx++)
    msd += (2 * (Nframes * sq_sum[idx] - (c_sum[idx] * c_sum[idx])));
  msd /= ((double)Nframes * (double(Nframes-1)));
  //norm_msd = msd / N_atoms
  msd /= Natoms;
  return msd;
}

/** Calculate 1-similarity, 0-similarity, and dissimilarity counters.
  * \param c_total Column sum of the data (c_sum)
  * \param n_objects Number of samples (frames)
  * \param c_threshold 
  */
//int ExtendedSimilarity::calculate_counters(Darray const& c_total, unsigned int n_objects,
//                                           int c_threshold, double cutoff, double w_factor)
//{

