#include "ExtendedSimilarity.h"
#include "DataSet_Coords.h"
#include "CpptrajStdio.h"
#include <cmath> // ceil, pow

/** CONSTRUCTOR for MSD - sum of squares array, number of atoms. */
ExtendedSimilarity::Opts::Opts(Darray const& sq_sum, unsigned int natoms) :
  metric_(MSD),
  sq_sum_ptr_(&sq_sum),
  natoms_(natoms)
{ }

/** CONSTRUCTOR - metric, c. threshold type, c. threshold value, weight type, weight value */
ExtendedSimilarity::Opts::Opts(MetricType mt, CoincidenceThresholdType ct, double cv,
                               WeightFactorType wt, double wv) :
  metric_(mt),
  cthreshType_(ct),
  c_threshold_(cv),
  wfactorType_(wt),
  power_(wv)
{ }

/** CONSTRUCTOR - metric only, defaults for the rest. */
ExtendedSimilarity::Opts::Opts(MetricType mt) :
  metric_(mt),
  cthreshType_(NO_THRESHOLD),
  c_threshold_(0),
  wfactorType_(FRACTION),
  power_(0)
{}

/** \return True if options are valid. */
bool ExtendedSimilarity::Opts::IsValid(unsigned int n_objects) const {
  if (metric_ == MSD) {
    if (sq_sum_ptr_ == 0 || sq_sum_ptr_->empty()) {
      mprinterr("Error: Similarity options are set up for MSD metric but sum of squares array is empty.\n");
      return false;
    }
    if (natoms_ < 1) {
      mprinterr("Error: Similarity options are set up for MSD metric but # atoms < 1.\n");
      return false;
    }
  } else if (metric_ == NO_METRIC) {
    mprinterr("Error: Similarity options metric not set.\n");
    return false;
  } else {
    if (cthreshType_ == N_OBJECTS) {
      if ((unsigned int)c_threshold_ >= n_objects) {
        mprinterr("Error: c_threshold cannot be equal or greater to n_objects.\n");
        return false;
      }
    } else if (cthreshType_ == FRAC_OBJECTS) {
      bool in_range = (c_threshold_ > 0 && c_threshold_ < 1);
      if (!in_range) {
        mprinterr("Error: c_threshold fraction must be between 0 and 1.\n");
        return false;
      }
    }
  }
  return true;
}

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
    case BUB : calculate_counters(c_sum, Nframes, opts); break; //FIXME
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

// Power weight functions
ExtendedSimilarity::Darray ExtendedSimilarity::f_s_power(Darray const& in, unsigned int n_objects, double p) {
  Darray out;
  out.reserve(in.size());
  for (Darray::const_iterator d = in.begin(); d != in.end(); ++d)
    out.push_back( pow(p, (double)(-(n_objects - *d))) );
  return out;
}

ExtendedSimilarity::Darray ExtendedSimilarity::f_d_power(Darray const& in, unsigned int n_objects, double p) {
  Darray out;
  out.reserve(in.size());
  for (Darray::const_iterator d = in.begin(); d != in.end(); ++d)
    out.push_back( pow(p, (double)(-(*d - (n_objects % 2)))) );
  return out;
}

ExtendedSimilarity::Darray ExtendedSimilarity::f_s_frac(Darray const& in, unsigned int n_objects, double p) {
  Darray out;
  out.reserve(in.size());
  for (Darray::const_iterator d = in.begin(); d != in.end(); ++d)
    out.push_back( *d / n_objects );
  return out;
}

ExtendedSimilarity::Darray ExtendedSimilarity::f_d_frac(Darray const& in, unsigned int n_objects, double p) {
  Darray out;
  out.reserve(in.size());
  for (Darray::const_iterator d = in.begin(); d != in.end(); ++d)
    out.push_back( (*d - (n_objects % 2)) / n_objects );
  return out;
}

ExtendedSimilarity::Darray ExtendedSimilarity::f_one(Darray const& in, unsigned int n_objects, double p) {
  Darray out(in.size(), 1.0);
  return out;
}

static inline void printDarray(std::vector<double> const& arr) {
  int col = 0;
  mprintf("[");
  for (std::vector<double>::const_iterator it = arr.begin(); it != arr.end(); ++it) {
    mprintf(" %10.8g", *it);
    col++;
    if (col == 6) {
      mprintf("\n");
      col = 0;
    }
  }
  mprintf("]\n");
}

static inline void printBarray(std::vector<bool> const& arr) {
  int col = 0;
  mprintf("[");
  for (std::vector<bool>::const_iterator it = arr.begin(); it != arr.end(); ++it) {
    if (*it)
      mprintf(" True");
    else
      mprintf(" False");
    col++;
    if (col == 12) {
      mprintf("\n");
      col = 0;
    }
  }
  mprintf("]\n");
}

static inline unsigned int Bsum(std::vector<bool> const& arr) {
  unsigned int count = 0;
  for (std::vector<bool>::const_iterator it = arr.begin(); it != arr.end(); ++it)
    if (*it) count++;
  return count;
}

ExtendedSimilarity::Darray ExtendedSimilarity::subArray(Darray const& d, Barray const& b, unsigned int n_objects)
{
  Darray out;
  out.reserve(d.size());
  for (unsigned int idx = 0; idx != d.size(); idx++)
    if (b[idx]) out.push_back(2 * d[idx] - n_objects);
  return out;
}

/** Calculate 1-similarity, 0-similarity, and dissimilarity counters.
  * \param c_total Column sum of the data (c_sum)
  * \param n_objects Number of samples (frames)
  * \param opts Extended similarity options
  */
int ExtendedSimilarity::calculate_counters(Darray const& c_total, unsigned int n_objects,
                                           Opts const& opts)
const
{
  // Assign c_threshold
  unsigned int c_threshold;
  switch (opts.CoincidenceThreshold()) {
    case NO_THRESHOLD : c_threshold = n_objects % 2; break;
    case DISSIMILAR   : c_threshold = ceil(n_objects / 2); break;
    case N_OBJECTS    : c_threshold = (unsigned int)opts.CoincidenceThresholdVal(); break;
    case FRAC_OBJECTS : c_threshold = (unsigned int)(opts.CoincidenceThresholdVal() * n_objects); break;
  }
  // Set w_factor
  double power = opts.WeightFactorPower();
  typedef Darray (*WgtFxnType)(Darray const&, unsigned int, double);
  WgtFxnType f_s; // Similarity function
  WgtFxnType f_d; // Dissimilarity function

  switch(opts.WeightFactor()) {
    case POWER:
      f_s = f_s_power;
      f_d = f_d_power;
      break;
    case FRACTION:
      f_s = f_s_frac;
      f_d = f_d_frac;
      break;
    case OTHER:
      f_s = f_one;
      f_d = f_one;
      break;
  }


  Barray a_indices;
  a_indices.reserve(c_total.size());
  for (Darray::const_iterator it = c_total.begin(); it != c_total.end(); ++it)
    a_indices.push_back( 2 * *it - n_objects > c_threshold );
  //printBarray( a_indices );

  Barray d_indices;
  d_indices.reserve(c_total.size());
  for (Darray::const_iterator it = c_total.begin(); it != c_total.end(); ++it)
    d_indices.push_back( n_objects - 2 * *it > c_threshold );
  //printBarray( d_indices );

  Barray dis_indices;
  dis_indices.reserve(c_total.size());
  for (Darray::const_iterator it = c_total.begin(); it != c_total.end(); ++it)
    dis_indices.push_back( fabs( 2 * *it - n_objects) <= c_threshold );
  //printBarray( dis_indices );

  unsigned int a_count = Bsum(a_indices);
  unsigned int d_count = Bsum(d_indices);
  unsigned int total_dis = Bsum(dis_indices);

  mprintf("%u %u %u\n", a_count, d_count, total_dis);

  Darray a_w_array = f_s( subArray(c_total, a_indices, n_objects), n_objects, power );
  printDarray( a_w_array );
  
  return 0;
}

