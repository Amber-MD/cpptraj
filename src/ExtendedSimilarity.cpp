#include "ExtendedSimilarity.h"
#include "DataSet_Coords.h"
#include "CpptrajStdio.h"
#include <cmath> // ceil, pow, sqrt

/** CONSTRUCTOR */
ExtendedSimilarity::ExtendedSimilarity() :
  metric_(NO_METRIC),
  cthreshType_(NO_THRESHOLD),
  c_threshold_(0),
  wfactorType_(FRACTION),
  power_(0)
{}

/** Set options - metric, c. threshold type, c. threshold value, weight type, weight value,
  * nframes, ncoords.
  */
int ExtendedSimilarity::SetOpts(MetricType mt, CoincidenceThresholdType ct, double cv,
                                WeightFactorType wt, double wv, unsigned int Nframes,
                                unsigned int Ncoords)
{
  metric_ = mt;
  cthreshType_ = ct;
  c_threshold_ = cv;
  wfactorType_ = wt;
  power_ = wv;
  c_sum_.assign( Ncoords, 0.0 );
  return isValid(Nframes);
}

/** Set options - metric, nframes, ncoords only; defaults for the rest. */
int ExtendedSimilarity::SetOpts(MetricType mt, unsigned int Nframes, unsigned int Ncoords)
{
  metric_ = mt;
  cthreshType_ = NO_THRESHOLD;
  c_threshold_ = 0;
  wfactorType_ = FRACTION;
  power_ = 0;
  c_sum_.assign( Ncoords, 0.0 );
  return isValid(Nframes);
}

/** Check options and do any common setup.
  * \return 0 if options are valid.
  */
int ExtendedSimilarity::isValid(unsigned int n_objects) {
  if (metric_ == MSD) {
    if (c_sum_.size() % 3 != 0) {
      mprinterr("Error: # of coords (%zu) is not divisible by 3; required for MSD.\n", c_sum_.size());
      return 1;
    }
    natoms_ = c_sum_.size() / 3;
    if (natoms_ < 1) {
      mprinterr("Error: Similarity options are set up for MSD metric but # atoms < 1.\n");
      return 1;
    }
    sq_sum_.assign( c_sum_.size(), 0.0 );
  } else if (metric_ == NO_METRIC) {
    mprinterr("Error: Similarity options metric not set.\n");
    return 1;
  } else {
    if (cthreshType_ == N_OBJECTS) {
      if ((unsigned int)c_threshold_ >= n_objects) {
        mprinterr("Error: c_threshold cannot be equal or greater to n_objects.\n");
        return 1;
      }
    } else if (cthreshType_ == FRAC_OBJECTS) {
      bool in_range = (c_threshold_ > 0 && c_threshold_ < 1);
      if (!in_range) {
        mprinterr("Error: c_threshold fraction must be between 0 and 1.\n");
        return 1;
      }
    }
  }
  return 0;
}

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

/** Keywords corresponding to MetricType. */
const char* ExtendedSimilarity::MetricKeys_[] = {
  "msd",
  "bub",
  "fai",
  "gle",
  "ja",
  "jt",
  "rt",
  "rr",
  "sm",
  "ss1",
  "ss2",
  0
};

/** \return Character string corresponding to given metric type. */
const char* ExtendedSimilarity::metricStr(MetricType m) {
  return MetricStr_[m];
}

/** \return MetricType corresponding to keyword. */
ExtendedSimilarity::MetricType ExtendedSimilarity::TypeFromKeyword(std::string const& key) {
  for (int i = 0; i < (int)NO_METRIC; i++) {
    if (key == std::string(MetricKeys_[i])) {
      return (MetricType)i;
    }
  }
  return NO_METRIC;
} 

/// For debug, print double array.
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

/// For debug, print Frame
static inline void printFrame(Frame const& arr) {
  int col = 0;
  mprintf("[");
  for (int idx = 0; idx != arr.size(); idx++) {
    mprintf(" %10.8g", arr[idx]);
    col++;
    if (col == 6) {
      mprintf("\n");
      col = 0;
    }
  }
  mprintf("]\n");
}
/** Calculate the comparitive similarity value between two frames. */
double ExtendedSimilarity::CalculateCompSim(Frame const& f1, Frame const& f2)
{
  // Sanity check
  if (metric_ == NO_METRIC || c_sum_.empty()) {
    mprinterr("Internal Error: ExtendedSimilarity::Comparison() called before SetOpts().\n");
    return 0;
  }
  //mprintf("[");
  //printFrame(f1);
  //printFrame(f2);
  //mprintf("\n");
  double val = 0;
  unsigned int Ncoords = c_sum_.size();
  //c_sum_.assign( Ncoords, 0.0 );
  for (unsigned int icrd = 0; icrd < Ncoords; icrd++)
    c_sum_[icrd] = f1[icrd] + f2[icrd];
  if (metric_ == MSD) {
    for (unsigned int icrd = 0; icrd < Ncoords; icrd++)
      sq_sum_[icrd] = (f1[icrd] * f1[icrd]) + (f2[icrd] * f2[icrd]);
    val = msd_condensed(c_sum_, sq_sum_, 2, natoms_);
  } else {
    val = Comparison(c_sum_, 2);
  }
  //mprintf("DEBUG val=%10.8g\n",val);

  return val;
}

/** Calculate the comparitive similarity values for COORDS set. */
ExtendedSimilarity::Darray ExtendedSimilarity::CalculateCompSim(DataSet_Coords& crdIn)
{
  // Sanity check
  if (metric_ == NO_METRIC || c_sum_.empty()) {
    mprinterr("Internal Error: ExtendedSimilarity::Comparison() called before SetOpts().\n");
    return Darray();
  }

  unsigned int Nframes = crdIn.Size();
  unsigned int Ncoords = c_sum_.size();
  Frame frmIn = crdIn.AllocateFrame();
  c_sum_.assign( Ncoords, 0.0 );
  if (metric_ == MSD) {
    sq_sum_.assign( Ncoords, 0.0 );
    // Get sum and sum squares for each coordinate
    for (unsigned int idx = 0; idx < Nframes; idx++) {
      crdIn.GetFrame(idx, frmIn);
      for (unsigned int icrd = 0; icrd < Ncoords; icrd++) {
        c_sum_[icrd]  += frmIn[icrd];
        sq_sum_[icrd] += frmIn[icrd] * frmIn[icrd];
      }
    }
  } else {
    // Get sum for each coordinate
    for (unsigned int idx = 0; idx < Nframes; idx++) {
      crdIn.GetFrame(idx, frmIn);
      for (unsigned int icrd = 0; icrd < Ncoords; icrd++)
        c_sum_[icrd] += frmIn[icrd];
    }
  }

  // For each frame, get the comp. similarity.
  Darray comp_sims;
  comp_sims.reserve( Nframes );
  Darray c_arr(c_sum_.size(), 0.0);
  if (metric_ == MSD) {
    Darray sq_arr(sq_sum_.size(), 0.0);
    for (unsigned int idx = 0; idx < Nframes; idx++) {
      crdIn.GetFrame(idx, frmIn);
      for (unsigned int icrd = 0; icrd < Ncoords; icrd++) {
        c_arr[icrd]  = c_sum_[icrd] - frmIn[icrd];
        sq_arr[icrd] = sq_sum_[icrd] - (frmIn[icrd]*frmIn[icrd]);
      }
      //printDarray(sq_arr);
      //printDarray(c_sum);
      double val = msd_condensed(c_arr, sq_arr, Nframes-1, natoms_);
      //mprintf("%8u %16.8f\n", idx, val);
      comp_sims.push_back( val );
    }
  } else {
    for (unsigned int idx = 0; idx < Nframes; idx++) {
      crdIn.GetFrame(idx, frmIn);
      for (unsigned int icrd = 0; icrd < Ncoords; icrd++)
        c_arr[icrd]  = c_sum_[icrd] - frmIn[icrd];
      //printDarray(c_sum);
      double val = Comparison(c_arr, Nframes-1);
      //mprintf("%8u %16.8f\n", idx, val);
      comp_sims.push_back( val );
    }
  }
  // Record the medioid (max dissimilarity)
  max_dissim_val_ = -1;
  max_dissim_idx_ = -1;
  for (Darray::const_iterator it = comp_sims.begin(); it != comp_sims.end(); ++it) {
    if (*it > max_dissim_val_) {
      max_dissim_val_ = *it;
      max_dissim_idx_ = it - comp_sims.begin();
    }
  }
  mprintf("DEBUG: Max dissim. val %g at idx %li\n", max_dissim_val_, max_dissim_idx_);

  return comp_sims;
}

/** \return Extended comparison value. Should NOT be called for MSD.
  * \param c_sum Column sum of the data
  * \param Nframes Number of samples (frames)
  */
double ExtendedSimilarity::Comparison(Darray const& c_sum, unsigned int Nframes) 
const
{
  double val = 0;
  Counters count;
  count = calculate_counters(c_sum, Nframes);
  //mprintf("%10.8g %10.8g %10.8g %10.8g %10.8g\n", count.w_a_, count.w_d_, count.a_, count.d_, count.total_dis_);
  switch (metric_) {
    case MSD :
      mprinterr("Internal Error: ExtendedSimilarity::Comparison() called for MSD metric.\n");
      break;
    case BUB :
      val = (sqrt(count.w_a_ * count.w_d_) + count.w_a_) / 
            (sqrt(count.a_ * count.d_) + count.a_ + count.total_dis_);
      break;
    case FAI :
      val = (count.w_a_ + 0.5 * count.w_d_) / (count.p_); break;
    case GLE :
      val = (2 * count.w_a_) / (2 * count.a_ + count.total_dis_); break;
    case JA :
      val = (3 * count.w_a_) / (3 * count.a_ + count.total_dis_); break;
    case JT :
      val = (count.w_a_) / (count.a_ + count.total_dis_); break;
    case RT :
      val = (count.total_w_sim_) / (count.p_ + count.total_dis_); break;
    case RR :
      val = (count.w_a_) / (count.p_); break;
    case SM :
      val = (count.total_w_sim_) / (count.p_); break;
    case SS1 :
      val = (count.w_a_) / (count.a_ + 2 * count.total_dis_); break;
    case SS2 :
      val = (2 * count.total_w_sim_) / (count.p_ + count.total_sim_); break;
    default:
      mprinterr("Internal Error: ExtendedSimilarity::Comparison(): Metric '%s' is unhandled.\n",
                MetricStr_[metric_]);
  }
  // NOTE 1 - val is invalid for MSD, but we should not be here for MSD
  return 1 - val;
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

// -------------------------------------
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

// -------------------------------------
/// For debug, print boolean array
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

/// \return count of true elements in boolean array
static inline unsigned int Bsum(std::vector<bool> const& arr) {
  unsigned int count = 0;
  for (std::vector<bool>::const_iterator it = arr.begin(); it != arr.end(); ++it)
    if (*it) count++;
  return count;
}

/// \return sum over double array
static inline double Dsum(std::vector<double> const& arr) {
  double sum = 0;
  for (std::vector<double>::const_iterator it = arr.begin(); it != arr.end(); ++it)
    sum += *it;
  return sum;
}

/** \return Array containing elements of d for which b is true, times 2 minus n_objects */
ExtendedSimilarity::Darray ExtendedSimilarity::subArray(Darray const& d, Barray const& b, unsigned int n_objects)
{
  Darray out;
  out.reserve(d.size());
  for (unsigned int idx = 0; idx != d.size(); idx++)
    if (b[idx]) out.push_back(2 * d[idx] - n_objects);
  return out;
}

/** \return Array containing absolute value of elements of d for which b is true, times 2 minus n_objects */
ExtendedSimilarity::Darray ExtendedSimilarity::absSubArray(Darray const& d, Barray const& b, unsigned int n_objects)
{
  Darray out;
  out.reserve(d.size());
  for (unsigned int idx = 0; idx != d.size(); idx++)
    if (b[idx]) out.push_back( fabs(2 * d[idx] - n_objects) );
  return out;
}

/** Calculate 1-similarity, 0-similarity, and dissimilarity counters.
  * \param c_total Column sum of the data (c_sum)
  * \param n_objects Number of samples (frames)
  */
ExtendedSimilarity::Counters
   ExtendedSimilarity::calculate_counters(Darray const& c_total, unsigned int n_objects)
const
{
  // Assign c_threshold
  unsigned int c_threshold;
  switch (cthreshType_) {
    case NO_THRESHOLD : c_threshold = n_objects % 2; break;
    case DISSIMILAR   : c_threshold = ceil(n_objects / 2); break;
    case N_OBJECTS    : c_threshold = (unsigned int)c_threshold_; break;
    case FRAC_OBJECTS : c_threshold = (unsigned int)(c_threshold_ * n_objects); break;
  }
  // Set w_factor
  typedef Darray (*WgtFxnType)(Darray const&, unsigned int, double);
  WgtFxnType f_s; // Similarity function
  WgtFxnType f_d; // Dissimilarity function

  switch(wfactorType_) {
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
  // A indices
  Barray a_indices;
  a_indices.reserve(c_total.size());
  for (Darray::const_iterator it = c_total.begin(); it != c_total.end(); ++it)
    a_indices.push_back( 2 * *it - n_objects > c_threshold );
  //printBarray( a_indices );
  // D indices
  Barray d_indices;
  d_indices.reserve(c_total.size());
  for (Darray::const_iterator it = c_total.begin(); it != c_total.end(); ++it)
    d_indices.push_back( n_objects - 2 * *it > c_threshold );
  //printBarray( d_indices );
  // Dissimilarity indices
  Barray dis_indices;
  dis_indices.reserve(c_total.size());
  for (Darray::const_iterator it = c_total.begin(); it != c_total.end(); ++it)
    dis_indices.push_back( fabs( 2 * *it - n_objects) <= c_threshold );
  //printBarray( dis_indices );

  Counters count;
  count.a_         = Bsum(a_indices);
  count.d_         = Bsum(d_indices);
  count.total_dis_ = Bsum(dis_indices);

  //mprintf("%g %g %g\n", count.a_, count.d_, count.total_dis_);

  Darray a_w_array = f_s( subArray(c_total, a_indices, n_objects), n_objects, power_ );
  //printDarray( a_w_array );
  Darray d_w_array = f_s( absSubArray(c_total, d_indices, n_objects), n_objects, power_ );
  //printDarray( d_w_array );
  Darray total_w_dis_array = f_d( absSubArray(c_total, dis_indices, n_objects), n_objects, power_ );
  //printDarray( total_w_dis_array );

  count.w_a_         = Dsum( a_w_array );
  count.w_d_         = Dsum( d_w_array );
  count.total_w_dis_ = Dsum( total_w_dis_array );
  //mprintf("%10.8f %10.8f %10.8f\n", count.w_a_, count.w_d_, count.total_w_dis_);

  count.total_sim_   = count.a_ + count.d_;
  count.total_w_sim_ = count.w_a_ + count.w_d_;
  count.p_           = count.total_sim_ + count.total_dis_;
  count.w_p_         = count.total_w_sim_ + count.total_w_dis_;
  //mprintf("%8g %10.8f %8g %10.8f\n", count.total_sim_, count.total_w_sim_, count.p_, count.w_p_);
  
  return count;
}
