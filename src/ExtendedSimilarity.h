#ifndef INC_EXTENDED_SIMILARITY_H
#define INC_EXTENDED_SIMILARITY_H
#include <vector>
// Fwd declares
class DataSet_Coords;
/// Implements extended similarity comparisons
class ExtendedSimilarity {
  public:
    typedef std::vector<double> Darray;
    /// Metric types
    enum MetricType {
      MSD = 0, ///< Mean-squared deviation
      BUB,     ///< Bhattacharyya's U coefficient
      FAI,     ///< Faiman's coefficient
      GLE,     ///< Gleason's coefficient
      JA,      ///< Jaccard's coefficient
      JT,      ///< Jaccard-Tanimoto coefficient
      RT,      ///< Rogers-Tanimoto coefficient
      RR,      ///< Russell-Rao coefficient
      SM,      ///< Simpson's coefficient
      SS1,     ///< Sokal-Sneath 1 coefficient
      SS2,     ///< Sokal-Sneath 2 coefficient
      NO_METRIC
    };
    /// Coincidence threshold types
    enum CoincidenceThresholdType {
      NO_THRESHOLD = 0, ///< Default, c_threshold = n_objects % 2
      DISSIMILAR,       ///< Dissimilar, c_threshold = ceil(n_objects/2)
      N_OBJECTS,        ///< Target number of objects (< total number of objects)
      FRAC_OBJECTS      ///< Fraction of total number of objects
    };
    /// Weight factor types
    enum WeightFactorType {
      FRACTION = 0, ///< similarity = d[k]/n_objects, dissimilarity = 1 - (d[k] - n_objects % 2)/n_objects
      POWER,        ///< similarity = n^-(n_objects - d[k]), dissimilarity = n^-(d[k] - n_objects % 2)
      OTHER         ///< similarity = dissimilarity = 1
    };
    /// CONSTRUCTOR
    ExtendedSimilarity();
    /// Hold extended similarity options
    class Opts;
    /// \return Char string corresponding to given MetricType
    static const char* metricStr(MetricType);
    /// \return Extended comparison value for given COORDS set TODO c_threshold, w_factor
    //double Comparison(DataSet_Coords&, MetricType) const;
    /// \return Extended comparison value for given arrays
    double Comparison(Darray const&, unsigned int, Opts const&) const;
  private:
    typedef std::vector<bool> Barray;

    static const char* MetricStr_[];
    /// Calculate MSD from sum and squared sum arrays
    double msd_condensed(Darray const&, Darray const&, unsigned int, unsigned int) const;
    /// \return Sub-array based on values of given boolean array
    static inline Darray subArray(Darray const&, Barray const&, unsigned int);
    /// Calculate 1-similarity, 0-similarity, and dissimilarity counters from sum array
    int calculate_counters(Darray const&, unsigned int, Opts const&) const;

    static Darray f_s_power(Darray const&, unsigned int, double);
    static Darray f_d_power(Darray const&, unsigned int, double);
    static Darray f_s_frac(Darray const&, unsigned int, double);
    static Darray f_d_frac(Darray const&, unsigned int, double);
    static Darray f_one(Darray const&, unsigned int, double);
};
// -----------------------------------------------------------------------------
/** Hold options for extended similarity. */
class ExtendedSimilarity::Opts {
  public:
    /// Blank constructor
    Opts() : metric_(NO_METRIC) {}
    /// MSD constructor - Sum of squares, number of atoms
    Opts(Darray const&, unsigned int);
    /// Constructor - metric, c. threshold type, c. threshold value, weight type, weight value
    Opts(MetricType, CoincidenceThresholdType, double, WeightFactorType, double);
    /// Constructor - metric only (not MSD!)
    Opts(MetricType);
    /// \return Metric type
    MetricType Metric() const { return metric_; }
    /// \return Sum of squares array (MSD)
    Darray const& Sq_sum() const { return *sq_sum_ptr_; }
    /// \return Number of atoms (MSD)
    unsigned int Natoms() const { return natoms_; }

    /// \return Coincidence threshold type
    CoincidenceThresholdType CoincidenceThreshold() const { return cthreshType_; }
    /// \return Coincidence threshold value
    double CoincidenceThresholdVal() const { return c_threshold_; }
    /// \return Weight factor type
    WeightFactorType WeightFactor() const { return wfactorType_; }
    /// \return Weight factor power (for type POWER)
    double WeightFactorPower() const { return power_; }

    /// \return True if options are valid. Takes total number of objects (frames) to check
    bool IsValid(unsigned int) const;
  private:
    MetricType metric_;        ///< Desired metric
    Darray const* sq_sum_ptr_; ///< Pointer to sum of squares array (MSD)
    unsigned int natoms_;      ///< Number of atoms (MSD)
    CoincidenceThresholdType cthreshType_; ///< Coincidence threshold type
    double c_threshold_;                   ///< Coincidence threshold value
    WeightFactorType wfactorType_;         ///< Weight factor type
    double power_;                         ///< Power for weight factor power type
};
#endif
