#ifndef INC_EXTENDED_SIMILARITY_H
#define INC_EXTENDED_SIMILARITY_H
#include <vector>
#include <string>
// Fwd declares
class DataSet_Coords;
class Frame;
/// Implements extended similarity comparisons
class ExtendedSimilarity {
  public:
    typedef std::vector<double> Darray;
    /// Metric types. Sync with MetricStr_ and MetricKeys_
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
    /// Set options - metric, c. threshold type, c. threshold value, weight type, weight value, #frames, #coords
    int SetOpts(MetricType, CoincidenceThresholdType, double, WeightFactorType, double, unsigned int, unsigned int);
    /// Set options - metric, #frames, #coords
    int SetOpts(MetricType, unsigned int, unsigned int);

    /// \return Char string corresponding to given MetricType
    static const char* metricStr(MetricType);
    /// \return Type corresponding to given keyword
    static MetricType TypeFromKeyword(std::string const&);
    /// \return string containing all metric keywords
    static std::string MetricKeys();

    /// Calculate complimentary similarity over given COORDS set 
    Darray CalculateCompSim(DataSet_Coords&);
    /// Calculate complimentary similarity between two frames 
    double CalculateCompSim(Frame const&, Frame const&);

    /// \return Medoid frame index from last call to Darray CalculateCompSim
    long int MedoidIndex() const { return max_dissim_idx_; }
  private:
    typedef std::vector<bool> Barray;
    /// Hold counters from calculate_counters
    class Counters {
      public:
        Counters() : a_(0), w_a_(0), d_(0), w_d_(0), total_sim_(0), total_w_sim_(0),
                     total_dis_(0), total_w_dis_(0), p_(0), w_p_(0) {}

        double a_;
        double w_a_;
        double d_;
        double w_d_;
        double total_sim_;
        double total_w_sim_;
        double total_dis_;
        double total_w_dis_;
        double p_;
        double w_p_;
    };
    /// Descriptive strings corresponding to MetricType
    static const char* MetricStr_[];
    /// Keywords corresponding to MetricType
    static const char* MetricKeys_[];

    /// \return 1 if set up is invalid
    int isValid(unsigned int);
    /// \return Extended comparison value for given arrays
    double Comparison(Darray const&, unsigned int) const;
    /// Calculate MSD from sum and squared sum arrays
    double msd_condensed(Darray const&, Darray const&, unsigned int, unsigned int) const;
    /// \return Sub-array based on values of given boolean array
    static inline Darray subArray(Darray const&, Barray const&, unsigned int);
    /// \return Absolute value sub-array based on values of given boolean array
    static inline Darray absSubArray(Darray const&, Barray const&, unsigned int);
    /// Calculate 1-similarity, 0-similarity, and dissimilarity counters from sum array
    Counters calculate_counters(Darray const&, unsigned int) const;

    static Darray f_s_power(Darray const&, unsigned int, double);
    static Darray f_d_power(Darray const&, unsigned int, double);
    static Darray f_s_frac(Darray const&, unsigned int, double);
    static Darray f_d_frac(Darray const&, unsigned int, double);
    static Darray f_one(Darray const&, unsigned int, double);

    MetricType metric_;        ///< Desired metric
    Darray c_sum_;             ///< Sum array
    Darray sq_sum_;            ///< Sum of squares array (MSD)
    unsigned int natoms_;      ///< Number of atoms (MSD)
    CoincidenceThresholdType cthreshType_; ///< Coincidence threshold type
    double c_threshold_;                   ///< Coincidence threshold value
    WeightFactorType wfactorType_;         ///< Weight factor type
    double power_;                         ///< Power for weight factor power type
    double max_dissim_val_;                ///< Max dissimilarity value (medioid)
    long int max_dissim_idx_;              ///< Max dissimilarity index (medioid)
};
#endif
