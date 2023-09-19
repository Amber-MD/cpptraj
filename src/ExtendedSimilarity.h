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
    enum MetricType { MSD = 0, ///< Mean-squared deviation
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
                      NO_METRIC };
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

    static const char* MetricStr_[];

    double msd_condensed(Darray const&, Darray const&, unsigned int, unsigned int) const;

    //Darray c_sum_;      ///< Hold sum over samples of each feature
    //Darray sq_sum_;     ///< Hold sum of squares over samples of each feature
};
// -----------------------------------------------------------------------------
/** Hold options for extended similarity. */
class ExtendedSimilarity::Opts {
  public:
    /// Blank constructor
    Opts() : metric_(NO_METRIC) {}
    /// MSD constructor - Sum of squares, number of atoms
    Opts(Darray const&, unsigned int);
    /// \return Metric type
    MetricType Metric() const { return metric_; }
    /// \return Sum of squares array (MSD)
    Darray const& Sq_sum() const { return *sq_sum_ptr_; }
    /// \return Number of atoms (MSD)
    unsigned int Natoms() const { return natoms_; }
  private:
    MetricType metric_;        ///< Desired metric
    Darray const* sq_sum_ptr_; ///< Pointer to sum of squares array (MSD)
    unsigned int natoms_;      ///< Number of atoms (MSD)
};
#endif
