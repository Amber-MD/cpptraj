#ifndef INC_EXTENDED_SIMILARITY_H
#define INC_EXTENDED_SIMILARITY_H
#include <vector>
/// Implements extended similarity comparisons
class ExtendedSimilarity {
  public:
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
  private:
    typedef std::vector<double> Darray;

    static const char* MetricStr_[];

    Darray c_sum_;      ///< Hold sum over samples of each feature
    Darray sq_sum_;     ///< Hold sum of squares over samples of each feature
};
#endif
