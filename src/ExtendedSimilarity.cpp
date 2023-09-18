#include "ExtendedSimilarity.h"

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
  "",
  "Sokal-Sneath 1 coefficient",
  "Sokal-Sneath 2 coefficient",
  "No metric"
};
