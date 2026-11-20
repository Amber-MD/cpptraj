#ifndef INC_CPPTRAJ_MDANCE_BTS_H
#define INC_CPPTRAJ_MDANCE_BTS_H
#ifdef HAS_EIGEN
//#include <cstdlib>
#include "esim.h"
#include "types.h"
namespace Cpptraj {
namespace Mdance {

double meanSqDev(const ArrayXXd& data, int nAtoms = 1);
double msdCondensed(const ArrayXd& cSum, const ArrayXd& sqSum, Index N, int nAtoms = 1);
double extendedComparison(const ArrayXXd& data, Index N = 0, int nAtoms = 1, bool isCondensed = false, MD::Metric mt = MD::Metric::MSD, MD::Threshold cThreshold = MD::Threshold(MD::ThresholdType::None, 1), int wPower = 0);
ArrayXd calculateCompSim(const ArrayXXd& data, int nAtoms = 1, MD::Metric mt = MD::Metric::MSD);
Index calculateMedoid(const ArrayXXd& data, int nAtoms = 1, MD::Metric mt = MD::Metric::MSD);
Index calculateMedoid(const ArrayXd& data);
Index calculateOutlier(const ArrayXXd& data, int nAtoms = 1, MD::Metric mt = MD::Metric::MSD);
Index calculateOutlier(const ArrayXd& data);
ArrayXXd trimOutliers(const ArrayXXd& data, int nTrimmed, int nAtoms = 1, bool isMedoid = false, MD::Metric mt = MD::Metric::MSD);
ArrayXXd trimOutliers(const ArrayXXd& data, float nTrimmed, int nAtoms = 1, bool isMedoid = false, MD::Metric mt = MD::Metric::MSD);
vector<Index> diversitySelection(const ArrayXXd& data, int percentage, MD::Metric mt = MD::Metric::MSD, int nAtoms = 1, bool isCompSim = false, MD::StartSeed start = MD::StartSeed::Medoid);
vector<Index> diversitySelection(const ArrayXXd& data, int percentage, MD::Metric mt, int nAtoms, vector<Index>& start);
Index getNewIndexN(const ArrayXXd& data, MD::Metric mt, ArrayXXd& selectedCondensed, int N, set<Index>& selectFromN, int nAtoms);
ArrayXi repSample(const ArrayXXd& data, MD::Metric mt = MD::Metric::MSD, int nAtoms = 1, int nBins = 10, int nSamples = 100, bool hardCap = true);
ArrayXi repSample(const ArrayXXd& data, MD::Metric mt, int nAtoms, int nBins, double nSamples, bool hardCap = true);
ArrayXXd refineDisMatrix(const ArrayXXd& data);
ArrayXXd alignTraj(const ArrayXXd& data, int nAtoms, MD::AlignMethod = MD::AlignMethod::None);
} /* END namespace Mdance */
} /* END namespace Cpptraj */
#endif /* HAS_EIGEN */
#endif
