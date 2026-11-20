#include "bts.h"
#ifdef HAS_EIGEN
using namespace Cpptraj::Mdance;
/* O(N) Mean square deviation(MSD) calculation for n-ary objects.
 *  
 * Parameters:
 *  - data: A feature array of size (nSamples, nFeatures)
 *  - nAtoms: Number of atoms in the Molecular Dynamics (MD) system. nAtoms=1 for non-MD systems.
 * 
 * Returns: normalized MSD value
 * 
 * Reference: https://github.com/mqcomplab/MDANCE/blob/9a895e72d71fee1d1a4fad1700a806473dff2f71/src/mdance/tools/bts.py#L14
*/ 
double Cpptraj::Mdance::meanSqDev(const ArrayXXd& data, int nAtoms){
    Index N = data.rows();
    if (N == 1)
        return 0;
    ArrayXd cSum = data.colwise().sum();
    ArrayXd sqSum = data.square().colwise().sum();
    return msdCondensed(cSum, sqSum, N, nAtoms);
}

/* Condensed version of Mean square deviation (MSD) calculation for n-ary objects
 * 
 * Parameters:
 *  - cSum: A feature array of the column-wise sum of the data (nFeatures)
 *  - sqSum: A feature array of the column-wise sum of the squared data (nFeatures)
 *  - N: Number of data points
 *  - nAtoms: Number of atoms in the Molecular Dynamics (MD) system. nAtoms=1 for non-MD systems.
 * 
 * Returns: normalized MSD value.
 * 
 * Reference: https://github.com/mqcomplab/MDANCE/blob/9a895e72d71fee1d1a4fad1700a806473dff2f71/src/mdance/tools/bts.py#L54
*/
double Cpptraj::Mdance::msdCondensed(const ArrayXd& cSum, const ArrayXd& sqSum, Index N, int nAtoms){
    if (N == 1)
        return 0;
    // The following is a step-by-step explanation of what we are returning. May need to use this instead if we run into overflow issues?
    // ArrayXd meanSqSum = sqSum / N;
    // ArrayXd meanCSum = cSum / N;
    // return msd = (meanSqSum - meanCSum.square()).sum() * 2.0 / nAtoms;
    if (N <= 0) {
        throw std::invalid_argument("N must be positive and non-zero.");
    }
    if (nAtoms <= 0) {
        throw std::invalid_argument("nAtoms must be positive and non-zero.");
    }
    if (N > 46340){
        return (double)2.0 * (sqSum * N - cSum.square()).sum() / ((double)N * N * nAtoms);
    }
    else{
        return (double)2.0 * (sqSum * N - cSum.square()).sum() / (N * N * nAtoms);
    }
}

/* O(N) Extended comparison function for n-ary objects.
 *
 * Parameters:
 * - data: A feature array which can take on multiple formats
 *    --> if (!isCondensed): an array of size (nSamples, nFeatures)
 *    --> if (isCondensed): a ArrayXXd with 1 row (cSum)
 *    --> if (isCondensed): a ArrayXXd with 2 rows (cSum, sqSum)
 * - isCondensed: Controls type of data (see above)
 * - mt: The metric to use when calculating distance between n objects in an array
 * - N: Number of data points.
 * - nAtoms: Number of atoms in the Molecular Dynamics (MD) system. nAtoms=1 for non-MD systems
 * - cThreshold: CCoincidence threshold for calculatinig extended similarity
 * - wPower: Controls the type of weight function for calculating extended similarity
 *    --> if (wPower): use fraction method
 *    --> else: use powerN method where N=wPower
 * 
 * Returns: Extended comparison value.
 * 
 * Reference: https://github.com/mqcomplab/MDANCE/blob/main/src/mdance/tools/bts.py#L96
*/
double Cpptraj::Mdance::extendedComparison(const ArrayXXd& data, Index N, int nAtoms, bool isCondensed, MD::Metric mt, MD::Threshold cThreshold, int wPower) {
    // Handle default initialization of MD::Threshold
    if (cThreshold.type == MD::ThresholdType::None)
        cThreshold.value = N % 2;
    // Data check
    if (isCondensed){
        if (data.rows() > 2){
            throw std::runtime_error("Data must have at most two rows: either (cSum) or (cSum, sqSum)");
        }
        ArrayXd cSum = data.row(0);
        if (mt == MD::Metric::MSD){
            ArrayXd sqSum = data.row(1);
            return msdCondensed(cSum, sqSum, N, nAtoms);
        } else {
            MD::Indices idx = genSimIdx(cSum, N, cThreshold, wPower);
            return 1 - idx.getIndex(mt);
        }
    } else {
        if (mt == MD::Metric::MSD) {
            return meanSqDev(data, nAtoms);
        } else {
            ArrayXd cSum = data.colwise().sum();
            MD::Indices idx = genSimIdx(cSum, N, cThreshold, wPower);
            return 1 - idx.getIndex(mt);
        }
    }

}

/* O(N) Complementary similarity calculation for n-ary objects.
 *
 * Parameters:
 *  - data: A feature array of shape (nSamples, nFeatures).
 *  - nAtoms: Number of atoms in the Molecular Dynamics (MD) system. nAtoms=1 for non-MD systems.
 *  - mt: The metric to use when calculating distance between n objects in an array.
 * 
 * Returns: Vector (N) of complementary similarities for each object
 * 
 * Reference: https://github.com/mqcomplab/MDANCE/blob/main/src/mdance/tools/bts.py#L190
*/
ArrayXd Cpptraj::Mdance::calculateCompSim(const ArrayXXd& data, int nAtoms, MD::Metric mt) {
    Index N = data.rows();

    ArrayXXd sqData = data.square();
    ArrayXd cSum = data.colwise().sum();
    ArrayXd sqSum = sqData.colwise().sum();

    ArrayXd compSims(N);

    if (mt == MD::Metric::MSD){
        ArrayXXd compC = ((-data).rowwise()+cSum.transpose()) / (N-1);
        ArrayXXd compSq = ((-sqData).rowwise()+sqSum.transpose()) / (N-1);
        compSims = (2 * (compSq - compC.square())/ nAtoms).rowwise().sum();
    } else {
        for (int i=0; i<N; ++i){
            ArrayXXd compData (2,data.cols());
            compData.row(0) = cSum.transpose() - data.row(i);
            compData.row(1) = sqSum - sqData.row(i);
            compSims[i] = extendedComparison(compData, N-1, nAtoms, true, mt);
        }
    }

   return compSims;
}

/* O(N) medoid calculation for n-ary objects.
 *
 * Parameters:
 *  - data: A feature array of shape (nSamples, nFeatures)
 *     --> Can also be a vector (N) of complementary similarities for each object. Useful when calculating medoid and outlier so you only calculate compSims once.
 *  - nAtoms: Number of atoms in the Molecular Dynamics (MD) system. nAtoms=1 for non-MD systems.
 *  - mt: The metric to use when calculating distance between n objects in an array. 
 * 
 * Returns: The index of the medoid in the dataset.
 * 
 * Reference: https://github.com/mqcomplab/MDANCE/blob/main/src/mdance/tools/bts.py#L241
*/
Index Cpptraj::Mdance::calculateMedoid(const ArrayXXd& data, int nAtoms, MD::Metric mt) {
    ArrayXd compSims = calculateCompSim(data, nAtoms, mt);
    return calculateMedoid(compSims);
}
Index Cpptraj::Mdance::calculateMedoid(const ArrayXd& data) {
    Index maxIdx;
    data.maxCoeff(&maxIdx);
    return maxIdx;
}

/* O(N) outlier calculation for n-ary objects.
 *
 * Parameters:
 *  - data: A feature array of shape (nSamples, nFeatures)
 *     --> Can also be a vector (N) of complementary similarities for each object. Useful when calculating medoid and outlier so you only calculate compSims once.
 *  - nAtoms: Number of atoms in the Molecular Dynamics (MD) system. nAtoms=1 for non-MD systems.
 *  - mt: The metric to use when calculating distance between n objects in an array. 
 * 
 * Returns: The index of the medoid in the dataset.
 * 
 * Reference: https://github.com/mqcomplab/MDANCE/blob/main/src/mdance/tools/bts.py#L271
*/
Index Cpptraj::Mdance::calculateOutlier(const ArrayXXd& data, int nAtoms, MD::Metric mt) {
    ArrayXd compSims = calculateCompSim(data, nAtoms, mt);
    return calculateOutlier(compSims);
}
Index Cpptraj::Mdance::calculateOutlier(const ArrayXd& data) {
    Index minIdx;
    data.minCoeff(&minIdx);
    return minIdx;
}

/* O(N * log(nTrimmed)) method of trimming a desired percentage of outliers (most dissimilar) from a feature array.
 *
 * Parameters:
 *  - data: A feature array of shape (nSamples, nFeatures).
 *  - nTrimmed: The desired # of outliers to be removed. Can be a number (int), or a percentage (float).
 *  - nAtoms: Number of atoms in the Molecular Dynamics (MD) system. nAtoms=1 for non-MD systems.
 *  - isMedoid: Criterion to use for data trimming
 *     --> if (!isMedoid): remove most dissimilar objects based on complement similarity.
 *     --> if (isMedoid): remove most dissimilar objects based on similarity to the medoid.
 * 
 *  Returns: A feature array with the desired outliers removed.
 * 
 *  Reference: https://github.com/mqcomplab/MDANCE/blob/main/src/mdance/tools/bts.py#L301
*/
ArrayXXd Cpptraj::Mdance::trimOutliers(const ArrayXXd& data, int nTrimmed, int nAtoms, bool isMedoid, MD::Metric mt) {
    Index N = data.rows();
    if (nTrimmed <= 0) {
        return data;
    }
    if (nTrimmed >= N) {
        return ArrayXXd(0, data.cols());
    }

    Index medoidIdx = -1;
    ArrayXd cSum;
    ArrayXd sqSumTotal;
    if (isMedoid) {
        medoidIdx = calculateMedoid(data, nAtoms, mt);
    } else {
        cSum = data.colwise().sum();
        sqSumTotal = data.square().colwise().sum();
    }

    ArrayXXd compData(2, data.cols());
    auto trimValue = [&](Index i) {
        if (isMedoid) {
            // Dissimilarity between frame i and the medoid (a 2-object
            // comparison; the medoid itself scores 0 and is never trimmed).
            compData.row(0) = data.row(i);
            compData.row(1) = data.row(medoidIdx);
            return extendedComparison(compData, 2, nAtoms, false, mt);
        }
        // Complementary similarity: dissimilarity of the data without frame i.
        compData.row(0) = cSum.transpose() - data.row(i);
        compData.row(1) = sqSumTotal.transpose() - data.row(i).square();
        return extendedComparison(compData, N - 1, nAtoms, true, mt);
    };

    // Select the nTrimmed outliers in O(N log nTrimmed) with a bounded heap.
    // comp_sim trims the LOWEST complementary similarities (a max-heap keeps
    // the nTrimmed smallest values); sim_to_medoid trims the HIGHEST
    // distances to the medoid (a min-heap keeps the nTrimmed largest).
    vector<pair<double, Index>> heap;
    heap.reserve(nTrimmed);
    auto cmp = [isMedoid](const pair<double, Index>& a, const pair<double, Index>& b) {
        return isMedoid ? a.first > b.first : a.first < b.first;
    };
    for (Index i = 0; i < N; ++i) {
        double val = trimValue(i);
        if ((Index)heap.size() < nTrimmed) {
            heap.emplace_back(val, i);
            std::push_heap(heap.begin(), heap.end(), cmp);
        } else if (isMedoid ? (val > heap.front().first) : (val < heap.front().first)) {
            std::pop_heap(heap.begin(), heap.end(), cmp);
            heap.back() = pair<double, Index>(val, i);
            std::push_heap(heap.begin(), heap.end(), cmp);
        }
    }

    vector<bool> trimmed(N, false);
    for (auto& p : heap) {
        trimmed[p.second] = true;
    }
    vector<Index> indices;
    indices.reserve(N - nTrimmed);
    for (Index i = 0; i < N; ++i) {
        if (!trimmed[i]) {
            indices.push_back(i);
        }
    }
    return data(indices, Eigen::all);
}
ArrayXXd Cpptraj::Mdance::trimOutliers(const ArrayXXd& data, float nTrimmed, int nAtoms, bool isMedoid, MD::Metric mt) {
    int num = std::floor(data.rows() * nTrimmed);
    if (num == 0)
        return data;
    return trimOutliers(data, num, nAtoms, isMedoid, mt);
}

/* O(N) method of selecting the most diverse subset of a data matrix using the complementary similarity.
 *
 * Parameters:
 *  - data: A feature array of shape (nSamples, nFeatures).
 *  - percentage: Indicates the percentage of data to be selected.
 *  - mt: The metric to use when calculating distance between n objects in an array
 *  - nAtoms: Number of atoms in the Molecular Dynamics (MD) system. nAtoms=1 for non-MD systems.
 *  - isCompSim: The method to use for diversity selection.
 *     --> if (!isCompSim): Uses stratified sampling.
 *     --> if (isCompSim): Maximizes the MSD between the selected objects and the rest of the data.
 *  - start: The initial seed for initiating diversity selection.
 *     --> You can also specify the seed indices as a vector<Index>
 * 
 * Returns: A vector of indices of the diversity selected data (in order selected).
 * 
 * Reference: https://github.com/mqcomplab/MDANCE/blob/016bd9aff30d1c2add26b36bfcf64aa665a34a1d/src/mdance/tools/bts.py#L376
*/
vector<Index> Cpptraj::Mdance::diversitySelection(const ArrayXXd& data, int percentage, MD::Metric mt, int nAtoms, bool isCompSim, MD::StartSeed start){
    if (isCompSim) {
        vector<Index> seed;
        switch(start) {
            case MD::StartSeed::Medoid: seed.push_back(calculateMedoid(data, nAtoms, mt)); break;
            case MD::StartSeed::Outlier: seed.push_back(calculateOutlier(data, nAtoms, mt)); break;
            case MD::StartSeed::Random: seed.emplace_back(rand() % data.rows()); break;
        }
        return diversitySelection(data, percentage, mt, nAtoms, seed);
    }
    Index N = data.rows();
    int nMax = N * percentage / 100;
    if (nMax > N) {
        throw std::runtime_error("Percentage is too high for the given matrix size");
    }
    vector<Index> indices (nMax);
    if (nMax == 1)
        indices[0] = 0;
    else {
        double step = (static_cast<double>(N - 1)) / (nMax -1);
        for (int i = 0; i < nMax; ++i){
            indices[i] = std::round(i * step);
        }
    } 
    ArrayXd compSims = calculateCompSim(data, nAtoms, mt);
    vector<pair<double,int>> compSimArray;
    compSimArray.reserve(compSims.size());
    for (int i=0; i<compSims.size(); ++i){
        compSimArray.emplace_back(-compSims[i],i);
    }
    std::sort(compSimArray.begin(), compSimArray.end());

    for(int i=0; i<nMax; ++i){
        indices[i] = compSimArray[indices[i]].second;
    }
    return indices;

}
vector<Index> Cpptraj::Mdance::diversitySelection(const ArrayXXd& data, int percentage, MD::Metric mt, int nAtoms, vector<Index>& indices){
    ArrayXXd selection = data(indices, Eigen::all);
    ArrayXXd selected (mt == MD::Metric::MSD ? 2 : 1, data.row(0).cols());
    selected.row(0) = selection.colwise().sum();
    if (mt == MD::Metric::MSD) {
        selected.row(1) = selection.square().colwise().sum();
    }
    
    int nTotal = data.rows();
    int nMax = nTotal * percentage / 100;

    set<Index> selectFromN;
    set<Index> selectedSet;
    for (size_t i=0; i<indices.size(); ++i) {
        selectedSet.insert(indices[i]);
    }
    for (int i=0; i<nTotal; ++i) {
        if (selectedSet.find(i) == selectedSet.end()){
            selectFromN.insert(i);
        }
    }

    indices.reserve(nMax);
    while ((int)indices.size() < nMax && !selectFromN.empty()) {
        Index newIndexN = getNewIndexN(data, mt, selected, indices.size(), selectFromN, nAtoms);

        selected.row(0) += data.row(newIndexN);
        if (mt == MD::Metric::MSD)
            selected.row(1) = selected.row(1) + data.row(newIndexN).square();

        selectFromN.erase(newIndexN);
        indices.push_back(newIndexN);
    }

    return indices;
}

/* Extract the new index to add to the list of selected indices.
 *
 * Parameters:
 *  - data: A feature array of shape (nSamples, nFeatures).
 *  - mt: The metric to use when calculating distance between n objects in an array
 *  - selectedCondensed: A fingerprint feature array that can take on multiple shapes:
 *     --> if (mt == MSD): a ArrayXXd with 2 rows (cSum, sqSum)
 *     --> else: a ArrayXXd with 1 row (cSum)
 *  - N: number of selected objects
 *  - selectFromN: Array of indices to select from
 *  - nAtoms: Number of atoms in the Molecular Dynamics (MD) system. nAtoms=1 for non-MD systems.
 * 
 * Returns: index of the new fingerprint to add to the selected indices.
 * 
 * Reference: https://github.com/mqcomplab/MDANCE/blob/016bd9aff30d1c2add26b36bfcf64aa665a34a1d/src/mdance/tools/bts.py#L489
*/
Index Cpptraj::Mdance::getNewIndexN(const ArrayXXd& data, MD::Metric mt, ArrayXXd& selectedCondensed, int N, set<Index>& selectFromN, int nAtoms) {
    // Number of fingerprints already selected and the new one to add
    int nTotal = N + 1;

    double maxVal = -1;
    Index idx = data.row(0).size();

    ArrayXXd temp = selectedCondensed;

    for (auto i=selectFromN.begin(); i!=selectFromN.end(); ++i){
        temp.row(0) = selectedCondensed.row(0) + data.row(*i);
        double simIdx;
        if (mt == MD::Metric::MSD) {
            temp.row(1) = selectedCondensed.row(1) + data.row(*i).square();
            simIdx = extendedComparison(temp, nTotal, nAtoms, true, mt);
        } else {
            simIdx = extendedComparison(temp, nTotal, nAtoms, true, mt);
        }

        if (simIdx > maxVal){
            maxVal = simIdx;
            idx = *i;
        }
    }
    return idx;
}

/* Representative sampling according to compSim values: Divides the range of comp_sim values in nbins and then uniformly selects n_samples molecules, consecutively taking one from each bin
 *
 * Parameters:
 *  - data: A feature array of shape (nSamples, nFeatures).
 *  - mt: The metric to use when calculating distance between n objects in an array
 *  - nAtoms: Number of atoms in the Molecular Dynamics (MD) system. nAtoms=1 for non-MD systems.
 *  - nBins: Number of bins to divide the compSim values.
 *  - nSamples: Number of samples to be selected.
 *  - hardCap: whether the number of samples will be *exactly* nSamples
 * 
 * Returns: List of indices of the sampled objects in the original data
 * 
 * Reference: https://github.com/mqcomplab/MDANCE/blob/016bd9aff30d1c2add26b36bfcf64aa665a34a1d/src/mdance/tools/bts.py#L652
*/
ArrayXi Cpptraj::Mdance::repSample(const ArrayXXd& data, MD::Metric mt, int nAtoms, int nBins, double nSamples, bool hardCap) {
    // nSamples in (0,1) is a fraction of the data; anything >= 1 is a count.
    if (nSamples < 1) {
        return repSample(data, mt, nAtoms, nBins, (int)(data.rows() * nSamples), hardCap);
    }
    return repSample(data, mt, nAtoms, nBins, (int)nSamples, hardCap);
}
ArrayXi Cpptraj::Mdance::repSample(const ArrayXXd& data, MD::Metric mt, int nAtoms, int nBins, int nSamples, bool hardCap) {
    Index N = data.rows();
    if (nSamples <= 0 || N == 0) {
        return ArrayXi();
    }
    if (nBins < 1) {
        throw std::invalid_argument("nBins must be at least 1");
    }

    ArrayXd compSims = calculateCompSim(data, nAtoms, mt);
    // (compSim value, original index), sorted by value
    vector<pair<double,int>> compSimArray;
    compSimArray.reserve(N);
    for (int i=0; i<N; ++i){
        compSimArray.emplace_back(compSims[i], i);
    }
    std::sort(compSimArray.begin(), compSimArray.end());

    double mi = compSimArray.front().first;
    double ma = compSimArray.back().first;
    double binWidth = (ma - mi) / nBins;

    // Distribute the objects across nBins bins of equal compSim width. When
    // every value is identical the reference puts everything in the last bin.
    vector<vector<int>> bins(nBins);
    for (auto& p : compSimArray) {
        int b = binWidth > 0 ? (int)((p.first - mi) / binWidth) : nBins - 1;
        if (b < 0) b = 0;
        if (b >= nBins) b = nBins - 1;  // the maximum value lands in the last bin
        bins[b].push_back(p.second);
    }

    // Round-robin over the bins, taking the next unused object of each bin per
    // pass. With hardCap the collection stops exactly at nSamples; otherwise
    // the current pass is completed first (so slightly more than nSamples may
    // be returned), as in the Python reference.
    vector<int> sampled;
    sampled.reserve(nSamples);
    for (int depth = 0; (int)sampled.size() < nSamples; ++depth) {
        bool added = false;
        for (int b = 0; b < nBins; ++b) {
            if ((int)bins[b].size() > depth) {
                sampled.push_back(bins[b][depth]);
                added = true;
                if (hardCap && (int)sampled.size() >= nSamples)
                    break;
            }
        }
        if (!added) {
            break;  // every bin is exhausted: fewer than nSamples objects exist
        }
    }

    return Eigen::Map<ArrayXi>(sampled.data(), sampled.size());
}

/* Refine a distance matrix by setting the diagonal to zero and symmetrizing the matrix
 *
 * Parameters:
 *  - data: A feature array of shape (nSamples, nFeatures).
 * 
 * Returns: A refined 2D matrix.
 * 
 * Reference: https://github.com/mqcomplab/MDANCE/blob/016bd9aff30d1c2add26b36bfcf64aa665a34a1d/src/mdance/tools/bts.py#L720
*/
ArrayXXd Cpptraj::Mdance::refineDisMatrix(const ArrayXXd& data) {
    if (data.rows() == 1 || data.cols() == 1) {
        throw std::invalid_argument("Matrix must be 2D.");
    }
    if (data.rows() != data.cols()) {
        throw std::invalid_argument("Matrix must be square.");
    }

    ArrayXXd distances = (data + data.transpose()) / 2;
    distances = distances - distances.minCoeff();
    for (int i=0; i<distances.rows(); ++i) {
        distances.row(i)[i] = 0;
    }
    return distances;
}


//not finished
ArrayXXd Cpptraj::Mdance::alignTraj(const ArrayXXd& data, int nAtoms, MD::AlignMethod alignMeth){
    /*
        Aligns trajectory using uniforms or kronecker alignment

        Parameters
        --------------
        data: dimensions: nSamples x nFeatures
        nAtoms: number of atoms in the system
        alignMethod: {Uni, Kron}, default = None

        Returns
        --------------
        matrix of aligned data

        References
        --------------
        Klem, H., Hocky, G. M., and McCullagh M., `"Size-and-Shape Space Gaussian 
        Mixture Models for Structural Clustering of Molecular Dynamics Trajectories"`_.
        *Journal of Chemical Theory and Computation* **2022** 18 (5), 3218-3230

        .. _"Size-and-Shape Space Gaussian Mixture Models for Structural Clustering of Molecular Dynamics Trajectories":
        https://pubs.acs.org/doi/abs/10.1021/acs.jctc.1c01290
    */

    (void)nAtoms;  // stub: unused until the alignment methods are implemented

    if(alignMeth == MD::AlignMethod::None){
        return data;
    }

    return data;
}
#endif /* HAS_EIGEN */
