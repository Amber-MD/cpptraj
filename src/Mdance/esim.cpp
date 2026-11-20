#include "esim.h"
#ifdef HAS_EIGEN
using namespace Cpptraj::Mdance;
/* Helper functions for calculateCounters.
 * Reference: https://github.com/mqcomplab/MDANCE/blob/main/src/mdance/tools/esim.py#L65
*/
double fs(double d, int wFactor, int nObjects) {
    if (wFactor)
        return pow(wFactor, d - nObjects);
    return d/nObjects;
}
double fd(double d, int wFactor, int nObjects) {
    if (wFactor)
        return pow(wFactor, nObjects % 2 - d);
    return 1-(d-nObjects % 2)/nObjects;
}

/* Calculate 1-similariy (a), 0-similarity (d), and dissimilarity (dis) counters
 *
 * Parameters:
 *  - cTotal: Vector (nFeatures) containing the sums of each column of the fingerprint matrix.
 *  - nObjects: Number of objects to be compared.
 *  - cThreshold: Coincidence threshold
 *  - wFactor: Type of weight function that will be used
 *     --> 0 = "fraction": similarity = d[k]/n, 
 *                      dissimilarity = 1 - (d[k] - nObjects % 2)/nObjects
 *     --> n = "power_n" : similarity = n**(d[k] - nObjects)
 *                      dissimilarity = n**(nObjects % 2 - d[k])
 * 
 * Returns: Dictionary with the weighted and non-weighted counters.
*/
MD::Counters Cpptraj::Mdance::calculateCounters(const ArrayXd& cTotal, int nObjects, MD::Threshold& cThreshold, int wFactor){
    int a=0;
    int d=0;
    int dis=0;
    double wa=0;
    double wd=0;
    double wDis=0;

    for(double w : cTotal){
        double diff = 2 * w - nObjects;
        if(diff > cThreshold.value){
            ++a;
            wa += fs(diff, wFactor, nObjects);
        } 
        else if (-diff > cThreshold.value){
            ++d;
            wd += fs(-diff, wFactor, nObjects);
        }
        else{
            ++dis;
            wDis += fd(abs(diff), wFactor, nObjects);
        }
    }

    return MD::Counters(a, wa, d, wd, dis, wDis);
}

/* Generate a dictionary with the similarity indices.
 *
 * Parameters:
 *  - cTotal: Vector (nFeatures) containing the sums of each column of the fingerprint matrix.
 *  - nObjects: Number of objects to be compared.
 *  - cThreshold: Coincidence threshold
 *  - wFactor: Type of weight function that will be used
 *     --> 0 = "fraction": similarity = d[k]/n, 
 *                      dissimilarity = 1 - (d[k] - nObjects % 2)/nObjects
 *     --> n = "power_n" : similarity = n**(d[k] - nObjects)
 *                      dissimilarity = n**(nObjects % 2 - d[k])
 * 
 * Returns: Dictionary with the similarity indices.
 * 
 * Reference: https://github.com/mqcomplab/MDANCE/blob/main/src/mdance/tools/esim.py#L122
 * 
 * TODO: implement other indices https://github.com/mqcomplab/MDANCE/blob/main/src/mdance/tools/esim.py#L253
*/
MD::Indices Cpptraj::Mdance::genSimIdx(const ArrayXd& cTotal, int nObjects, MD::Threshold& cThreshold, int wt) {
    MD::Counters cnt = Cpptraj::Mdance::calculateCounters(cTotal, nObjects, cThreshold, wt);
    double bub = (sqrt(cnt.wa * cnt.wd) + cnt.wa) / (sqrt(cnt.a * cnt.d) + cnt.a + cnt.totalDis);
    double fai = (cnt.wa + 0.5 * cnt.wd) / cnt.p;
    double gle = 2 * cnt.wa / (2 * cnt.a + cnt.totalDis);
    double ja = 3 * cnt.wa / (3 * cnt.a + cnt.totalDis);
    double jt = cnt.wa / (cnt.a + cnt.totalDis);
    double rt = cnt.totalWsim / (cnt.p + cnt.totalDis);
    double rr = cnt.wa / cnt.p;
    double sm = cnt.totalWsim / cnt.p;
    double ss1 = cnt.wa / (cnt.a + 2 * cnt.totalDis);
    double ss2 = 2 * cnt.totalWsim / (cnt.p + cnt.totalSim);
    return MD::Indices(bub, fai, gle, ja, jt, rt, rr, sm , ss1, ss2);
}
#endif /* HAS_EIGEN */
