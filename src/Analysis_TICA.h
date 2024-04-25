#ifndef INC_ANALYSIS_TICA_H
#define INC_ANALYSIS_TICA_H
#include "Analysis.h"
#include "Array1D.h"
class DataSet_2D;
class DataSet_Modes;
/// Perform time-independent correlation analysis 
class Analysis_TICA : public Analysis {
  public:
    Analysis_TICA();
    ~Analysis_TICA();
    DispatchObject* Alloc() const { return (DispatchObject*)new Analysis_TICA(); }
    void Help() const;

    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
    Analysis::RetType Analyze();
  private:
    enum EvectorScaleType {
      NO_SCALING = 0, ///< Do not scale eigenvectors
      KINETIC_MAP,    ///< Scale eigenvectors by eigenvalues
      COMMUTE_MAP     ///< Scale eigenvectors by regularized time scales
    };
    enum CalcType {
      COORDS = 0, ///< Use XYZ coordinates
      DATA,       ///< Use 1D data sets
      PERIODIC    ///< Use 1D periodic data sets
    };

    typedef std::vector<DataSet_1D*> DSarray;
    typedef std::vector<double> Darray;

    /// Calculate total weight
    static double calc_total_weight(Darray const&, unsigned int);
    /// Calculate sums of X and Y for 1D data sets in a single pass
    static void calc_sums_from1Dsets(Darray&, std::vector<DataSet_1D*> const&, unsigned int, unsigned int);
    /// Calculate sums of X and Y for COORDS sets in a single pass
    void calc_sums_fromCoordsSet(Darray&, unsigned int) const;
    /// Create XXYY and XYYX matrices from an array of 1D sets
    static void create_matrices_from1Dsets(std::vector<DataSet_1D*> const&, DataSet_2D*, DataSet_2D*,
                                           Darray const&, unsigned int, unsigned int);
    /// Create XXYY and XYYX matrices from an array of 1D sets in a single pass
    void create_matrices_from1Dsets(DataSet_2D*, DataSet_2D*, Darray const&, unsigned int) const;
    /// Create XXYY and XYYX matrices from an array of periodic 1D sets in a single pass
    void create_matrices_fromPeriodicSets(DataSet_2D*, DataSet_2D*, Darray const&, unsigned int) const;
    /// Create XXYY and XYYX matrices from COORDS set in a single pass
    void create_matrices_fromCoordsSet(DataSet_2D*, DataSet_2D*, Darray const&, unsigned int) const;
    /// Calculate TICA matrices
    int calcMatrices(unsigned int) const;
    /// Calculate TICA modes
    int calculateTICA(Darray const&, DataSet_2D const&, DataSet_2D const&) const;
#   ifdef CPPTRAJ_DEBUG_TICA
    /// Analyze using coordinates data set (TgtTraj_)
    Analysis::RetType analyze_crdset();
    /// Analyze using 1D data sets (sets_)
    Analysis::RetType analyze_datasets();
    /// Calculate instantaneous and lagged covariance matrices
    int calculateCovariance_C0CT(DSarray const&) const;
#   endif

    Array1D sets_;                   ///< Input 1D data sets (data)
    DSarray cossin_;                 ///< Input 1D periodic sets converted to cos/sin pairs (data)
    DataSet_Coords* TgtTraj_;        ///< Input trajectory (crdset)
    AtomMask mask1_;                 ///< Atoms to use in matrix calc (crdset)
    int lag_;                        ///< TICA time lag
    int debug_;
    CalcType calcType_;              ///< Type of calculation (TODO allow mixed data types)
    bool useMass_;                   ///< Control whether to mass-weight TODO enable
    EvectorScaleType evectorScale_;  ///< Eigenvector scaling type
    DataSet_Modes* ticaModes_;       ///< Output TICA modes
    DataSet_1D* cumulativeVariance_; ///< Hold output cumulative variance
#   ifdef CPPTRAJ_DEBUG_TICA
    AtomMask mask2_;                 ///< Second atom mask for debugging full covar matrix
    CpptrajFile* debugC0_;           ///< Debug output for C0
    CpptrajFile* debugCT_;           ///< Debug output for CT
#   endif
};
#endif
