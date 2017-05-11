#ifndef INC_ANALYSIS_TI_H
#define INC_ANALYSIS_TI_H
#include "Analysis.h"
#include "Array1D.h"
class Analysis_TI : public Analysis {
  public:
    Analysis_TI();
    DispatchObject* Alloc() const { return (DispatchObject*)new Analysis_TI(); }
    void Help() const;

    Analysis::RetType Setup(ArgList&, AnalysisSetup&, int);
    Analysis::RetType Analyze();
  private:
    /// Averaging type: normal, skip #s of points, incremental avg
    enum AvgType { AVG = 0, SKIP, INCREMENT, BOOTSTRAP };
    /// Integration type
    enum ModeType { GAUSSIAN_QUAD = 0, TRAPEZOID };
    typedef std::vector<int> Iarray;
    typedef std::vector<double> Darray;
    typedef std::vector<DataSet*> DSarray;

    int SetQuadAndWeights(int);
    void Integrate_Trapezoid(Darray&) const;
    int Calc_Bootstrap();
    int Calc_Nskip();
    int Calc_Increment();
    int Calc_Avg();

    Array1D input_dsets_; ///< Input DV/DL data sets
    Iarray nskip_;        ///< Numbers of data points to skip in calculating <DV/DL>
    DataSet* dAout_;      ///< Free energy data set
    DataSet* dA_SD_;      ///< Standard deviation of free energy
    DataFile* curveout_;  ///< File to write DV/DL curves to.
    DataSetList* masterDSL_;
    DSarray curve_;       ///< TI curve data set for each skip/increment value
    Darray xval_;         ///< Hold abscissas corresponding to data sets.
    Darray wgt_;          ///< Hold Gaussian quadrature weights
    double bootstrap_fac_;///< Fraction of total points to use when # pts not specified.
    ModeType mode_;       ///< Integration mode
    AvgType avgType_;     ///< Type of averaging to be performed.
    int debug_;
    int n_bootstrap_pts_; ///< # points for bootstrap error analysis
    int n_bootstrap_samples_; ///< # of times to resample for bootstrap analysis
    int bootstrap_seed_;  ///< # RNG seed. Input data set index added for each set.
    int avg_increment_;   ///< # points to skip between each average calc
    int avg_max_;         ///< Max number of points to use in average (default all) 
    int avg_skip_;        ///< Number of points to skip when calculating the average.
};
#endif
