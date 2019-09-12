#ifndef INC_DATASET_1D_H
#define INC_DATASET_1D_H
#include "DataSet.h"
#include "CpptrajFile.h"
/// Class that all 1D scalar DataSets will inherit.
class DataSet_1D : public DataSet {
  public:
    DataSet_1D() {}
    DataSet_1D(DataSet::DataType tIn, TextFormat const& fIn) : DataSet(tIn, SCALAR_1D, fIn, 1) {}
    virtual ~DataSet_1D() {}
    // ---- DataSet_1D functions -----------------
    /// \return data from set at position as double precision.
    virtual double Dval(size_t) const = 0;
    /// \return the value of the X coordinate at position.
    virtual double Xcrd(size_t) const = 0;
    /// \return Memory address at position cast to void *.
    virtual const void* VoidPtr(size_t) const = 0;
    // -------------------------------------------
    /// \return Average over set Y values
    double Avg()           const { return Avg( 0 ); }
    /// \return Average over set Y values; calculate standard deviation.
    double Avg(double& sd) const { return Avg(&sd); }
    /// \return Set minimum Y value.
    double Min() const;
    /// \return Set maximum Y value.
    double Max() const;
    /// Calculate cross-correlation to another set.
    int CrossCorr(DataSet_1D const&, DataSet_1D&, int, bool, bool) const;
    /// Calculate auto-correlation
    double CorrCoeff(DataSet_1D const&) const;
    /// Calculate linear regression; report slope, intercept, and correlation.
    int LinearRegression(double&, double&, double&, CpptrajFile*) const;
    /// Integration types.
    enum IntegrationType { TRAPEZOID = 0 };
    /// \return sum of integration over DataSet.
    double Integrate(IntegrationType) const;
    /// \return sum of integration over DataSet; compute cumulative sum.
    double Integrate(IntegrationType, std::vector<double>&) const;
  private:
    double Avg(double*) const;
};
#endif
