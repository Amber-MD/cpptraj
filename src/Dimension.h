#ifndef INC_DIMENSION_H
#define INC_DIMENSION_H
#include <string>
#include <vector>
/// Holds information about a coordinate dimension.
class Dimension {
  public:
    enum DimIdxType { X = 0, Y, Z };
    /// DEFAULT CONSTRUCTOR
    Dimension() : min_(0.0), step_(0.0) {}
    /// CONSTRUCTOR - Min, step
    Dimension(double m, double s) : min_(m), step_(s) {}
    /// CONSTRUCTOR - Min, step, label
    Dimension(double m, double s, std::string const& l) :
      label_(l), min_(m), step_(s) {}
    /// Set dimension from given values. \return true if monotonic
    inline bool SetDimension(std::vector<double> const&, std::string const&);
    /// COPY CONSTRUCTOR
    Dimension(const Dimension& rhs) :
      label_(rhs.label_), min_(rhs.min_), step_(rhs.step_) {}
    /// ASSIGNMENT OP
    Dimension& operator=(const Dimension& rhs) {
      if (this != &rhs) {
        label_ = rhs.label_;
        min_ = rhs.min_;
        step_ = rhs.step_;
      }
      return *this;
    }
    /// \return true if this dim min/step not equal to given dim min/step.
    bool operator!=(const Dimension& rhs) const {
      if (min_ != rhs.min_ || step_ != rhs.step_) return true;
      return false;
    }
    /// \return true if this dim min/step equal to given min/step.
    bool operator==(const Dimension& rhs) const {
      if (min_ != rhs.min_ || step_ != rhs.step_) return false;
      return true;
    }
    /// Set dimension with given min, step, and label.
    void SetDimension(double m, double s, std::string const& l) {
      label_ = l;
      min_ = m;
      step_ = s;
    }
    /// Only set dimension label FIXME is this needed? Useful for DataSet_Mesh
    void SetLabel(std::string const& l) { label_ = l; }
    void ChangeStep(double s) { step_ = s; } // Used by DataFile
    void ChangeMin(double m) { min_ = m; } // Used by DataFile

    std::string const& Label() const { return label_;         }
    const char* label()        const { return label_.c_str(); }
    double Min()               const { return min_;           }
    double Step()              const { return step_;          }
    double Coord(size_t i)     const { return ((step_ * (double)i) + min_); }
  private:
    std::string label_;
    double min_;
    double step_;
};
// ----- INLINE FUNCTIONS ------------------------
bool Dimension::SetDimension(std::vector<double> const& Vals, std::string const& labelIn)
{
  bool isMonotonic = true;
  double step = 1.0;
  if (Vals.size() > 1) {
    step = (Vals.back() - Vals.front()) / (double)(Vals.size() - 1);
    for (std::vector<double>::const_iterator X = Vals.begin()+2; X != Vals.end(); ++X)
      if ((*X - *(X-1)) - step > 0.00000000000001) {
        isMonotonic = false;
        break;
      }
    // Set dim even for non-monotonic sets so Label is correct. FIXME is this ok?
    SetDimension( Vals.front(), step, labelIn );
  } else
    // No X values. set generic X dim.
    SetDimension(1.0, 1.0, labelIn);
  return isMonotonic;
}
#endif
