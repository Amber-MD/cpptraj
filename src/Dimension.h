#ifndef INC_DIMENSION_H
#define INC_DIMENSION_H
#include <string>
/// Holds information about a histogram dimension
class Dimension {
  public:
    Dimension();
    Dimension(const Dimension&);
    Dimension& operator=(const Dimension&);
 
    void SetLabel(std::string const&);
    void SetMin(double);
    void SetMax(double);
    void SetStep(double);
    void SetBins(int);
    void SetOffset(int);

    const char* c_str() { return label_.c_str(); }
    std::string const& Label() { return label_; }
    double Min() { return min_; }
    double Max() { return max_; }
    double Step() { return step_; }
    int Bins() { return bins_; }
    int Offset() { return offset_; }

    int CalcBinsOrStep();
    void PrintDim();
    
  private:
    std::string label_;
    double min_;
    double max_;
    double step_;
    int bins_;
    int offset_;
};
#endif
