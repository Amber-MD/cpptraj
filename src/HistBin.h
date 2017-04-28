#ifndef INC_HISTBIN_H
#define INC_HISTBIN_H
#include "Dimension.h"
/// Holds information for histogram binning in a dimension.
class HistBin : public Dimension {
  public:
    HistBin();
    /// CONSTRUCTOR - Bins, Min, Step, Label
    HistBin(int, double, double, std::string const&);
    HistBin(const HistBin&);
    HistBin& operator=(const HistBin&);
 
    double Max()               const { return max_;      }
    int Bins()                 const { return bins_;     }
    /// Attempt to set up bins/step from given min max and either step/bins.
    int CalcBinsOrStep(double,double,double,int,std::string const&);
    void PrintHistBin() const;
  private:
    double max_;
    int bins_;
};
#endif
