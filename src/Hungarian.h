#ifndef INC_HUNGARIAN_H
#define INC_HUNGARIAN_H
#include "Matrix_2D.h"
/// Use Hungarian algorithm to perform minimum-matching on a matrix.
class Hungarian {
  public:
    Hungarian() {}
    Hungarian(Matrix_2D const&);
    /// \return Array containing Map[ref] = tgt
    std::vector<int> Optimize();
  private:
    int Assign();
    int CoverZeroElements();
    void UpdateMatrix();
    void PrintLines(const char*);

    Matrix_2D matrix_;
    std::vector<bool> lineThroughRow_;
    std::vector<bool> lineThroughCol_;
    std::vector<int> assignRowToCol_;
    std::vector<int> assignColToRow_;
};
#endif
