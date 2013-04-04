#ifndef INC_HUNGARIAN_H
#define INC_HUNGARIAN_H
#include "Matrix_2D.h"
/// Use Hungarian algorithm to perform minimum-matching on a matrix.
class Hungarian {
  public:
    Hungarian() {}
    Hungarian(Matrix_2D const&);
    /// \return Array containing Map[col] = row
    std::vector<int> Optimize();
  private:
    int Assign();
    int CoverZeroElements();
    void UpdateMatrix();
#   ifdef DEBUG_HUNGARIAN
    void PrintLines(const char*);
    void PrintMatrix(const char*);
#   endif
    Matrix_2D matrix_;                 // Working matrix. 
    std::vector<bool> lineThroughRow_; // True if specified row is marked.
    std::vector<bool> lineThroughCol_; // True if specified col is marked.
    std::vector<int> assignRowToCol_;  // map[col] = row
    std::vector<int> assignColToRow_;  // map[row] = col
};
#endif
