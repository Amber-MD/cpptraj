#ifndef INC_HUNGARIAN_H
#define INC_HUNGARIAN_H
#include "DataSet_MatrixDbl.h"
/// Use Hungarian algorithm to perform minimum-matching on a matrix.
class Hungarian {
  public:
    Hungarian() {}
    Hungarian(DataSet_MatrixDbl const&);
    int Initialize(size_t);
    void AddElement(double d) { matrix_.AddElement( d ); }
    /// \return Array containing Map[col] = row
    std::vector<int> Optimize();
  private:
    int AssignRowsToColumns();
    int CoverZeroElements();
    void UpdateMatrix();
#   ifdef DEBUG_HUNGARIAN
    void PrintLines(const char*);
    void PrintMatrix(const char*);
#   endif
    DataSet_MatrixDbl matrix_;         ///< Working matrix. 
    std::vector<bool> lineThroughRow_; ///< True if specified row is marked.
    std::vector<bool> lineThroughCol_; ///< True if specified col is marked.
    std::vector<int> assignRowToCol_;  ///< map[col] = row
    std::vector<int> assignColToRow_;  ///< map[row] = col
    // TODO: Make these size_t
    int nrows_;                        ///< # of rows in matrix.
    int ncols_;                        ///< # of cols in matrix.
};
#endif
