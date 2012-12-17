#ifndef INC_TRIANGLEMATRIX_H
#define INC_TRIANGLEMATRIX_H
#include <cstddef> // size_t
#include "DataSet.h"
// Class: TriangleMatrix
/// Store the upper half of a symmetric matrix
/** Useful when calculating e.g. all N^2 distances between all atoms, 
  * where the diagonal elements would be zero and element i,j == element j,i.
  * Accepts doubles, but internal storage is float to reduce memory footprint
  */
class TriangleMatrix : public DataSet {
  public :
    int Nrows()            { return nrows_;     }
    int Nelements()        { return (int)nelements_; }

    TriangleMatrix();
    TriangleMatrix(int);
    TriangleMatrix(const TriangleMatrix&);
    ~TriangleMatrix();
    TriangleMatrix & operator=(const TriangleMatrix &);

    int SaveFile(const char*);
    int LoadFile(const char*,int);
    int Setup(int);
    void Ignore(int);
    int AddElement(double);
    int AddElement(float);
    void SetElement(int,int,double);
    void SetElementF(int,int,float);
    double GetElement(int,int) const;
    float GetElementF(int,int) const;
    double FindMin(int *, int *);
    void PrintElements();

    double operator[](int idx) { return (double)elements_[idx]; }

    int Xmax() { return nrows_ - 1; }
    int Size() { return (int)nelements_; }
    void Write2D( CpptrajFile&, int, int);
    void GetDimensions(std::vector<int>&); 
  private:
    float *elements_;       ///< Hold all elements
    int nrows_;             ///< Number of elements in one row
    size_t nelements_;      ///< Total number of elements
    size_t currentElement_; ///< Current element, used in AddElement only.
    bool *ignore_;          ///< If true, ignore the row/col when printing/searching etc

    int calcIndex(int,int) const;
};
#endif
