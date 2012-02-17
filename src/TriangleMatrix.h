#ifndef INC_TRIANGLEMATRIX_H
#define INC_TRIANGLEMATRIX_H
// Class: TriangleMatrix
/// Store the upper half of a symmetric matrix
/** Useful when calculating e.g. all N^2 distances between all atoms, 
  * where the diagonal elements would be zero and element i,j == element j,i.
  * Accepts doubles, but internal storage is float to reduce memory footprint
  */
class TriangleMatrix {
    float *elements;       ///< Hold all elements
    int nrows;             ///< Number of elements in one row
    size_t nelements;      ///< Total number of elements
    size_t currentElement; ///< Current element, used in AddElement only.
    bool *ignore;          ///< If true, ignore the row/col when printing/searching etc

    int calcIndex(int,int);
  public :
    int Nrows()            { return nrows;     }
    int Nelements()        { return nelements; }

    TriangleMatrix();
    TriangleMatrix(const TriangleMatrix&);
    ~TriangleMatrix();

    TriangleMatrix & operator=(const TriangleMatrix &);
    int SaveFile(char*);
    int LoadFile(char*,int);
    int Setup(int);
    void Ignore(int);
    int AddElement(double);
    int AddElement(float);
    void SetElement(int,int,double);
    void SetElementF(int,int,float);
    double GetElement(int,int);
    float GetElementF(int,int);
    double FindMin(int *, int *);
    void PrintElements();
};
#endif
