#ifndef INC_TRIANGLEMATRIX_H
#define INC_TRIANGLEMATRIX_H
/// Class: TriangleMatrix
/// Store the upper half of a symmetric matrix, useful when calculating
/// e.g. all N^2 distances between all atoms, where the diagonal elements
/// would be zero and element i,j == element j,i
class TriangleMatrix {
    double *elements; // Hold all elements
    int Nrows;        // Number of elements in one row
    int Nelements;    // Total number of elements
    int currentElement;
  public :
    TriangleMatrix();
    ~TriangleMatrix();

    int Setup(int);
    int AddElement(double);
    double GetElement(int,int);
};
#endif
