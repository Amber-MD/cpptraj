#ifndef INC_DATASET_2D_H
#define INC_DATASET_2D_H
class DataSet_2D {
  public:
    /// Set up matrix for given # rows and columns.
    virtual int Allocate2D(size_t, size_t) = 0;
    // TODO: Add AllocateHalf and AllocateTri
    /// Write 2D data to file (2D)
    virtual void Write2D(CpptrajFile&,int,int) const = 0;
    /// \return Data from matrix at row/col
    virtual double GetElement(size_t, size_t) const = 0;
    /// \return the number of rows.
    virtual size_t Nrows() const = 0;
    /// \return the number of columns.
    virtual size_t Ncols() const = 0;
};
#endif
