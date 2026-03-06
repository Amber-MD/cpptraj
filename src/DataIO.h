#ifndef INC_DATAIO_H
#define INC_DATAIO_H
#include "ArgList.h"
#include "DataSetList.h"
#include "CpptrajFile.h"
#include "BaseIOtype.h"
/// Base class that all DataIO objects inherit from.
class DataIO : public BaseIOtype {
  public:
    DataIO();
    /// Valid for 1d, 2d, and/or 3d data sets.
    DataIO(bool, bool, bool);
    virtual ~DataIO() {}
    // ----- Inherited Functions -----------------
    virtual int processReadArgs(ArgList&) = 0;
    virtual int ReadData(FileName const&, DataSetList&, std::string const&) = 0;
    virtual int processWriteArgs(ArgList&) = 0;
    virtual int WriteData(FileName const&, DataSetList const&) = 0;
    virtual bool ID_DataFormat(CpptrajFile&) = 0; // TODO: -> BaseIOtype?
    /// \return True if this DataIO valid for given DataSet
    bool CheckValidFor(DataSet const&) const;
    /// Set DataIO debug level.
    void SetDebug(int d) { debug_ = d; }
    /// Set x column format.
    void SetXcolFmt(TextFormat::FmtType t) { xcol_fmt_ = t; }
    /// Set x column width, and precision.
    void SetXcolPrec(int w, int p) { xcol_width_ = w; xcol_prec_ = p; x_prec_set_ = true; }
    /// \return Current x column format
    TextFormat::FmtType XcolFmt() const { return xcol_fmt_; }
    /// \return Current x column width
    int XcolWidth()               const { return xcol_width_; }
    /// \return Current x column precision
    int XcolPrec()                const { return xcol_prec_;  }

    typedef std::vector<DataSet*>::const_iterator set_iterator;
    /// \return iterator to beginning of sets added by last ReadData() call
    set_iterator added_begin() const { return mySets_.begin(); }
    /// \return iterator to end of sets added by last ReadData() call
    set_iterator added_end()   const { return mySets_.end(); }
    /// \return the last set added by last ReadData() call
    DataSet* added_back()      const { return mySets_.back(); }
    /// \return Number of sets added by last ReadData() call
    unsigned int Nadded()      const { return mySets_.size(); }
  protected:
    /// Note a set that has been added to a DataSetList by this DataIO
    void AddedByMe(DataSet*);
    /// Clear list of sets added by this DataIO
    void ClearAddedByMe();
    /// Indicate this DataIO is valid for given DataSet type
    void SetValid(DataSet::DataType t) { valid_.push_back( t ); }
    /// Check that all sets in given list have given dimension.
    static int CheckAllDims(DataSetList const&, unsigned int);
    /// \return True if the X coordinates for both sets match exactly
    static bool xDimMatch(DataSet const&, DataSet const&);
    /// Typedef for an array of DataSetLists
    typedef std::vector<DataSetList> DSLarray;
    /// Check that X dim for all sets in given list match; all must be 1D.
    static DSLarray CheckXDimension(DataSetList const&);
    /// \return max size of DataSets in given list.
    static size_t DetermineMax(DataSetList const&);
    /// Convert flattened matrix array to matrix in DataSetList.
    static DataSet* DetermineMatrixType(std::vector<double> const&, int, int,
                                        DataSetList&, std::string const&);
    /// \return true if x column format/width/precision previously set TODO always use these values?
    bool XcolPrecSet()            const { return x_prec_set_; }
    int debug_;
  private:
    std::vector<DataSet::DataType> valid_; ///< Data sets for which DataIO is valid writer.
    std::vector<DataSet*> mySets_; ///< Can be used by child class to note sets that have been added
    TextFormat::FmtType xcol_fmt_; ///< X column format type
    int xcol_width_;               ///< X column width
    int xcol_prec_;                ///< X column precision
    bool x_prec_set_;              ///< True if X column width/precision have been explicitly set.
    bool valid1d_; ///< Valid for all 1D data sets. //TODO Remove
    bool valid2d_; ///< Valid for all 2D data sets.
    bool valid3d_; ///< Valid for all 3D data sets.
};
#endif
