#ifndef INC_DATAFILTER_H
#define INC_DATAFILTER_H
#include <vector>
#include <cstddef> // size_t
class ArgList;
class DataSetList;
class DataFileList;
class DataSet;
class DataSet_1D;
class DataSet_integer;
/// Can be used to count/filter out data set elements based on user-defined criteria.
class DataFilter {
    enum ResultType { PASSED = 0, FILTERED };
  public:
    DataFilter();
    /// \return Keywords recognized by InitFilter();
    static const char* Keywords();
    /// Process arguments, get sets to filter
    int InitFilter(ArgList&, DataSetList&, DataFileList&, int);
    /// \return 1 if specified index was filtered, 0 otherwise
    int FilterIndex(unsigned int);
    /// \return Minimum number of elements among all input data sets
    size_t MinNumElements() const;
    /// Perform any actions necessary to finish filtering.
    int Finalize() const;

    unsigned int Npassed()   const { return Nresult_[PASSED]; }
    unsigned int Nfiltered() const { return Nresult_[FILTERED]; }
  private:
    typedef std::vector<double> Darray;
    typedef std::vector<DataSet*> SetArray;

    inline double GetInpValue(unsigned int, unsigned int) const;

    // The data set integer value for ResultType; PASSED is 1, FILTERED is 0
    static const int ResultValue[];

    DataSet_integer* filterSet_;  ///< Output DataSet containing for each index 1 for OK, 0 for filtered out (not multi).
    DataSet_1D* SetToBeFiltered_; ///< Optional 1D data set to be filtered.
    DataSet_1D* FilteredSet_;     ///< Output 1D data set resulting from filter.
    DataSetList* masterDSL_;      ///< Pointer to the master DataSetList.
    Darray Xvals_;                ///< X values for FilteredSet_.
    Darray Max_;                  ///< Only allow values less than these
    Darray Min_;                  ///< Only allow values greater than these
    SetArray inpSets_;            ///< Sets to base criteria on
    SetArray outSets_;            ///< Sets containing for each index 1 for OK, 0 for filtered out (multi).
    unsigned int Nresult_[2];     ///< Number of indices not filtered/filtered out.
    int debug_;                   ///< Debug level
    int outIdx_;                  ///< Current output index.
    bool multi_;                  ///< If false, filter based on each set. If true, filter each set.
};
#endif
