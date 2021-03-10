#ifndef INC_DATAFILTER_H
#define INC_DATAFILTER_H
#include <vector>
#include "DataSet_double.h"
class ArgList;
class DataSetList;
class DataFileList;
class DataSet;
class DataSet_integer;
/// Can be used to count/filter out data set elements based on user-defined criteria.
class DataFilter {
  public:
    DataFilter();
    /// \return Keywords recognized by InitFilter();
    static const char* Keywords();
    /// Process arguments, get sets to filter
    int InitFilter(ArgList&, DataSetList&, DataFileList&, int);
    /// \return 1 if specified index was filtered, 0 otherwise
    int FilterIndex(int);
  private:
    typedef std::vector<double> Darray;
    typedef std::vector<DataSet*> SetArray;

    DataSet_integer* filterSet_;  ///< Output DataSet containing for each index 1 for OK, 0 for filtered out (not multi).
    DataSet_1D* SetToBeFiltered_; ///< Optional 1D data set to be filtered.
    DataSet_1D* FilteredSet_;     ///< Output 1D data set resulting from filter.
    DataSet_double Xvals_;        ///< X values for FilteredSet_.
    Darray Max_;                  ///< Only allow values less than these
    Darray Min_;                  ///< Only allow values greater than these
    SetArray inpSets_;            ///< Sets to base criteria on
    SetArray outSets_;            ///< Sets containing for each index 1 for OK, 0 for filtered out (multi).
    unsigned int Npassed_;        ///< Number of indices not filtered out.
    unsigned int Nfiltered_;      ///< Number of indices filtered out.
    int debug_;                   ///< Debug level
    bool multi_;                  ///< If false, filter based on each set. If true, filter each set.
};
#endif
