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
class DataFile;
/// Can be used to count/filter out data set elements based on user-defined criteria.
class DataFilter {
  public:
    /// Result of FilterIndex()
    enum ResultType { PASSED = 0, FILTERED };
    /// CONSTRUCTOR
    DataFilter();
    /// Print keywords recognized by InitFilter();
    static void PrintKeywords();
    /// Print a more detail description of keywords.
    static void PrintKeywordDescriptions();
    /// Process arguments, get sets to filter
    int InitFilter(ArgList&, DataSetList&, DataFileList&, int);
    /// \return 1 if specified index was filtered, 0 otherwise
    ResultType FilterIndex(unsigned int);
    /// \return Minimum number of elements among all input data sets
    size_t MinNumElements() const;
    /// Print input sets and min/max values to stdout.
    void PrintInputSets() const;
    /// Perform any actions necessary to finish filtering.
    int Finalize() const;

    /// \return True if only creating output filter sets for multiple input sets.
    bool IsMulti()            const { return multi_; }
    /// \return Number of frames passed (not multi only)
    unsigned int Npassed()    const { return Nresult_[PASSED]; }
    /// \return Number of frames filtered out (not multi only)
    unsigned int Nfiltered()  const { return Nresult_[FILTERED]; }
    /// \return Number of input sets
    unsigned int NinputSets() const { return inpSets_.size(); }
    /// \return Pointer to file that filter sets will be written to (0 if no file)
    DataFile* OutputFile()    const { return maxminfile_; }
  private:
    typedef std::vector<double> Darray;
    typedef std::vector<DataSet*> SetArray;

    inline double GetInpValue(unsigned int, unsigned int) const;

    // The data set integer value for ResultType; PASSED is 1, FILTERED is 0
    static const int ResultValue[];

    DataSet_integer* filterSet_;  ///< Output DataSet containing for each index 1 for OK, 0 for filtered out (not multi).
    DataSet_1D* SetToBeFiltered_; ///< Optional 1D data set to be filtered.
    DataSet_1D* FilteredSet_;     ///< Output 1D data set resulting from filter.
    DataSet* npassedSet_;         ///< Output data set holding number of frames passed (!multi)
    DataSet* nfilteredSet_;       ///< Output data set holding number of frames filtered out (!multi)
    DataFile* maxminfile_;        ///< File to write filter (integer) sets to.
    DataFile* countFile_;         ///< File to write final counts to.
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
