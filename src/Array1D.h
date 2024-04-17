#ifndef INC_ARRAY1D_H
#define INC_ARRAY1D_H
#include <vector>
class DataSetList;
class DataSet;
class DataSet_1D;
class ArgList;
/// Hold 1D Scalar DataSets
class Array1D {
  public:
    Array1D() {}
    Array1D(const Array1D&);
    Array1D& operator=(const Array1D&);
    Array1D(DataSetList const&);
    int push_back( DataSet* );
    DataSet_1D* const& operator[](int idx) const { return array_[idx];       }
    DataSet_1D* operator[](int idx)              { return array_[idx];       }
    bool empty()                           const { return array_.empty();    }
    typedef std::vector<DataSet_1D*>::const_iterator const_iterator;
    const_iterator begin()                 const { return array_.begin();    }
    const_iterator end()                   const { return array_.end();      }
    size_t size()                          const { return array_.size();     }
    void clear()                                 { array_.clear();           }
    void SortArray1D();
    int AddDataSets(DataSetList const&);
    int AddTorsionSets(DataSetList const&);
    /// Add sets to the list, clear first
    int AddSetsFromArgs(ArgList const&, DataSetList const&);
    /// Append sets to the list
    int AppendSetsFromArgs(ArgList const&, DataSetList const&);
    /// \return the internal array
    std::vector<DataSet_1D*> const& Array() const { return array_; }
  private:
    int addSetsFromArgs(ArgList const&, DataSetList const&, bool);

    std::vector<DataSet_1D*> array_;
};
#endif
