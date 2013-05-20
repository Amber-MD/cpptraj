#ifndef INC_ARRAY1D_H
#define INC_ARRAY1D_H
#include "DataSet_1D.h"
/// Hold 1D DataSets
class Array1D {
  public:
    Array1D() : errorMsg_("") {}
    Array1D(DataSetList const&);
    size_t DetermineMax();
    int push_back( DataSet_1D* const& );
    DataSet_1D* const& operator[](int idx) const { return array_[idx];       }
    bool empty()                           const { return array_.empty();    }
    const char* Error()                    const { return errorMsg_.c_str(); }
    typedef std::vector<DataSet_1D*>::const_iterator const_iterator;
    const_iterator begin()                 const { return array_.begin();    }
    const_iterator end()                   const { return array_.end();      }
    size_t size()                          const { return array_.size();     }
    void clear()                                 { array_.clear();           }
    int AddDataSets(DataSetList const&);
  private:
    std::vector<DataSet_1D*> array_;
    std::string errorMsg_;
};
// CONSTRUCTOR
Array1D::Array1D(DataSetList const& SetList) : errorMsg_("") {
  AddDataSets( SetList );
  if (array_.empty())
    errorMsg_.assign("Internal Error: DataIO called with no 1D data sets.");
}
// push_back()
int Array1D::push_back( DataSet_1D* const& val ) {
  if (val->Ndim() == 1)
    array_.push_back( val );
  else
    return 1;
  return 0;
}
// Array1D::AddDataSets()
int Array1D::AddDataSets(DataSetList const& SetList) {
  for (DataSetList::const_iterator ds = SetList.begin(); ds != SetList.end(); ++ds)
    if ( push_back( (DataSet_1D*)*ds ) ) {
      errorMsg_.assign("Internal Error: DataIO not set up for mixed dim writes.");
      array_.clear();
      return 1;
    }
  return 0;
}
// DetermineMax()
size_t Array1D::DetermineMax() {
  size_t maxFrames = 0L;
  for (std::vector<DataSet_1D*>::const_iterator set = array_.begin(); set != array_.end(); ++set)
    if ( (*set)->Size() > maxFrames )
      maxFrames = (*set)->Size();
  return maxFrames;
}
#endif
