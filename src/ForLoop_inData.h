#ifndef INC_FORLOOP_INDATA_H
#define INC_FORLOOP_INDATA_H
#include "ForLoop.h"
/// Loop over elements of a DataSet
class ForLoop_inData : public ForLoop {
  public:
    ForLoop_inData();

    static void helpText();

    int SetupFor(CpptrajState&, ArgList&);
    int BeginFor(DataSetList const&);
    bool EndFor(DataSetList&);
  private:
    DataSet* set_;     ///< Set with elements to iterate over
    unsigned int idx_; ///< Current index into set_
};
#endif
