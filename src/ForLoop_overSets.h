#ifndef INC_FORLOOP_OVERSETS_H
#define INC_FORLOOP_OVERSETS_H
#include <vector>
#include <string>
#include "ForLoop.h"
class ForLoop_overSets : public ForLoop {
  public:
    ForLoop_overSets() {}

    static void helpText();

    int SetupFor(CpptrajState&, ArgList&);
    int BeginFor(DataSetList const&);
    bool EndFor(DataSetList&);
  private:
    typedef std::vector<std::string> Sarray;
    Sarray Names_;               ///< List of non-expanded names
    Sarray List_;                ///< List of data set names to iterate over.
    Sarray::const_iterator sdx_; ///< Iterator to current list item.
};
#endif
