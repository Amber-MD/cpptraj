#ifndef INC_FORLOOP_LIST_H
#define INC_FORLOOP_LIST_H
#include <vector>
#include <string>
#include "ForLoop.h"
class ForLoop_list : public ForLoop {
  public:
    ForLoop_list() {}

    static void helpText();

    int SetupFor(CpptrajState&, ArgList&);
    int BeginFor(DataSetList const&);
    bool EndFor(DataSetList&);
  private:
    typedef std::vector<std::string> Sarray;
    Sarray Names_;               ///< List of strings, potentially need expansion/replacement
    Sarray List_;                ///< (LIST only) List of strings to iterate over.
    Sarray::const_iterator sdx_; ///< (LIST only) Iterator to current list item.
};
#endif
