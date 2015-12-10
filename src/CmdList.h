#ifndef INC_CMDLIST_H
#define INC_CMDLIST_H
#include "Cmd.h"
class CmdList {
    typedef std::vector<Cmd> Carray;
  public:
    CmdList() {}
    ~CmdList();
    void Clear();
    void Add(Cmd const& cIn) { cList_.push_back( cIn ); }
    typedef Carray::const_iterator const_iterator;
    const_iterator begin() const { return cList_.begin(); }
    const_iterator end()   const { return cList_.end();   }
    Cmd const& Back()      const { return cList_.back();  }
  private:
    Carray cList_;
};
#endif
