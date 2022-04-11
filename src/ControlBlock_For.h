#ifndef INC_CONTROLBLOCK_FOR_H
#define INC_CONTROLBLOCK_FOR_H
#include "ControlBlock.h"
class ForLoop;
/// Hold one or more 'for' loops of various types 
class ControlBlock_For : public ControlBlock {
  public:
    ControlBlock_For() {}
    ~ControlBlock_For();
    void Help() const;
    void Help(ArgList&) const;
    DispatchObject* Alloc() const { return (DispatchObject*)new ControlBlock_For(); }

    int SetupBlock(CpptrajState&, ArgList&);
    bool EndBlock(ArgList const&) const;
    void AddCommand(ArgList const& c) { commands_.push_back(c); }
    const_iterator begin() const { return commands_.begin(); }
    const_iterator end()   const { return commands_.end();   }
    int Start(DataSetList const&);
    DoneType CheckDone(DataSetList&);
  private:
    typedef std::vector<ForLoop*> Marray;
    Marray Vars_; ///< Hold all for loops for this block
    ArgArray commands_;
};
#endif
