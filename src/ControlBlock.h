#ifndef INC_CONTROLBLOCK_H
#define INC_CONTROLBLOCK_H
#include <vector>
#include <string>
#include "DispatchObject.h"
class CpptrajState;
class DataSetList;
class ArgList;
/// Abstract base class for control block structures.
class ControlBlock : public DispatchObject {
  public:
    typedef std::vector<ArgList> ArgArray;
    typedef ArgArray::const_iterator const_iterator;
    /// Control block states
    enum DoneType { DONE = 0, NOT_DONE, ERROR };

    ControlBlock() : DispatchObject(CONTROL) {}
    virtual ~ControlBlock() {}
    /// \return Description of control block.
    std::string const& Description() const { return description_; }
    /// Set up control block.
    virtual int SetupBlock(CpptrajState&, ArgList&) = 0;
    /// Check for control block end command.
    virtual bool EndBlock(ArgList const&) const = 0;
    /// Add command to control block.
    virtual void AddCommand(ArgList const&) = 0;
    /// \return iterator to first command in the block.
    virtual const_iterator begin() const = 0;
    /// \return iterator to last command in the block.
    virtual const_iterator end() const = 0;
    /// Start control block. Init internal variables if necessary.
    virtual int Start(DataSetList const&) = 0;
    /// Add/update variables and increment, check block state, create sets if needed
    virtual DoneType CheckDone(DataSetList&) = 0;
  protected:
    std::string description_; ///< Describe control TODO private?
};
#endif
