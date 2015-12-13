#ifndef INC_CMD_H
#define INC_CMD_H
#include <vector>
#include <string>
#include "DispatchObject.h"
/** Class for holding a command and associated keywords.
  * NOTE: This class does NOT contain a destructor because otherwise object_
  *       might be freed when copying or assigning. Currently CmdList is
  *       responsible for freeing object_ memory.
  */
class Cmd {
  public:
    typedef std::vector< std::string > Sarray; // TODO put in common header?
    /// Command destinations. EXEcute, ACTion, ANAlysis, DEPrecated.
    enum DestType { EXE = 0, ACT, ANA, DEP };
    /// CONSTRUCTOR
    Cmd() : object_(0), dest_(EXE) {}
    /// CONSTRUCTOR - takes destination, DispatchObject pointer, and keywords.
    Cmd(DispatchObject* o, Sarray k, DestType d) : object_(o), keywords_(k), dest_(d) {}
    /// \return command destination 
    DestType Destination() const { return dest_; }
    /// Iterator for command keywords.
    typedef Sarray::const_iterator key_iterator;
    /// \return iterator to beginning of command keywords.
    key_iterator keysBegin() const { return keywords_.begin(); }
    /// \return iterator to end of command keywords.
    key_iterator keysEnd()   const { return keywords_.end();   }
    /// \return true if given key matches any of this commands keywords.
    bool KeyMatches(const char*) const;
    /// Free DispatchObject
    void Clear();
    /// \return true if no DispatchObject
    bool Empty() const { return object_ == 0; }
    /// \return const reference to underlying DispatchObject
    DispatchObject const& Obj() const { return *object_; }
    /// Execute Help for underlying DispatchObject
    void Help() const { object_->Help(); }
    /// \return Copy of underlying DispatchObject
    DispatchObject* Alloc() const { return object_->Alloc(); }
  private:
    DispatchObject* object_; ///< Pointer to DispatchObject
    Sarray keywords_;        ///< Keywords for this command
    DestType dest_;          ///< The command destination.
};
#endif
