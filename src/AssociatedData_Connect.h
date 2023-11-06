#ifndef INC_ASSOCIATEDDATA_CONNECT_H
#define INC_ASSOCIATEDDATA_CONNECT_H
#include <vector>
#include "AssociatedData.h"
/// For holding unit connect atoms (COORDS data sets) 
class AssociatedData_Connect : public AssociatedData {
    typedef std::vector<int> Iarray;
  public:
    /// Empty CONSTRUCTOR
    AssociatedData_Connect() : AssociatedData(CONNECT) {}
    /// CONSTRUCTOR - head and tail atom
    AssociatedData_Connect(int, int);
    // ----- Inherited functions -------
    static const char* HelpText;
    int ProcessAdataArgs(ArgList&);
    AssociatedData* Copy() const { return new AssociatedData_Connect(*this); }
    void Ainfo() const;
    // ---------------------------------
    void AddConnectAtom(int at) { connect_.push_back( at ); }
    unsigned int NconnectAtoms() const { return connect_.size(); }
    Iarray const& Connect() const { return connect_; }

  private:
    Iarray connect_;
};
#endif
