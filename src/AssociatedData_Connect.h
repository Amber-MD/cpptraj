#ifndef INC_ASSOCIATEDDATA_CONNECT_H
#define INC_ASSOCIATEDDATA_CONNECT_H
#include <vector>
#include "AssociatedData.h"
/// For holding unit connect atoms (COORDS data sets) 
class AssociatedData_Connect : public AssociatedData {
  public:
    AssociatedData_Connect() : AssociatedData(CONNECT) {}
    static const char* HelpText;

    void AddConnectAtom(int at) { connect_.push_back( at ); }

    AssociatedData* Copy() const { return new AssociatedData_Connect(*this); }
    void Ainfo() const;
  private:
    typedef std::vector<int> Iarray;
    Iarray connect_;
};
#endif
