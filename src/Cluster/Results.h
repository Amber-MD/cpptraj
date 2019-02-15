#ifndef INC_CLUSTER_RESULTS_H
#define INC_CLUSTER_RESULTS_H
#include "../ArgList.h"
namespace Cpptraj {
namespace Cluster {

/// Abstract base class for handling results specific to input data type.
class Results {
  public:
    enum Type { COORDS = 0 };

    Results(Type t) : type_(t) {}
    virtual ~Results() {}

    virtual int GetOptions(ArgList&) = 0;
    virtual void Info() const = 0;
    virtual int DoOutput() const = 0;
  private:
    Type type_;
};

}
}
#endif
