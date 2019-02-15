#ifndef INC_CLUSTER_RESULTS_H
#define INC_CLUSTER_RESULTS_H
#include "../ArgList.h"
#include "List.h"
namespace Cpptraj {
namespace Cluster {

//FIXME Should this be allocated and kept inside the Metric?

/// Abstract base class for handling results specific to input data type.
class Results {
  public:
    enum Type { COORDS = 0 };

    Results(Type t) : type_(t) {}
    virtual ~Results() {}

    virtual int GetOptions(ArgList&) = 0;
    virtual void Info() const = 0;
    virtual int DoOutput(List const&) const = 0;
  private:
    Type type_;
};

}
}
#endif
