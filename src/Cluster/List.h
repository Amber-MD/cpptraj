#ifndef INC_CLUSTER_LIST_H
#define INC_CLUSTER_LIST_H
#include <list>
#include "Node.h"
namespace Cpptraj {
namespace Cluster {

/// Hold all individual clusters.
/** Currently implemented as an STL list since sorting and erasing are more
  * efficient.
  */
class List {
  public:
    List() {}
    /// Iterator over clusters
    typedef std::list<Node>::iterator cluster_it;
    cluster_it begin() { return clusters_.begin(); }
    cluster_it end()   { return clusters_.end();   }
    /// Const Iterator over clusters
    typedef std::list<Node>::const_iterator cluster_iterator;
    const cluster_iterator begincluster() const { return clusters_.begin(); }
    const cluster_iterator endcluster()   const { return clusters_.end();   }
  private:
    typedef std::list<Node> Narray;
    Narray clusters_; ///< Hold all clusters.
    Node noise_;      ///< Hold any frames classified as noise.
};

} /* END namespace Cluster */
} /* END namespace Cpptraj */

#endif
