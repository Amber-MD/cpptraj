#include "CentroidArray.h"
#include "Centroid.h"

/** CONSTRUCTOR */
Cpptraj::Cluster::CentroidArray::CentroidArray() {}

/** DESTRUCTOR */
Cpptraj::Cluster::CentroidArray::~CentroidArray() {
  Clear();
}

/** COPY CONSTRUCTOR */
Cpptraj::Cluster::CentroidArray::CentroidArray(CentroidArray const& rhs) {
  centroids_.reserve( rhs.centroids_.size() );
  for (std::vector<Centroid*>::const_iterator it = rhs.centroids_.begin(); it != rhs.centroids_.end(); ++it)
    centroids_.push_back( (*it)->Copy() );
}

/** ASSIGNMENT */
Cpptraj::Cluster::CentroidArray&
  Cpptraj::Cluster::CentroidArray::operator=(CentroidArray const& rhs)
{
  if (this == &rhs) return *this;
  centroids_.clear();
  centroids_.reserve( rhs.centroids_.size() );
  for (std::vector<Centroid*>::const_iterator it = rhs.centroids_.begin(); it != rhs.centroids_.end(); ++it)
    centroids_.push_back( (*it)->Copy() );
  return *this;
}

/** Clear the centroid array. */
void Cpptraj::Cluster::CentroidArray::Clear() {
  for (std::vector<Centroid*>::iterator it = centroids_.begin(); it != centroids_.end(); ++it)
    delete *it;
}
