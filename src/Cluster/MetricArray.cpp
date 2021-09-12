#include "MetricArray.h"
#include "Metric.h"

/** CONSTRUCTOR */
Cpptraj::Cluster::MetricArray::MetricArray() {}

/** COPY CONSTRUCTOR */
Cpptraj::Cluster::MetricArray::MetricArray(MetricArray const& rhs) {
  metrics_.reserve( rhs.metrics_.size() );
  for (std::vector<Metric*>::const_iterator it = rhs.metrics_.begin(); it != rhs.metrics_.end(); ++it)
    metrics_.push_back( (*it)->Copy() );
}

/** ASSIGNMENT */
Cpptraj::Cluster::MetricArray&
  Cpptraj::Cluster::MetricArray::operator=(MetricArray const& rhs)
{
  if (this == &rhs) return *this;
  metrics_.clear();
  metrics_.reserve( rhs.metrics_.size() );
  for (std::vector<Metric*>::const_iterator it = rhs.metrics_.begin(); it != rhs.metrics_.end(); ++it)
    metrics_.push_back( (*it)->Copy() );
  return *this;
}

