#include "DataSet_integer_disk.h"

DataSet_integer_disk::DataSet_integer_disk() :
  ncid_(-1),
  nvals_(0)
{
  start_[0] = 0;
  count_[0] = 0;
}

