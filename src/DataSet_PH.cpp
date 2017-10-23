#include "DataSet_PH.h"

DataSet_PH::DataSet_PH() :
  // 0 dim indicates DataSet-specific write
  DataSet(PH, GENERIC, TextFormat(TextFormat::DOUBLE, 10, 4), 0),
  nframes_(0)
{}
