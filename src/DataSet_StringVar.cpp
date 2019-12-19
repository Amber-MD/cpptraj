#include "DataSet_StringVar.h"
#include "CpptrajStdio.h"

/// CONSTRUCTOR
DataSet_StringVar::DataSet_StringVar() :
  DataSet(STRINGVAR, GENERIC, TextFormat(), 0)
{ }

size_t DataSet_StringVar::MemUsageInBytes() const {
  return sizeof(std::string) + (var_.size() * sizeof(char));
}
