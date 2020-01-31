#include "DataSet_StringVar.h"
#include "CpptrajStdio.h"

/// CONSTRUCTOR
DataSet_StringVar::DataSet_StringVar() :
  DataSet(STRINGVAR, GENERIC, TextFormat(), 0)
{ }

size_t DataSet_StringVar::MemUsageInBytes() const {
  return sizeof(std::string) + (var_.size() * sizeof(char));
}

void DataSet_StringVar::Add(size_t n, const void* ptr) {
  const char* cptr = (const char*)ptr;
  var_ = std::string(cptr);
}
