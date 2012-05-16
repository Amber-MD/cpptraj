#include "PDBfile.h"

PDBfile::PDBfile() :
  anum_(1)
{}

PDBfile::~PDBfile() {
  CloseFile();
}

int PDBfile::OpenPDB(std::string const& nameIn) {
  if (OpenWrite( nameIn )) return 1;
  return 0;
}

int PDBfile::OpenPDB(const char* nameIn) {
  if (SetupWrite( nameIn, 0)) return 1;
  if (OpenFile()) return 1;
  return 0;
}

void PDBfile::WriteHET(int res, double x, double y, double z) {
  pdb_write_ATOM(IO, PDBtype::PDBHET, anum_, "XX", "XXX", ' ', res, x, y, z);
  ++anum_;
}

