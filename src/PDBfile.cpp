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
  if (OpenWrite( nameIn )) return 1;
  return 0;
}

void PDBfile::WriteHET(int res, double x, double y, double z) {
  pdb_write_ATOM(IOptr(), PDBtype::PDBHET, anum_++, "XX", "XXX", ' ', res, x, y, z);
}

void PDBfile::WriteATOM(int res, double x, double y, double z, const char* resnameIn,
                        double Occ)
{
  pdb_write_ATOM(IOptr(), PDBtype::PDBATOM, anum_++, "XX", resnameIn, ' ',
                 res, x, y, z, (float)Occ, 0, "", false);
}
