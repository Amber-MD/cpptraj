#include "Deprecated.h"
#include "CpptrajStdio.h"

void Deprecated_MinDist::Help() const {
  mprinterr("  Use the 'nativecontacts' action instead.\n");
}

void Deprecated_Hbond::Help() const {
  mprinterr("  Hydrogen bond acceptors and donors are defined within the 'hbond' action.\n");
}

void Deprecated_TopSearch::Help() const {
  mprinterr("  Bonds and/or molecules are automatically searched for if needed.\n");
}

void Deprecated_ParmBondInfo::Help() const {
  mprinterr("  Use bonds, bondinfo, or printbonds instead.\n");
}

void Deprecated_ParmResInfo::Help() const {
  mprinterr("  Use resinfo instead.\n");
}

void Deprecated_ParmMolInfo::Help() const {
  mprinterr("  Use molinfo instead.\n");
}

void Deprecated_AvgCoord::Help() const {
  mprinterr("  Use 'vector center' (optionally with keyword 'magnitude') instead.\n");
}

void Deprecated_DihScan::Help() const {
  mprinterr("  Use the 'permutedihedrals' command instead.\n"
            "  See 'help permutedihedrals' for more details.\n");
}
