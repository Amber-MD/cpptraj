#include "Pucker_PuckerSearch.h"
#include "CpptrajStdio.h"

using namespace Cpptraj;

Pucker::PuckerSearch::PuckerSearch() {}

/** Indicate we want to search for the specified pre-defined pucker. */
int Pucker::PuckerSearch::SearchFor(Type ptype) {
  PuckerToken::NameArray names;
  if (ptype == NUCLEIC) {
    // Amber nucleic acid
    names.push_back("C1'");
    names.push_back("C2'");
    names.push_back("C3'");
    names.push_back("C4'");
    names.push_back("O4'");
  } else if (ptype == FURANOSE) {
    // Glycam furanose
    names.push_back("C1");
    names.push_back("C2");
    names.push_back("C3");
    names.push_back("C4");
    names.push_back("O4");
  } else if (ptype == PYRANOSE) {
    names.push_back("C1");
    names.push_back("C2");
    names.push_back("C3");
    names.push_back("C4");
    names.push_back("C5");
    names.push_back("O5");
  } else {
    mprinterr("Internal Error: PuckerSearch::SearchFor(): Unhandled pucker type.\n");
    return 1;
  }
  puckersToSearchFor_.push_back( PuckerToken(names) );
  return 0;
}


