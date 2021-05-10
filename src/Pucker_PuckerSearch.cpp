#include "Pucker_PuckerSearch.h"
#include "CpptrajStdio.h"
#include "ArgList.h"
#include "Range.h"
#include "Topology.h"

using namespace Cpptraj;

Pucker::PuckerSearch::PuckerSearch() {}

/** Recognized pucker keywords. */
const char* Pucker::PuckerSearch::Keywords_[] = {
  "nucleic",  // NUCLEIC
  "furanose", // FURANOSE
  "pyranose"  // PYRANOSE
};

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
  puckersToSearchFor_.push_back( PuckerToken(Keywords_[ptype], names) );
  return 0;
}

/** See if ArgList has any recognized pucker type keywords. */
int Pucker::PuckerSearch::SearchForArgs(ArgList& argIn) {
  for (int i = 0; i != (int)NTYPES; i++) {
    if (argIn.hasKey( Keywords_[i] ))
      SearchFor( (Type)i );
  }
  return 0;
}

/** Define a custom pucker argument from ArgList:
  *   'puckertype <name>:<a0>:...:<aN>'
  */
int Pucker::PuckerSearch::SearchForNewTypeArgs(ArgList& argIn) {
  std::string puckertype_arg = argIn.GetStringKey("puckertype");
  while (!puckertype_arg.empty()) {
    ArgList puckertype(puckertype_arg, ":");
    if (puckertype.Nargs() < 6) {
      mprinterr("Error: Malformed puckertype arg '%s': expected at least 6 args, got %i\n",
                puckertype_arg.c_str(), puckertype.Nargs());
      return 1;
    }
    PuckerToken::NameArray atomNames;
    atomNames.reserve(puckertype.Nargs()-1);
    for (int iarg = 1; iarg != puckertype.Nargs(); iarg++)
      atomNames.push_back( puckertype[iarg] );
    SearchForNewType(puckertype[0], atomNames);
    puckertype_arg = argIn.GetStringKey("puckertype");
  }
  return 0;
}

/** Define a custom pucker */
int Pucker::PuckerSearch::SearchForNewType(std::string const& name, PuckerToken::NameArray const& atomNames)
{
  for (std::vector<PuckerToken>::const_iterator tkn = puckersToSearchFor_.begin();
                                                tkn != puckersToSearchFor_.end(); ++tkn)
    if (tkn->Name() == name) {
      mprintf("Warning: Pucker type %s already defined.\n", name.c_str());
      return 1;
    }
  puckersToSearchFor_.push_back( PuckerToken(name, atomNames) );
  return 0;
}

/** If no puckers selected yet, select all. */
int Pucker::PuckerSearch::SearchForAll() {
  if (!puckersToSearchFor_.empty()) return 0;
  for (int ptype = 0; ptype != (int)NTYPES; ptype++)
    SearchFor( (Type)ptype );
  return 0;
}

/** Print all set up types to stdout. */
void Pucker::PuckerSearch::PrintTypes() const {
  for (std::vector<PuckerToken>::const_iterator it = puckersToSearchFor_.begin();
                                                it != puckersToSearchFor_.end(); ++it)
    mprintf(" %s", it->Name().c_str());
}

/** Find all defined puckers in a residue range */
int Pucker::PuckerSearch::FindPuckers(Topology const& currentParm, Range const& rangeIn)
{
  foundPuckers_.clear();
  for (Range::const_iterator res = rangeIn.begin(); res != rangeIn.end(); ++res)
  {
    for (std::vector<PuckerToken>::const_iterator tkn = puckersToSearchFor_.begin();
                                                  tkn != puckersToSearchFor_.end(); ++tkn)
    {
      PuckerMask puckerMask = tkn->FindPuckerAtoms(currentParm, *res);
      if (!puckerMask.None()) {
        foundPuckers_.push_back( puckerMask );
      }
    }
  }
  if (foundPuckers_.empty()) {
    mprintf("Warning: No puckers selected for topology %s\n", currentParm.c_str());
    return 1;
  }
  return 0;
}
