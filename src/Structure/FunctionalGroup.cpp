#include "FunctionalGroup.h"

using namespace Cpptraj::Structure;

FunctionalGroup::FunctionalGroup() {}

/** CONSTRUCTOR - resname */
FunctionalGroup::FunctionalGroup(NameType const& rname) :
  resname_(rname)
{}

/** Add atom name and ID */
void FunctionalGroup::AddAtom(NameType const& aname, std::string const& id) {
  anames_.push_back( aname );
  atomIDs_.push_back( id );
}

/** Corresponds to FunctionalGroupType */
const char* FunctionalGroup::FunctionalGroupStr_[] = {
  "SO3", "CH3", "Acetyl", "OH", "OCH3", "Unrecognized"
};

