#include "FunctionalGroup.h"
#include "../AtomMap.h"
#include "../CpptrajStdio.h"
#include "../Frame.h"
#include "../Topology.h"

using namespace Cpptraj::Structure;

FunctionalGroup::FunctionalGroup() {}

/** CONSTRUCTOR - resname */
/*FunctionalGroup::FunctionalGroup(NameType const& rname) :
  resname_(rname)
{}*/

/** Add atom name and ID */
void FunctionalGroup::AddAtom(NameType const& aname, std::string const& id) {
  anames_.push_back( aname );
  atomIDs_.push_back( id );
}

/** Corresponds to FunctionalGroupType */
const char* FunctionalGroup::FunctionalGroupStr_[] = {
  "SO3", "CH3", "Acetyl", "OH", "OCH3", "Unrecognized"
};

/** Clear all information. */
void FunctionalGroup::Clear() {
  resname_ = NameType("");
  anames_.clear();
  atomIDs_.clear();
  linkAtomElt_ = Atom::UNKNOWN_ELEMENT;
  chargeAtom_ = Atom::UNKNOWN_ELEMENT;
  chargeOffset_ = 0;
}

/** Set up functional group from given topology. */
int FunctionalGroup::SetupFromTop(Topology const& groupTop, Atom::AtomicElementType linkAtomEltIn) {
  Clear();
  if (groupTop.Nres() != 1) {
    mprinterr("Internal Error: FunctionalGroup::SetupFromTop: Expected 1 res, got %i\n", groupTop.Nres());
    return 1;
  }
  resname_ = groupTop.Res(0).Name();
  linkAtomElt_ = linkAtomEltIn;
  AtomMap groupmap;
  //groupmap.SetDebug(10); // DEBUG
  groupmap.Setup( groupTop, Frame() );
  groupmap.DetermineAtomIDs();

  for (int idx = 0; idx != groupmap.Natom(); idx++)
    AddAtom( groupmap[idx].Name(), groupmap[idx].Unique() );

  return 0;
}

/** Print info to stdout*/
void FunctionalGroup::PrintInfo() const {
  mprintf("FxnGroup '%s' :", *resname_);
  for (unsigned int idx = 0; idx != anames_.size(); idx++)
    mprintf(" '%s'[%s]", *(anames_[idx]), atomIDs_[idx].c_str());
  mprintf(" Linked via '%s'\n", Atom::ElementName(linkAtomElt_));
}
