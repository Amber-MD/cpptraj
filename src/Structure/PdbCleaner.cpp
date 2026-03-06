#include "PdbCleaner.h"
#include "../ArgList.h"
#include "../CpptrajStdio.h"
#include "../Topology.h"
#include <map>

using namespace Cpptraj::Structure;

/** CONSTRUCTOR */
PdbCleaner::PdbCleaner() :
  debug_(0),
  remove_water_(false),
  remove_h_(false)
{}

/** Print help text. */
void PdbCleaner::Help() {
  mprintf("\t[nowat [watermask <watermask>] [noh]\n"
          "\t[keepaltloc {<alt loc ID>|highestocc}]\n"
          "\t[stripmask <stripmask>]\n");
}

/** Initialize */
int PdbCleaner::InitPdbCleaner(ArgList& argIn, std::string const& solventResNameIn,
                               Iarray const& pdbResToRemoveIn)
{
  resnumsToRemove_ = pdbResToRemoveIn;
  remove_water_ = argIn.hasKey("nowat");
  if (solventResNameIn.empty())
    waterMask_ = ":HOH";
  else
    waterMask_ = argIn.GetStringKey("watermask", ":" + solventResNameIn);
  remove_h_ = argIn.hasKey("noh");
  altLocArg_ = argIn.GetStringKey("keepaltloc");
  if (!altLocArg_.empty()) {
    if (altLocArg_ != "highestocc" &&
        altLocArg_.size() > 1)
    {
      mprinterr("Error: Invalid keyword for 'keepaltloc' '%s'; must be 'highestocc' or 1 character.\n",
                altLocArg_.c_str());
      return 1;
    }
  }
  stripMask_ = argIn.GetStringKey("stripmask");

  return 0;
}

/** Setup */
int PdbCleaner::SetupPdbCleaner(Topology const& topIn) {
  // If keeping highest alt loc, check that alt locs and occupancies are present.
  if (altLocArg_ == "highestocc") {
    if (topIn.AtomAltLoc().empty()) {
      mprintf("Warning: 'highestocc' specified but no atom alternate location info.\n");
      altLocArg_.clear();
    } else if (topIn.Occupancy().empty()) {
      mprintf("Warning: 'highestocc' specified but no atom occupancy info.\n");
      altLocArg_.clear();
    }
  }
  // Check if alternate atom location IDs are present 
  if (!topIn.AtomAltLoc().empty()) {
    // For LEaP, must have only 1 atom alternate location
    char firstAltLoc = ' ';
    for (std::vector<char>::const_iterator altLocId = topIn.AtomAltLoc().begin();
                                           altLocId != topIn.AtomAltLoc().end();
                                         ++altLocId)
    {
      if (firstAltLoc == ' ') {
        // Find first non-blank alternate location ID
        if (*altLocId != ' ')
          firstAltLoc = *altLocId;
      } else if (*altLocId != ' ' && *altLocId != firstAltLoc) {
        // Choose a default if necessary
        if (altLocArg_.empty()) {
          altLocArg_.assign(1, firstAltLoc);
          mprintf("Warning: '%s' has atoms with multiple alternate location IDs.\n"
                  "Warning: Keeping only '%s'.\n"
                  "Warning: To choose a specific location to keep use the 'keepaltloc <char>'\n"
                  "Warning:  keyword.\n", topIn.c_str(), altLocArg_.c_str());
         }
         break;
      }
    }
  }
  return 0;
}

/** Print info to stdout. */
void PdbCleaner::PdbCleanerInfo() const {
  if (remove_water_)
    mprintf("\tRemoving solvent. Solvent mask= '%s'\n", waterMask_.c_str());
  if (remove_h_)
    mprintf("\tRemoving hydrogens.\n");
  if (!altLocArg_.empty())
    mprintf("\tIf present, keeping only alternate atom locations denoted by '%s'\n", altLocArg_.c_str());
  if (!stripMask_.empty())
    mprintf("\tRemoving atoms in mask '%s'\n", stripMask_.c_str());
}

/** Modify coords according to user wishes. */
int PdbCleaner::ModifyCoords( Topology& topIn, Frame& frameIn )
const
{
  // Create a mask denoting which atoms will be kept.
  std::vector<bool> atomsToKeep( topIn.Natom(), true );
  // Previously-determined array of residues to remove
  for (Iarray::const_iterator rnum = resnumsToRemove_.begin();
                              rnum != resnumsToRemove_.end(); ++rnum)
  {
    Residue const& res = topIn.Res( *rnum );
    mprintf("\tRemoving %s\n", topIn.TruncResNameOnumId( *rnum ).c_str());
    for (int at = res.FirstAtom(); at != res.LastAtom(); at++)
      atomsToKeep[at] = false;
  }
  // User-specified strip mask
  if (!stripMask_.empty()) {
    AtomMask mask;
    if (mask.SetMaskString( stripMask_ )) {
      mprinterr("Error: Invalid mask string '%s'\n", stripMask_.c_str());
      return 1;
    }
    if (topIn.SetupIntegerMask( mask )) return 1;
    mask.MaskInfo();
    if (!mask.None()) {
      for (AtomMask::const_iterator atm = mask.begin(); atm != mask.end(); ++atm)
        atomsToKeep[*atm] = false;
    }
  }
  if (remove_water_) {
    // Do not use cpptraj definition of solvent in case we have e.g. bound HOH
    AtomMask mask;
    if (mask.SetMaskString( waterMask_ )) {
      mprinterr("Error: Invalid solvent mask string '%s'\n", waterMask_.c_str());
      return 1;
    }
    if (topIn.SetupIntegerMask( mask )) return 1;
    mask.MaskInfo();
    if (!mask.None()) {
      for (AtomMask::const_iterator atm = mask.begin(); atm != mask.end(); ++atm)
        atomsToKeep[*atm] = false;
    }

  }
  // Identify alternate atom location groups.
  if (!altLocArg_.empty()) {
    if (topIn.AtomAltLoc().empty()) {
      mprintf("\tNo alternate atom locations.\n");
    } else {
      // Map atom name to atom indices
      typedef std::map<NameType, std::vector<int>> AlocMapType;
      AlocMapType alocMap;
      for (int rnum = 0; rnum != topIn.Nres(); rnum++) {
        alocMap.clear();
        for (int at = topIn.Res(rnum).FirstAtom(); at != topIn.Res(rnum).LastAtom(); at++) {
          if (topIn.AtomAltLoc()[at] != ' ') {
            AlocMapType::iterator it = alocMap.find( topIn[at].Name() );
            if (it == alocMap.end()) {
              alocMap.insert( std::pair<NameType, std::vector<int>>( topIn[at].Name(),
                                                                     std::vector<int>(1, at) ));
            } else {
              it->second.push_back( at );
            }
          }
        } // END loop over atoms in residue
        if (!alocMap.empty()) {
          if (debug_ > 0)
            mprintf("DEBUG: Alternate loc. for %s\n", topIn.TruncResNameOnumId(rnum).c_str());
          // Loop over atoms with alternate locations
          for (AlocMapType::const_iterator it = alocMap.begin(); it != alocMap.end(); ++it) {
            if (debug_ > 0) {
              // Print all alternate atoms
              mprintf("\t'%s'", *(it->first));
              for (std::vector<int>::const_iterator at = it->second.begin();
                                                    at != it->second.end(); ++at)
                mprintf(" %s[%c]", *(topIn[*at].Name()), topIn.AtomAltLoc()[*at]);
              mprintf("\n");
            }
            // For each, choose which location to keep.
            if (altLocArg_.size() == 1) {
              // Keep only specified character
              char altLocChar = altLocArg_[0];
              for (std::vector<int>::const_iterator at = it->second.begin();
                                                    at != it->second.end(); ++at)
                if (topIn.AtomAltLoc()[*at] != altLocChar)
                  atomsToKeep[*at] = false;
            } else {
              // Keep highest occupancy
              if (topIn.Occupancy().empty()) {
                mprintf("\tNo occupancy.\n"); // TODO error?
              } else {
                int highestOccAt = -1;
                float highestOcc = 0;
                for (std::vector<int>::const_iterator at = it->second.begin();
                                                      at != it->second.end(); ++at)
                {
                  if (highestOccAt == -1) {
                    highestOccAt = *at;
                    highestOcc = topIn.Occupancy()[*at];
                  } else if (topIn.Occupancy()[*at] > highestOcc) {
                    highestOccAt = *at;
                    highestOcc = topIn.Occupancy()[*at];
                  }
                }
                // Set everything beside highest occ to false
                for (std::vector<int>::const_iterator at = it->second.begin();
                                                      at != it->second.end(); ++at)
                  if (*at != highestOccAt)
                    atomsToKeep[*at] = false;
              }
            }
          } // END loop over atoms with alternate locations
        }
      } // END loop over residue numbers
    }
  }
  // Denote hydrogens if removing
  if (remove_h_) {
    for (int at = 0; at != topIn.Natom(); at++)
      if (topIn[at].Element() == Atom::HYDROGEN)
        atomsToKeep[at] = false;
  }

  // Set up mask of only kept atoms.
  AtomMask keptAtoms;
  keptAtoms.SetNatoms( topIn.Natom() );
  for (int idx = 0; idx != topIn.Natom(); idx++) {
    if (atomsToKeep[idx])
      keptAtoms.AddSelectedAtom(idx);
  }
  if (keptAtoms.Nselected() == topIn.Natom())
    // Keeping everything, no modifications
    return 0;
  // Modify top/frame
  Topology* newTop = topIn.modifyStateByMask( keptAtoms );
  if (newTop == 0) {
    mprinterr("Error: Could not create new topology.\n");
    return 1;
  }
  newTop->Brief("After removing atoms:");
  Frame newFrame;
  newFrame.SetupFrameV(newTop->Atoms(), frameIn.CoordsInfo());
  newFrame.SetFrame(frameIn, keptAtoms);

  topIn = *newTop;
  frameIn = newFrame;
  delete newTop;

  return 0;
}
