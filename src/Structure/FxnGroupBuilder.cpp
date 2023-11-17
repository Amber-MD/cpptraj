#include "FxnGroupBuilder.h"
#include "StructureRoutines.h"
#include "Sugar.h"
#include "../CpptrajStdio.h"
#include "../Frame.h"
#include "../Topology.h"

using namespace Cpptraj::Structure;

FxnGroupBuilder::FxnGroupBuilder(int debugIn) :
  debug_(debugIn)
{
  AddGroups();
}

/** Add functional groups to the array. */
int FxnGroupBuilder::AddGroups() {
  FunctionalGroup fg;
  // Sulfate
  //       O1
  //       |
  // (O) - S - O2
  //       |
  //       O3
  Topology so3top;
  so3top.AddTopAtom( Atom("S",  "S "), Residue("SO3", 1, ' ', ' ') );
  so3top.AddTopAtom( Atom("O1", "O "), Residue("SO3", 1, ' ', ' ') );
  so3top.AddTopAtom( Atom("O2", "O "), Residue("SO3", 1, ' ', ' ') );
  so3top.AddTopAtom( Atom("O3", "O "), Residue("SO3", 1, ' ', ' ') );
  so3top.AddBond(0, 1);
  so3top.AddBond(0, 2);
  so3top.AddBond(0, 3);
  if (fg.SetupFromTop( so3top, Atom::OXYGEN ))
    return 1;
  functionalGroups_.push_back( fg );
  // Hydroxyl
  // (C1) - O1
  Topology ohtop;
  ohtop.AddTopAtom( Atom("O1", "O"), Residue("ROH", 1, ' ', ' ') );
  if (fg.SetupFromTop( ohtop, Atom::CARBON ))
    return 1;
  functionalGroups_.push_back( fg );
  // Hydroxyl
  // (C1) - O1 - HO1
  ohtop.AddTopAtom( Atom("HO1", "H"), Residue("ROH", 1, ' ', ' ') );
  ohtop.AddBond(0, 1);
  if (fg.SetupFromTop( ohtop, Atom::CARBON ))
    return 1;
  functionalGroups_.push_back( fg );

  // Print all groups
  for (std::vector<FunctionalGroup>::const_iterator it = functionalGroups_.begin();
                                                    it != functionalGroups_.end(); ++it)
    it->PrintInfo();
  return 0;
} 

/// recursive function to visit all bonded atoms in a group
void visitGroupAtoms(std::vector<int>& groupAtoms, int currentAtom, Topology const& topIn, std::vector<bool>& visited)
{
  visited[currentAtom] = true;
  groupAtoms.push_back( currentAtom );
  for (Atom::bond_iterator bat = topIn[currentAtom].bondbegin();
                           bat != topIn[currentAtom].bondend(); ++bat)
  {
    if (!visited[*bat])
      visitGroupAtoms(groupAtoms, *bat, topIn, visited);
  }
}

/** Determine if atoms connected atoms to atIdx (which in turn is connected
  * to linkAtIdx), not including linkAtIdx or ignoreAtoms, form a
  * recognized functional group.
  */
FxnGroupBuilder::FGarray::const_iterator
  FxnGroupBuilder::GetGroup(Iarray& groupAtoms, Iarray const& ignoreAtoms, int atIdx, int linkAtIdx, Topology const& topIn)
const
{
  std::vector<FunctionalGroup>::const_iterator FG_it = functionalGroups_.end();
  Atom const& startAtom = topIn[atIdx];
  Atom const& linkAtom = topIn[linkAtIdx];
  // atIdx and linkAtIdx must be in the same residue
  if (startAtom.ResNum() != linkAtom.ResNum()) return FG_it;
  groupAtoms.clear();
  // Mark all atoms inside the residue (except the link atom) as not visited.
  std::vector<bool> visited( topIn.Natom(), true );
  for (int at = topIn.Res(startAtom.ResNum()).FirstAtom();
           at != topIn.Res(startAtom.ResNum()).LastAtom(); ++at)
  {
    if (at != linkAtIdx)
      visited[at] = false;
  }
  // Mark atoms to ignore as visited. If atIdx is to be ignored, exit.
  for (Iarray::const_iterator it = ignoreAtoms.begin(); it != ignoreAtoms.end(); ++it) {
    if (*it == atIdx) return FG_it;
    visited[*it] = true;
  }

  visitGroupAtoms( groupAtoms, atIdx, topIn, visited );
  mprintf("\tPotential functional group starting at [%s]-[%s]:", topIn.AtomMaskName(atIdx).c_str(), topIn.AtomMaskName(linkAtIdx).c_str());
  for (Iarray::const_iterator it = groupAtoms.begin(); it != groupAtoms.end(); ++it)
    mprintf(" %s", topIn.AtomMaskName( *it ).c_str());
  mprintf("\n");
  if (groupAtoms.empty()) return FG_it;

  // Create topology containing only groupAtoms
  Iarray oldToNew( topIn.Natom(), -1 );
  Topology groupTop;
  int newIdx = 0;
  for (Iarray::const_iterator it = groupAtoms.begin(); it != groupAtoms.end(); ++it, newIdx++)
  {
    groupTop.AddTopAtom( Atom(topIn[*it].Name(), topIn[*it].ElementName()), topIn.Res(topIn[*it].ResNum()) );
    oldToNew[*it] = newIdx;
  }
  // Add bonds to groupTop
  for (Iarray::const_iterator it = groupAtoms.begin(); it != groupAtoms.end(); ++it) {
    if (oldToNew[*it] > -1) {
      for (Atom::bond_iterator jt = topIn[*it].bondbegin(); jt != topIn[*it].bondend(); ++jt) {
        if (oldToNew[*jt] > -1) {
          groupTop.AddBond( oldToNew[*it], oldToNew[*jt], -1 );
        }
      }
    }
  }
  if (debug_ > 0)
    groupTop.Summary();
  FunctionalGroup fg;
  if (fg.SetupFromTop( groupTop, linkAtom.Element() )) {
    mprinterr("Internal Error: Set up functional group topology failed.\n");
    return FG_it;
  }
  if (debug_ > 0) fg.PrintInfo();

  // Try to find a match
  for (FG_it = functionalGroups_.begin(); FG_it != functionalGroups_.end(); ++FG_it)
  {
    if ( FG_it->Match( fg ) ) {
       break;
    }
  }
  if (FG_it != functionalGroups_.end()) {
    mprintf("\tFunctional group match found. ");
    FG_it->PrintInfo();
  }
  return FG_it;
}

/** \return Type of group represented by the atom atIdx. */
FunctionalGroup::Type
  FxnGroupBuilder::IdFunctionalGroup_Silent(Iarray& selected, int rnum,
                                            int atIdx, int linkAtIdx, Topology const& topIn)
const
{
  selected.clear();
  // ----- Sulfur ----------------------
  if (topIn[atIdx].Element() == Atom::SULFUR &&
      topIn[atIdx].ResNum() == rnum &&
      topIn[atIdx].Nbonds() == 4) {
    //so3_idx = atIdx;
    selected.push_back(atIdx);
    // All 4 bonds must be to oxygen
    int bonds_to_o = 0;
    for (Atom::bond_iterator bat = topIn[atIdx].bondbegin();
                             bat != topIn[atIdx].bondend(); ++bat)
    {
      if (topIn[*bat].Element() == Atom::OXYGEN &&
          topIn[*bat].ResNum() == rnum) {
        bonds_to_o++;
        if (*bat != linkAtIdx)
          selected.push_back( *bat );
      }
    }
    if (bonds_to_o == 4) {
      return FunctionalGroup::G_SO3;
    }
  // ----- Carbon ----------------------
  } else if (topIn[atIdx].Element() == Atom::CARBON &&
             topIn[atIdx].ResNum() == rnum) {
    //so3_idx = atIdx;
    selected.push_back(atIdx);
    if (topIn[atIdx].Nbonds() == 4) {
      // If 4 bonds, 3 must be to hydrogen for CH3
      int bonds_to_h = 0;
      for (Atom::bond_iterator bat = topIn[atIdx].bondbegin();
                               bat != topIn[atIdx].bondend(); ++bat)
      {
        if (topIn[*bat].Element() == Atom::HYDROGEN) {
          bonds_to_h++;
          selected.push_back(*bat);
        }
      }
      if (bonds_to_h == 3) {
        return FunctionalGroup::G_CH3;
      }
    } else if (topIn[atIdx].Nbonds() == 1)
      return FunctionalGroup::G_CH3;
    } else if (topIn[atIdx].Nbonds() == 3) {
      //             O
      //             ||
      // Check for - C - CH3
      bool has_o = false;
      bool has_ch3 = false;
      Iarray ch3_atoms;
      for (Atom::bond_iterator bat = topIn[atIdx].bondbegin();
                               bat != topIn[atIdx].bondend(); ++bat)
      {
        if (*bat != linkAtIdx) {
          if (topIn[*bat].Element() == Atom::OXYGEN) {
            has_o = true;
            selected.push_back(*bat);
          } else if (topIn[*bat].Element() == Atom::CARBON) {
            FunctionalGroup::Type fgt = IdFunctionalGroup_Silent(ch3_atoms, rnum, *bat, atIdx, topIn);
            has_ch3 = (fgt == FunctionalGroup::G_CH3);
          } else {
            // If bonded to anything else, not Acetyl
            has_o = false;
            has_ch3 = false;
            break;
          }
        }
      }
      if (has_o && has_ch3) {
        for (Iarray::const_iterator cit = ch3_atoms.begin(); cit != ch3_atoms.end(); ++cit)
          selected.push_back(*cit);
        return FunctionalGroup::G_ACX;
      }
  // ----- Oxygen ----------------------
  } else if (topIn[atIdx].Element() == Atom::OXYGEN &&
             topIn[atIdx].ResNum() == rnum) {
    selected.push_back( atIdx );
    if (topIn[atIdx].Nbonds() == 1)
      // If only 1 bond, -OH
      return FunctionalGroup::G_OH;
    else if (topIn[atIdx].Nbonds() == 2) {
      int bonded_atom;
      if (topIn[atIdx].Bond(0) == linkAtIdx)
        bonded_atom = topIn[atIdx].Bond(1);
      else
        bonded_atom = topIn[atIdx].Bond(0);
      if (topIn[bonded_atom].Element() == Atom::HYDROGEN) {
        // -OH
        selected.push_back(bonded_atom);
        return FunctionalGroup::G_OH;
      } else if (topIn[bonded_atom].Element() == Atom::CARBON) {
        // Might be -OCH3
        int c_idx = bonded_atom;
        selected.push_back(c_idx);
        // Check for only 1 bond or 3 bonds to hydrogen
        int bonds_to_h = 0;
        for (Atom::bond_iterator bat = topIn[c_idx].bondbegin();
                                 bat != topIn[c_idx].bondend(); ++bat)
        {
          if (topIn[*bat].Element() == Atom::HYDROGEN) {
            bonds_to_h++;
            selected.push_back(*bat);
          }
        }
        if (topIn[c_idx].Nbonds() == 1 || bonds_to_h == 3)
          return FunctionalGroup::G_OME;
      }
    }

  }
  selected.clear();      
  return FunctionalGroup::UNRECOGNIZED_GROUP;
}

/** Identify functional group and print to stdout. */
FunctionalGroup::Type
  FxnGroupBuilder::IdFunctionalGroup(Iarray& selected, int rnum, int atIdx, int linkAtIdx,
                                     Topology const& topIn)
const
{
  FunctionalGroup::Type groupType = IdFunctionalGroup_Silent(selected, rnum, atIdx, linkAtIdx, topIn);
  if (groupType != FunctionalGroup::UNRECOGNIZED_GROUP) {
    mprintf("\tFound %s group centered on atom '%s %s' bonded to '%s'\n",
            FunctionalGroup::typeString(groupType),
            topIn.TruncResNameOnumId(rnum).c_str(),
            *(topIn[atIdx].Name()), *(topIn[linkAtIdx].Name()));
  }

  return groupType;
}

/** Check for functional groups that need to be separate residues. */
int FxnGroupBuilder::CheckForFunctionalGroups(Sugar& sugar, Topology& topIn, Frame& frameIn)
const
{
  int rnum = sugar.ResNum(topIn);
  int original_at0 = topIn.Res(rnum).FirstAtom();
  int original_at1 = topIn.Res(rnum).LastAtom();
  std::string sugarName = topIn.TruncResNameOnumId(rnum);
  if (debug_ > 0) mprintf("DEBUG: Functional group check: %s\n", sugarName.c_str());

  bool atomsRemain = true;
  while (atomsRemain) {
    // Even if the residue is split, rnum will always refer to the original
    // sugar since the split SO3 will come AFTER this residue.
    // Find an oxygen that is both bound to a chain carbon and an SO3 group
    // in this residue.
    int o_idx = -1;
    int so3_idx = -1;
    Iarray selected;
    FunctionalGroup::Type groupType = FunctionalGroup::UNRECOGNIZED_GROUP;
    Iarray::const_iterator cat = sugar.ChainAtoms().begin();
    for (; cat != sugar.ChainAtoms().end(); ++cat)
    {
      //mprintf("\t%s\n", *(topIn[*cat].Name())); // DEBUG
      for (Atom::bond_iterator oat = topIn[*cat].bondbegin();
                               oat != topIn[*cat].bondend(); ++oat)
      {

        o_idx = -1;
        if (topIn[*oat].Element() == Atom::OXYGEN &&
            topIn[*oat].ResNum() == rnum &&
            topIn[*oat].Nbonds() > 1) {
          o_idx = *oat;
          // Is this oxygen bound to a recognized group?
          for (Atom::bond_iterator sat = topIn[*oat].bondbegin();
                                   sat != topIn[*oat].bondend(); ++sat)
          {
            // DEBUG
            Iarray groupAtoms;
            GetGroup(groupAtoms, sugar.ChainAtoms(), *sat, *oat, topIn);
            // DEBUG
            groupType = IdFunctionalGroup(selected, rnum, *sat, o_idx, topIn);
            if (groupType != FunctionalGroup::UNRECOGNIZED_GROUP) {
              so3_idx = *sat;
              break;
            }
          } // END loop over bonds to oxygen
          if (so3_idx != -1) break;
        } // END atom is oxygen
      } // END loop over bonds to carbon
      if (so3_idx != -1) break;
    } // END loop over chain atoms
    if (cat == sugar.ChainAtoms().end())
      atomsRemain = false;
    else {
      // sanity check
      if (so3_idx == -1 || o_idx == -1 || selected.empty()) {
        mprinterr("Internal Error: Functional group index is negative.\n");
        return 1;
      }
      
      // Change the atom names
      std::string newResName;
      if (groupType == FunctionalGroup::G_SO3) {
        ChangeAtomName(topIn.SetAtom(selected[0]), "S1");
        ChangeAtomName(topIn.SetAtom(selected[1]), "O1");
        ChangeAtomName(topIn.SetAtom(selected[2]), "O2");
        ChangeAtomName(topIn.SetAtom(selected[3]), "O3");
        newResName = "SO3";
      } else if (groupType == FunctionalGroup::G_CH3) {
        ChangeAtomName(topIn.SetAtom(selected[0]), "CH3");
        if (selected.size() > 1) {
          ChangeAtomName(topIn.SetAtom(selected[1]), "H1");
          ChangeAtomName(topIn.SetAtom(selected[2]), "H2");
          ChangeAtomName(topIn.SetAtom(selected[3]), "H3");
        }
        newResName = "MEX";
      } else if (groupType == FunctionalGroup::G_ACX) {
        ChangeAtomName(topIn.SetAtom(selected[0]), "C1A");
        ChangeAtomName(topIn.SetAtom(selected[1]), "O1A");
        ChangeAtomName(topIn.SetAtom(selected[2]), "C2A");
        if (selected.size() > 3) {
          ChangeAtomName(topIn.SetAtom(selected[3]), "H1A");
          ChangeAtomName(topIn.SetAtom(selected[4]), "H2A");
          ChangeAtomName(topIn.SetAtom(selected[5]), "H3A");
        }
      } else {
        mprinterr("Internal Error: Unhandled group in CheckForFunctionalGroups()\n");
        return 1;
      }
      // Create array with group selected
      AtomMask FxnGrpMask(selected, topIn.Natom());
      // Split the sulfate into a new residue named newResName for Glycam.
      // This may involve reordering atoms within the residue, but not
      // any other atoms, so we should not have to update other sugars.
      //mprintf("DEBUG: Before split: %s\n", topIn.AtomMaskName(sugar.RingOxygenAtom()).c_str());
      Iarray atomMap;
      if (topIn.SplitResidue(FxnGrpMask, newResName, atomMap)) {
        mprinterr("Error: Could not split functional group from residue '%s'.\n",
                  sugarName.c_str());
        return 1;
      }
      // Set the split residue as terminal
      topIn.SetRes(rnum+1).SetTerminal(true);
      // Reorder the frame to match
      Frame oldFrame = frameIn;
      frameIn.SetCoordinatesByMap( oldFrame, atomMap );
      // Remap the sugar indices
      //mprintf("DEBUG: Before remap: %s\n", topIn.AtomMaskName(sugar.RingOxygenAtom()).c_str());
      sugar.RemapIndices( atomMap, original_at0, original_at1 );
      //mprintf("DEBUG: After remap: %s\n", topIn.AtomMaskName(sugar.RingOxygenAtom()).c_str());
    }
    //atomsRemain = false; // DEBUG
  } // END while atoms remain

  return 0;
}

/** See if the sugar anomeric carbon is actually terminal and needs
  * to be a separate ROH residue.
  */
int FxnGroupBuilder::CheckIfSugarIsTerminal(Sugar& sugar, Topology& topIn, Frame& frameIn)
const
{
  int rnum = sugar.ResNum(topIn);
  int anomericAtom = sugar.AnomericAtom();
  int ringOxygen = sugar.RingOxygenAtom();
  std::string sugarName = topIn.TruncResNameOnumId(rnum);

  // Is the anomeric carbon bonded to an oxygen that is part of this residue.
  int o1_atom = -1;
  for (Atom::bond_iterator bat = topIn[anomericAtom].bondbegin();
                           bat != topIn[anomericAtom].bondend(); ++bat)
  {
    if (topIn[*bat].ResNum() == rnum &&
        *bat != ringOxygen &&
        topIn[*bat].Element() == Atom::OXYGEN) {
      if (o1_atom != -1) {
        mprintf("Warning: Anomeric atom '%s %s' bonded to more than 1 oxygen.\n",
                sugarName.c_str(), *(topIn[anomericAtom].Name()));
        o1_atom = -1;
        break;
      } else
        o1_atom = *bat;
    }
  }
  if (o1_atom == -1) return 0;
  //mprintf("DEBUG: Terminal check: %s O1 atom: '%s'\n", sugarName.c_str(), topIn.AtomMaskName(o1_atom).c_str());
  // DEBUG
  Iarray groupAtoms;
  GetGroup(groupAtoms, sugar.ChainAtoms(), o1_atom, anomericAtom, topIn);
  // DEBUG

  Iarray selected;
  FunctionalGroup::Type groupType = IdFunctionalGroup(selected, rnum, o1_atom, anomericAtom, topIn);

  std::string newResName;
  if (groupType == FunctionalGroup::G_OH) {
    newResName = "ROH";
    // Change atom names
    ChangeAtomName(topIn.SetAtom(selected[0]), "O1");
    if (selected.size() > 1)
      ChangeAtomName(topIn.SetAtom(selected[1]), "HO1");
  } else if (groupType == FunctionalGroup::G_OME) {
    newResName = "OME";
    // Change atom names
    ChangeAtomName(topIn.SetAtom(selected[0]), "O");
    ChangeAtomName(topIn.SetAtom(selected[1]), "CH3");
    if (selected.size() > 2) {
      ChangeAtomName(topIn.SetAtom(selected[2]), "H1");
      ChangeAtomName(topIn.SetAtom(selected[3]), "H2");
      ChangeAtomName(topIn.SetAtom(selected[4]), "H3");
    }
  } else if (groupType != FunctionalGroup::UNRECOGNIZED_GROUP) {
    mprinterr("Internal Error: Unhandled group in CheckIfSugarIsTerminal\n");
    return 1;
  }
  //mprintf("\t  Will split into %s group.\n", newResName.c_str());

  // If terminal group is unrecognized, this could just be a regular O1 linkage
  if (selected.empty()) {
    if (debug_ > 0)
      mprintf("DEBUG: Sugar '%s' has unrecognized terminal group.\n", sugarName.c_str());
    return 0;
  }
  
  // Create atom mask
  AtomMask ROH(selected, topIn.Natom());

  // Split the hydroxyl into a new residue named ROH for Glycam.
  // This may involve reordering atoms within the residue, but not
  // any other atoms, so we should not have to update SugarIndices.
  int original_at0 = topIn.Res(rnum).FirstAtom();
  int original_at1 = topIn.Res(rnum).LastAtom();
  Iarray atomMap;
  if (topIn.SplitResidue(ROH, newResName, atomMap)) {
    mprinterr("Error: Could not split the residue '%s'.\n", sugarName.c_str());
    return 1;
  }
  // Set the split residue as terminal
  topIn.SetRes(rnum+1).SetTerminal(true);
  // DEBUG
  //for (int at = original_at0; at != original_at1; at++)
  //  mprintf("DEBUG:\t\tAtomMap[%i] = %i\n", at, atomMap[at]);
  // Reorder the frame to match
  Frame oldFrame = frameIn;
  frameIn.SetCoordinatesByMap( oldFrame, atomMap );
  // Remap the sugar indices
  sugar.RemapIndices( atomMap, original_at0, original_at1 );
  return 0;
}
