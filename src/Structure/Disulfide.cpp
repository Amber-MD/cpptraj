#include "Disulfide.h"
#include "ResStatArray.h"
#include "StructureRoutines.h"
#include "../ArgList.h"
#include "../CpptrajStdio.h"
#include "../DistRoutines.h"
#include "../Frame.h"
#include "../Topology.h"
#include <algorithm> //sort

using namespace Cpptraj::Structure;

/** CONSTRUCTOR */
Disulfide::Disulfide() :
  disulfidecut_(2.5),
  searchForNewDisulfides_(true),
  addNewBonds_(NO_ADD_BONDS)
{}

/** Keywords recognized by InitDisulfide */
const char* Disulfide::keywords_ =
  "\t[{nodisulfides |\n"
  "\t  existingdisulfides |\n"
  "\t  [cysmask <cysmask>] [disulfidecut <cut>] [newcysname <name>]}]\n";

/** Init with args */
int Disulfide::InitDisulfide(ArgList& argIn, AddType addTypeIn, int debugIn)
{
  addNewBonds_ = addTypeIn;
  disulfidecut_ = argIn.getKeyDouble("disulfidecut", 2.5);
  newcysnamestr_ = argIn.GetStringKey("newcysname", "CYX");
  cysmaskstr_ = argIn.GetStringKey("cysmask", ":CYS@SG");
  searchForNewDisulfides_ = !argIn.hasKey("existingdisulfides");
  return 0;
}

/** Print info to stdout. */
//void Disulfide::DisulfideInfo() const {


/** Search for disulfide bonds. */
int Disulfide::SearchForDisulfides(ResStatArray& resStat, Topology& topIn,
                                   Frame const& frameIn,
                                   std::vector<BondType>& LeapBonds)
const
{
  return searchForDisulfides(resStat, disulfidecut_, newcysnamestr_,
                             cysmaskstr_, searchForNewDisulfides_,
                             topIn, frameIn, LeapBonds, addNewBonds_);
}

/** Search for disulfide bonds. */
int Disulfide::searchForDisulfides(ResStatArray& resStat,
                                             double disulfidecut,
                                             std::string const& newcysnamestr,
                                             std::string const& cysmaskstr,
                                             bool searchForNewDisulfides,
                                             Topology& topIn, Frame const& frameIn,
                                             std::vector<BondType>& LeapBonds,
                                             AddType addNewBonds)
{
  // Disulfide search
  typedef std::vector<int> Iarray;
  NameType newcysname(newcysnamestr);
  mprintf("\tCysteine residues involved in disulfide bonds will be changed to: %s\n", *newcysname);
  if (searchForNewDisulfides)
    mprintf("\tSearching for disulfide bonds with a cutoff of %g Ang.\n", disulfidecut);
  else
    mprintf("\tOnly using existing disulfide bonds, will not search for new ones.\n");

  AtomMask cysmask;
  if (cysmask.SetMaskString( cysmaskstr )) {
    mprinterr("Error: Could not set up CYS mask string %s\n", cysmaskstr.c_str());
    return 1;
  }
  if (topIn.SetupIntegerMask( cysmask )) return 1; 
  cysmask.MaskInfo();
  typedef std::vector<bool> Barray;
  Barray cysSulfurHasH( cysmask.Nselected(), false );
  if (cysmask.None())
    mprintf("Warning: No cysteine sulfur atoms selected by %s\n", cysmaskstr.c_str());
  else {
    // Sanity check - warn if non-sulfurs detected
    for (AtomMask::const_iterator at = cysmask.begin(); at != cysmask.end(); ++at)
    {
      if (topIn[*at].Element() != Atom::SULFUR)
        mprintf("Warning: Atom '%s' does not appear to be sulfur.\n",
                topIn.ResNameNumAtomNameNum(*at).c_str());
      // Check if the sulfur is bonded to a hydrogen.
      for (Atom::bond_iterator bat = topIn[*at].bondbegin();
                               bat != topIn[*at].bondend(); ++bat)
      {
        if ( topIn[*bat].Element() == Atom::HYDROGEN )
          cysSulfurHasH[at - cysmask.begin()] = true;
      }
    }

    int nExistingDisulfides = 0;
    int nDisulfides = 0;
    double cut2 = disulfidecut * disulfidecut;
    // Try to find potential disulfide sites.
    // Keep track of which atoms will be part of disulfide bonds.
    Iarray disulfidePartner( cysmask.Nselected(), -1 );
    // First, check for existing disulfides.
    for (AtomMask::const_iterator at1 = cysmask.begin(); at1 != cysmask.end(); ++at1)
    {
      for (AtomMask::const_iterator at2 = at1 + 1; at2 != cysmask.end(); ++at2)
      {
        if (topIn[*at1].IsBondedTo(*at2)) {
          if (StructureDebugLevel() > 0)
            mprintf("\tExisting disulfide: %s to %s\n",
                    topIn.ResNameNumAtomNameNum(*at1).c_str(),
                    topIn.ResNameNumAtomNameNum(*at2).c_str());
          nExistingDisulfides++;
          int idx1 = (int)(at1 - cysmask.begin());
          int idx2 = (int)(at2 - cysmask.begin());
          disulfidePartner[idx1] = idx2;
          disulfidePartner[idx2] = idx1;
        }
      }
    }
    mprintf("\t%i existing disulfide bonds.\n", nExistingDisulfides);
    // DEBUG - Print current array
    if (StructureDebugLevel() > 1) {
      mprintf("DEBUG: Disulfide partner array after existing:\n");
      for (Iarray::const_iterator it = disulfidePartner.begin(); it != disulfidePartner.end(); ++it)
      {
        mprintf("  S %i [%li]", cysmask[it-disulfidePartner.begin()]+1, it-disulfidePartner.begin());
        if (*it == -1)
          mprintf(" None.\n");
        else
          mprintf(" to S %i [%i]\n", cysmask[*it]+1, *it);
      }
    }
    // Second, search for new disulfides from remaining sulfurs.
    if (searchForNewDisulfides) {
      // Only search with atoms that do not have an existing partner.
      Iarray s_idxs; // Indices into cysmask/disulfidePartner
      for (int idx = 0; idx != cysmask.Nselected(); idx++) {
        if (disulfidePartner[idx] == -1) {
          s_idxs.push_back( idx );
        }
      }
      mprintf("\t%zu sulfur atoms do not have a partner.\n", s_idxs.size());
      if (!s_idxs.empty()) {
        // In some structures, there may be 2 potential disulfide partners
        // within the cutoff. In that case, only choose the shortest.
        // To try to do this as directly as possible, calculate all possible S-S
        // distances, save the ones below the cutoff, sort them from shortest to
        // longest, then assign each that to a disulfide if not already assigned.
        // In this way, shorter S-S distances will be prioritized.
        typedef std::pair<int,int> IdxPair;
        typedef std::pair<double, IdxPair> D2Pair;
        typedef std::vector<D2Pair> D2Array;
        D2Array D2;

        for (Iarray::const_iterator it1 = s_idxs.begin(); it1 != s_idxs.end(); ++it1)
        {
          int at1 = cysmask[*it1];
          for (Iarray::const_iterator it2 = it1 + 1; it2 != s_idxs.end(); ++it2)
          {
            int at2 = cysmask[*it2];
            double r2 = DIST2_NoImage(frameIn.XYZ(at1), frameIn.XYZ(at2));
            if (r2 < cut2)
              D2.push_back( D2Pair(r2, IdxPair(*it1, *it2)) );
          }
        }
        std::sort(D2.begin(), D2.end());
        if (StructureDebugLevel() > 1) {
          mprintf("DEBUG: Sorted S-S array:\n");
          for (D2Array::const_iterator it = D2.begin(); it != D2.end(); ++it)
          {
            int at1 = cysmask[it->second.first];
            int at2 = cysmask[it->second.second];
            mprintf("  %8i - %8i = %g Ang.\n", at1+1, at2+2, sqrt(it->first));
          }
        }
        // All distances in D2 are below the cutoff
        for (D2Array::const_iterator it = D2.begin(); it != D2.end(); ++it)
        {
          if (disulfidePartner[it->second.first] == -1 &&
              disulfidePartner[it->second.second] == -1)
          {
            // Neither index has a partner yet
            int at1 = cysmask[it->second.first];
            int at2 = cysmask[it->second.second];
            mprintf("\t  Potential disulfide: %s to %s (%g Ang.)\n",
                    topIn.ResNameNumAtomNameNum(at1).c_str(),
                    topIn.ResNameNumAtomNameNum(at2).c_str(), sqrt(it->first));
            if (cysSulfurHasH[it->second.first])
              mprintf("Warning: Sulfur %s could disulfide bond but already has a bond to hydrogen. Check for improperly added hydrogen.\n",
                      topIn.ResNameNumAtomNameNum(at1).c_str());
            if (cysSulfurHasH[it->second.second])
              mprintf("Warning: Sulfur %s could disulfide bond but already has a bond to hydrogen. Check for improperly added hydrogen.\n",
                      topIn.ResNameNumAtomNameNum(at2).c_str());
            disulfidePartner[it->second.first ] = it->second.second;
            disulfidePartner[it->second.second] = it->second.first;
            if (addNewBonds == ADD_BONDS)
              topIn.AddBond( at1, at2 );
          }
        } // END loop over sorted distances
      } // END s_idxs not empty()
    } // END search for new disulfides
    // For each sulfur that has a disulfide partner, generate a bond command
    // and change residue name.
    for (Iarray::const_iterator idx1 = disulfidePartner.begin(); idx1 != disulfidePartner.end(); ++idx1)
    {
      if (*idx1 != -1) {
        int at1 = cysmask[idx1-disulfidePartner.begin()];
        int at2 = cysmask[*idx1];
        if (at1 < at2) {
          nDisulfides++;
          // TODO remove bond from Top?
          LeapBonds.push_back( BondType(at1, at2, -1) );
          //outfile->Printf("%s\n", LeapInterface::LeapBond(at1, at2, leapunitname, topIn).c_str());
        }
        ChangeResName(topIn.SetRes(topIn[at1].ResNum()), newcysname);
        resStat[topIn[at1].ResNum()] = ResStatArray::VALIDATED;
      }
    }
    mprintf("\tDetected %i disulfide bonds.\n", nDisulfides);
  }
  return 0;
}

