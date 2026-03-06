#include "MetalCenterFinder.h"
#include "../ArgList.h"
#include "../CpptrajStdio.h"
#include "../DistRoutines.h"
#include "../Topology.h"
#include <cmath>

using namespace Cpptraj::Structure;

/** CONSTRUCTOR */
MetalCenterFinder::MetalCenterFinder() :
  dcut2_(9.0),
  debug_(0)
{}

/** Init with args */
int MetalCenterFinder::InitMetalCenters(ArgList& argIn, int debugIn)
{
  debug_ = debugIn;
  std::string metalMaskStr = argIn.GetStringKey("metalmask");
  if (metalMaskStr.empty()) {
    metalMaskStr.assign("@ZN,MG");
    mprintf("\tUsing default metal mask string.\n");
  }
  if (metalMask_.SetMaskString( metalMaskStr )) {
    mprinterr("Error: Invalid mask '%s' given for 'metalmask'\n", metalMaskStr.c_str());
    return 1;
  }

  std::string coordAtomMaskStr = argIn.GetStringKey("coordatommask");
  if (coordAtomMaskStr.empty()) {
    coordAtomMaskStr.assign("@/O,S");
    mprintf("\tUsing default coordinating atom mask string.\n");
  }
  if (coordAtomMask_.SetMaskString( coordAtomMaskStr )) {
    mprinterr("Error: Invalid mask '%s' given for 'coordatommask'\n", coordAtomMaskStr.c_str());
    return 1;
  }

  double distcut = argIn.getKeyDouble("mcdist", 3.0);
  if (distcut <= 0) {
    mprinterr("Error: Invalid distance cutoff: %g\n", distcut);
    return 1;
  }
  dcut2_ = distcut*distcut;

  return 0;
}

/** Print info to stdout */
void MetalCenterFinder::PrintMetalCenterInfo() const {
  mprintf("\tMetal center mask: %s\n", metalMask_.MaskString());
  mprintf("\tCoordinating atom mask: %s\n", coordAtomMask_.MaskString());
  mprintf("\tDistance cutoff: %g Ang.\n", sqrt(dcut2_));
}

/** Find metal centers. */
int MetalCenterFinder::FindMetalCenters(Topology const& topIn, Frame const& frameIn)
{
  metalCenters_.clear();

  if (topIn.SetupIntegerMask( metalMask_, frameIn )) {
    mprinterr("Error: Could not set up metal center mask '%s'\n", metalMask_.MaskString());
    return 1;
  }
  if (metalMask_.None()) {
    mprintf("Warning: Nothing selected by metal center mask '%s'\n", metalMask_.MaskString());
    return 0;
  }
  metalMask_.MaskInfo();

  if (topIn.SetupIntegerMask( coordAtomMask_, frameIn )) {
    mprinterr("Error: Could not set up metal center mask '%s'\n", coordAtomMask_.MaskString());
    return 1;
  }
  if (coordAtomMask_.None()) {
    mprintf("Warning: Nothing selected by coordinating atom mask '%s'\n", coordAtomMask_.MaskString());
    return 0;
  }
  coordAtomMask_.MaskInfo();

  // OUTER loop over coordinating atoms
  for (int idx0 = 0; idx0 != coordAtomMask_.Nselected(); idx0++)
  {
    int coordAt = coordAtomMask_[idx0];
    const double* xyz0 = frameIn.XYZ( coordAt );
    // INNER loop over metal center atoms
    for (int idx1 = 0; idx1 != metalMask_.Nselected(); idx1++)
    {
      int metalAt = metalMask_[idx1];
      const double* xyz1 = frameIn.XYZ( metalAt );

      double dist2 = DIST2_NoImage( xyz0, xyz1 );
      if (dist2 < dcut2_) {
        mprintf("\tPotential metal center at %s, coordinating atom %s, dist %f Ang\n",
                topIn.AtomMaskName(metalAt).c_str(),
                topIn.AtomMaskName(coordAt).c_str(),
                sqrt(dist2));
        MCmap::iterator it = metalCenters_.lower_bound( metalAt );
        if (it == metalCenters_.end() || it->first != metalAt) {
          if (debug_ > 0)
            mprintf("DEBUG: New metal center.\n");
          it = metalCenters_.insert(it, MCpair(metalAt, Iarray(1, coordAt)));
        } else {
          if (debug_ > 0)
            mprintf("DEBUG: Existing metal center.\n");
          it->second.push_back( coordAt );
        }
      }
    } // END inner loop over metal center atoms
  } // END outer loop over coordinating atoms

  return 0;
}

/** Print metal centers to stdout */
void MetalCenterFinder::PrintMetalCenters(Topology const& topIn) const {
  mprintf("\t%zu metal centers.\n", metalCenters_.size());
  for (MCmap::const_iterator it = metalCenters_.begin(); it != metalCenters_.end(); ++it)
  {
    mprintf("\t  %s (%zu coordinating atoms):",
            topIn.TruncAtomResNameOnumId(it->first).c_str(),
            it->second.size());
    for (Iarray::const_iterator jt = it->second.begin(); jt != it->second.end(); ++jt)
      mprintf(" %s", topIn.TruncAtomResNameOnumId(*jt).c_str());
    mprintf("\n");
  }
} 
