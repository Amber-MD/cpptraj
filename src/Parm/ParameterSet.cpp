#include "ParameterSet.h"
#include "../CpptrajFile.h"
#include "../CpptrajStdio.h"
#include "../Constants.h"
#include "../UpdateParameters.h"

using namespace Cpptraj::Parm;

size_t ParameterSet::DataSize() const {
  unsigned int cmapsize = 0;
  for (CmapGridArray::const_iterator it = CMAP_.begin(); it != CMAP_.end(); ++it)
    cmapsize += it->DataSize();
  return (atomTypes_.DataSize() +
          nbParm_.DataSize() +
          nb14Parm_.DataSize() +
          ljcParm_.DataSize() +
          bondParm_.DataSize() +
          angleParm_.DataSize() +
          ubParm_.DataSize() +
          impParm_.DataSize() +
          dihParm_.DataSize() +
          HBparm_.DataSize() +
          cmapsize +
          (hydrophilicAtomTypes_.size() * NameType::DataSize()) +
          sizeof(bool)
         );
}

/** Write parameters out to file with given name. */
void ParameterSet::Debug(const char* fnameIn) const {
  CpptrajFile Out;
  Out.OpenWrite( fnameIn );
  Print( Out );
  Out.CloseFile();
}

/** Write summary to stdout. */
void ParameterSet::Summary() const {
  mprintf("\tParameter set: %s\n", ParamSetName().c_str());
  mprintf("\t  %zu atom types:", atomTypes_.size());
  for (ParmHolder<AtomType>::const_iterator at = atomTypes_.begin(); at != atomTypes_.end(); ++at)
    mprintf(" %s", *(at->first[0]));
  mprintf("\n");
  // TODO 1-4 types?
  mprintf("\t  %zu LJ 6-12 parameters.\n", nbParm_.size());
  mprintf("\t  %zu LJ 6-12 1-4 parameters.\n", nb14Parm_.size());
  mprintf("\t  %zu LJ C parameters.\n", ljcParm_.size());
  mprintf("\t  %zu LJ 10-12 parameters.\n", HBparm_.size());
  mprintf("\t  %zu bond parameters.\n", bondParm_.size());
  mprintf("\t  %zu angle parameters.\n", angleParm_.size());
  mprintf("\t  %zu UB parameters.\n", ubParm_.size());
  mprintf("\t  %zu dihedral parameters.\n", dihParm_.size());
  mprintf("\t  %zu improper parameters.\n", impParm_.size());
  mprintf("\t  %zu CMAP parameters.\n", CMAP_.size());
}

/** Write parameters out to given file. */
void ParameterSet::Print(CpptrajFile& Out) const {
  if (!name_.empty())
    Out.Printf("Parameter set: %s\n", ParamSetName().c_str());
  if (!NBname_.empty())
    Out.Printf("Nonbond parameters name: %s\n", NBname_.c_str());
  static const char* hybStr[] = {"SP ", "SP2", "SP3", "UNK"};
  Out.Printf("Atom Types:\n");
  // Check if we have hybridizations or elements.
  bool has_hybridizations = false;
  bool has_elements = false;
  for (ParmHolder<AtomType>::const_iterator at = atomTypes_.begin(); at != atomTypes_.end(); ++at) {
    if (at->second.Hybridization() != AtomType::UNKNOWN_HYBRIDIZATION)
      has_hybridizations = true;
    if (at->second.EltStr()[0] != ' ' || at->second.EltStr()[1] != ' ')
      has_elements = true;
    if (has_hybridizations && has_elements)
      break;
  }
  // Atom Type header
  Out.Printf("\t%6s %8s %2s %12s %12s %12s %12s", "Name", "TypeIdx", "LJ", "Radius", "Depth", "Mass", "Polar.");
  if (has_hybridizations)
    Out.Printf(" %3s", "Hyb");
  if (has_elements)
    Out.Printf(" %2s", "El");
  Out.Printf("\n");
  unsigned int atidx = 0;
  for (ParmHolder<AtomType>::const_iterator at = atomTypes_.begin(); at != atomTypes_.end(); ++at, atidx++) {
    Out.Printf("\t%6s %8u %2i %12.4f %12.4f %12.4f %12.4f", *(at->first[0]), atidx, (int)at->second.HasLJ(), at->second.LJ().Radius(), at->second.LJ().Depth(), at->second.Mass(), at->second.Polarizability());
    if (has_hybridizations)
      Out.Printf(" %3s", hybStr[at->second.Hybridization()]);
    if (has_elements)
      Out.Printf(" %2s", at->second.EltStr());
    Out.Printf("\n");
  }
  if (!nbParm_.empty()) {
    Out.Printf("LJ 6-12 parameters:\n");
    Out.Printf("\t%6s %6s : %12s %12s\n", "Type1", "Type2", "A", "B");
    for (ParmHolder<NonbondType>::const_iterator nb = nbParm_.begin(); nb != nbParm_.end(); ++nb)
      Out.Printf("\t%6s %6s : %12.4E %12.4E\n", (*nb->first[0]), (*nb->first[1]), nb->second.A(), nb->second.B());
  }
  if (!nb14Parm_.empty()) {
    Out.Printf("LJ 1-4 parameters:\n");
    Out.Printf("\t%6s %6s : %12s %12s\n", "Type1", "Type2", "A", "B");
    for (ParmHolder<NonbondType>::const_iterator nb = nb14Parm_.begin(); nb != nb14Parm_.end(); ++nb)
      Out.Printf("\t%6s %6s : %12.4E %12.4E\n", (*nb->first[0]), (*nb->first[1]), nb->second.A(), nb->second.B());
  }
  if (!ljcParm_.empty()) {
    Out.Printf("LJ C parameters:\n");
    Out.Printf("\t%6s %6s : %12s\n", "Type1", "Type2", "C");
    for (ParmHolder<NonbondType>::const_iterator nb = nb14Parm_.begin(); nb != nb14Parm_.end(); ++nb)
      Out.Printf("\t%6s %6s : %12.4E\n", (*nb->first[0]), (*nb->first[1]), nb->second);
  }
  if (!HBparm_.empty()) {
    Out.Printf("HB LJ 10-12 parameters:\n");
    Out.Printf("\t%6s %6s : %12s %12s %12s\n", "Type1", "Type2", "Asol", "Bsol", "HBcut");
    for (ParmHolder<HB_ParmType>::const_iterator hb = HBparm_.begin(); hb != HBparm_.end(); ++hb)
      Out.Printf("\t%6s %6s : %12.4f %12.4f %12.4f\n", *(hb->first[0]), *(hb->first[1]), hb->second.Asol(), hb->second.Bsol(), hb->second.HBcut());
  }
  if (!bondParm_.empty()) {
    Out.Printf("Bond parameters:\n");
    Out.Printf("\t%6s %6s : %12s %12s\n", "Type1", "Type2", "Rk", "Req");
    for (ParmHolder<BondParmType>::const_iterator bp = bondParm_.begin(); bp != bondParm_.end(); ++bp)
      Out.Printf("\t%6s %6s : %12.4f %12.4f\n", *(bp->first[0]), *(bp->first[1]), bp->second.Rk(), bp->second.Req());
  }
  if (!angleParm_.empty()) {
    Out.Printf("Angle parameters:\n");
    Out.Printf("\t%6s %6s %6s : %12s %12s\n", "Type1", "Type2", "Type3", "Tk", "Teq");
    for (ParmHolder<AngleParmType>::const_iterator bp = angleParm_.begin(); bp != angleParm_.end(); ++bp)
      Out.Printf("\t%6s %6s %6s : %12.4f %12.4f\n", *(bp->first[0]), *(bp->first[1]), *(bp->first[2]), bp->second.Tk(), bp->second.Teq()*Constants::RADDEG);
  }
  if (!ubParm_.empty()) {
    Out.Printf("UB parameters:\n");
    Out.Printf("\t%6s %6s : %12s %12s\n", "Type1", "Type2", "Uk", "Ueq");
    for (ParmHolder<BondParmType>::const_iterator bp = ubParm_.begin(); bp != ubParm_.end(); ++bp)
      Out.Printf("\t%s %s : %f %f\n", *(bp->first[0]), *(bp->first[1]), bp->second.Rk(), bp->second.Req());
  }
  if (!dihParm_.empty()) {
    Out.Printf("Dihedral parameters:\n");
    Out.Printf("\t%6s %6s %6s %6s %12s %4s %8s %6s %6s\n", "Type1", "Type2", "Type3", "Type4", "Pk", "Pn", "Phase", "SCEE", "SCNB");
    for (DihedralParmHolder::const_iterator it0 = dihParm_.begin(); it0 != dihParm_.end(); ++it0)
      for (DihedralParmArray::const_iterator it1 = it0->second.begin();
                                             it1 != it0->second.end(); ++it1)
        Out.Printf("\t%6s %6s %6s %6s : %12.4f %4.1f %8.2f %6.2f %6.2f\n", *(it0->first[0]), *(it0->first[1]), *(it0->first[2]), *(it0->first[3]), it1->Pk(), it1->Pn(), it1->Phase()*Constants::RADDEG, it1->SCEE(), it1->SCNB());
  }
  if (!impParm_.empty()) {
    Out.Printf("Improper parameters:\n");
    Out.Printf("\t%6s %6s %6s %6s %12s %12s %12s\n", "Type1", "Type2", "Type3", "Type4", "Pk", "Pn", "Phase");
    for (ImproperParmHolder::const_iterator it0 = impParm_.begin(); it0 != impParm_.end(); ++it0)
      for (DihedralParmArray::const_iterator it1 = it0->second.begin();
                                             it1 != it0->second.end(); ++it1)
        Out.Printf("\t%6s %6s %6s %6s : %12.4f %12.4f %12.4f\n", *(it0->first[0]), *(it0->first[1]), *(it0->first[2]), *(it0->first[3]), it1->Pk(), it1->Pn(), it1->Phase()*Constants::RADDEG);
  }
  if (!hydrophilicAtomTypes_.empty()) {
    Out.Printf("Hydrophilic atom types:");
    for (NsetType::const_iterator it = hydrophilicAtomTypes_.begin(); it != hydrophilicAtomTypes_.end(); ++it)
      Out.Printf(" %s", it->Truncated().c_str());
    Out.Printf("\n");
  }
  if (!CMAP_.empty()) {
    Out.Printf("CMAP parameters:\n");
    for (CmapGridArray::const_iterator it = CMAP_.begin(); it != CMAP_.end(); ++it) {
      Out.Printf("\tCMAP %li '%s' (resolution %u) atoms {",
                 it - CMAP_.begin() + 1, it->Title().c_str(), it->Resolution());
      for (std::vector<std::string>::const_iterator an = it->AtomNames().begin();
                                                    an != it->AtomNames().end(); ++an)
        Out.Printf(" %s", an->c_str());
      Out.Printf("} residues:");
      for (std::vector<std::string>::const_iterator rn = it->ResNames().begin();
                                                    rn != it->ResNames().end(); ++rn)
        Out.Printf(" %s", rn->c_str());
      Out.Printf("\n");
      int col = 0;
      for (std::vector<double>::const_iterator gp = it->Grid().begin();
                                               gp != it->Grid().end(); ++gp)
      {
        if (col == 0)
          Out.Printf("%9.5f", *gp);
        else
          Out.Printf(" %9.5f", *gp);
        col++;
        if (col == 8) {
          Out.Printf("\n");
          col = 0;
        }
      }
      if (col != 0) Out.Printf("\n");
    }
  }
}

/** Update/add to CMAP parameters in this set with those from given set. */
int ParameterSet::updateCmapParams(CmapParmHolder const& cmap1, int debugIn, int verbose) {
  int updateCount = 0;
  for (CmapParmHolder::const_iterator newp = cmap1.begin(); newp != cmap1.end(); ++newp)
  {
    RetType ret = CMAP_.AddParm( *newp, true, debugIn );
    if (ret != ERR) {
      bool print = false;
      if (ret == ADDED) {
        if (verbose > 2) { mprintf("\tAdded NEW CMAP parameter:"); print = true; }
        updateCount++;
      } else if (ret == UPDATED) {
        if (verbose > 0) { mprintf("\tUpdated CMAP parameter:"); print = true; }
        updateCount++;
      } else if (ret == SAME) {
        if (verbose > 1) { mprintf("\tParameter for CMAP already present:"); print = true; }
      }
      if (print) {
        mprintf("\tCMAP '%s' (resolution %u) residues:", newp->Title().c_str(), newp->Resolution());
        for (std::vector<std::string>::const_iterator rn = newp->ResNames().begin();
                                                      rn != newp->ResNames().end(); ++rn)
          mprintf(" %s", rn->c_str());
        mprintf("\n");
      }
    }
  }
  return updateCount;
}

/** Update/add to parameters in this topology with those from given set. */
int ParameterSet::UpdateParamSet(ParameterSet const& set1, UpdateCount& uc, int debugIn, int verbose) {
  ParameterSet& set0 = *this;
  // Check
  if (debugIn > 0) {
    mprintf("DEBUG: Saving original parameters in originalp.dat, incoming parameters in incomingp.dat, new parameters in newp.dat.\n");
    set0.Debug("originalp.dat");
    set1.Debug("incomingp.dat");
  }

  // Bond parameters
  uc.nBondsUpdated_ = UpdateParameters< ParmHolder<BondParmType> >(set0.BP(), set1.BP(), "bond", verbose);
  // Angle parameters
  uc.nAnglesUpdated_ = UpdateParameters< ParmHolder<AngleParmType> >(set0.AP(), set1.AP(), "angle", verbose);
  // Dihedral/improper parameters
  uc.nDihedralsUpdated_ = UpdateParameters< DihedralParmHolder >(set0.DP(), set1.DP(), "dihedral", verbose);
  // Improper parameters
  uc.nImpropersUpdated_ = UpdateParameters< ImproperParmHolder >(set0.IP(), set1.IP(), "improper", verbose);
  // Urey-Bradley parameters
  uc.nUreyBradleyUpdated_ = UpdateParameters< ParmHolder<BondParmType> >(set0.UB(), set1.UB(), "Urey-Bradley", verbose);
  // Atom types
  uc.nAtomTypeUpdated_ = UpdateParameters< ParmHolder<AtomType> >(set0.AT(), set1.AT(), "atom type", verbose);
  // LJ Pairs
  uc.nLJparamsUpdated_ = UpdateParameters< ParmHolder<NonbondType> >(set0.NB(), set1.NB(), "LJ A-B", verbose);
  // LJ 1-4 Pairs
  uc.nLJ14paramsUpdated_ = UpdateParameters< ParmHolder<NonbondType> >(set0.NB14(), set1.NB14(), "LJ A-B 1-4", verbose);
  // HB LJ 10-12 Pairs
  uc.nHBparamsUpdated_ = UpdateParameters< ParmHolder<HB_ParmType> >(set0.HB(), set1.HB(), "LJ HB 10-12", verbose);
  // CMAP
  uc.nCmapUpdated_ = updateCmapParams(set1.CMAP(), debugIn, verbose);

  if (debugIn > 0) set0.Debug("newp.dat");
  return 0;
}

/** Add hydrophilic atom type. */
int ParameterSet::AddHydrophilicAtomType(NameType const& atype) {
  for (NsetType::const_iterator it = hydrophilicAtomTypes_.begin(); it != hydrophilicAtomTypes_.end(); ++it)
  {
    if (atype == *it) {
      mprintf("Warning: %s already defined as hydrophilic atom type.\n", atype.Truncated().c_str());
      return 0;
    }
  }
  // TODO check against existing types?
  hydrophilicAtomTypes_.push_back( atype );
  return 0;
}

/** \return Single string with total parameter set name */
std::string ParameterSet::ParamSetName() const {
  if (name_.size() == 1)
    return name_.front();
  else if (name_.size() > 1) {
    std::string out = name_.front();
    for (unsigned int idx = 1; idx < name_.size(); idx++)
      out.append(" + " + name_[idx]);
    return out;
  }
  return std::string();
}

/** Set parameter set name. */
void ParameterSet::SetParamSetName(std::string const& nameIn) {
  name_.push_back( nameIn );
}

/** Set nonbond parameter set name. */
void ParameterSet::SetNbParamName(std::string const& nameIn) {
  NBname_ = nameIn;
}
