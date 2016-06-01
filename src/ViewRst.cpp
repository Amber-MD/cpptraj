#include "ViewRst.h"
#include "CpptrajStdio.h"
#include "Trajout_Single.h"
#include "ParmFile.h"

// ViewRst::Init()
int ViewRst::Init(Topology const& topIn, OutputType typeIn)
{

  outType_ = typeIn;
  int Ntops = 1;
  if (outType_ == BY_STRENGTH)
    Ntops = 4; // strong, med, weak, v. weak
  // Add all atoms to each pseudo topology.
  Pseudo_.clear();
  Pseudo_.resize( Ntops );
  for (int nt = 0; nt != Ntops; nt++) {
    for (Topology::atom_iterator atm = topIn.begin(); atm != topIn.end(); ++atm)
      Pseudo_[nt].AddTopAtom( *atm, topIn.Res( atm->ResNum() ) );
  }
  return 0;
}

// ViewRst::AddRst()
void ViewRst::AddRst(int atom_i, int atom_j, NoeType strength) {
  int ntop = (int)strength;
  Pseudo_[ntop].AddBond(atom_i, atom_j);
}

// ViewRst::GenerateOutNames()
std::vector<FileName> ViewRst::GenerateOutNames(FileName const& fname) const {
  std::vector<FileName> OutNames;
  switch (outType_) {
    case BY_STRENGTH:
      OutNames.push_back( fname.PrependFileName("strong."  ) );
      OutNames.push_back( fname.PrependFileName("medium."  ) );
      OutNames.push_back( fname.PrependFileName("weak."    ) );
      OutNames.push_back( fname.PrependFileName("veryweak.") );
      break;
    case ALL:
      OutNames.push_back( fname );
      break;
  }
  return OutNames;
}

// ViewRst::WriteRstMol2()
// FIXME: Not const since PrepareTrajWrite needs Topology*
int ViewRst::WriteRstMol2(std::string const& mol2out, Frame const& frameIn) {
  if (Pseudo_.empty()) return 0;
  if (mol2out.empty()) {
    mprinterr("Internal Error: No mol2 output name given.\n");
    return 1;
  }
  std::vector< FileName > OutNames = GenerateOutNames(mol2out);
  // Check coords
  if (Pseudo_[0].Natom() != frameIn.Natom()) {
    mprinterr("Internal Error: Number of topology atoms (%i) != number frame atoms (%i)\n",
              Pseudo_[0].Natom(), frameIn.Natom());
    return 1;
  }
  if (frameIn.empty()) {
    mprinterr("Internal Error: Input frame is empty.\n");
    return 1;
  }

  for (unsigned int nt = 0; nt != Pseudo_.size(); nt++) {
    Trajout_Single trajout;
    if (trajout.PrepareTrajWrite(OutNames[nt], ArgList(), &(Pseudo_[nt]), CoordinateInfo(), 1,
                                 TrajectoryFile::MOL2FILE))
      return 1;
    if (trajout.WriteSingle(0, frameIn)) return 1;
    trajout.EndTraj();
  }

  return 0;
}

int ViewRst::WriteRstTop(std::string const& topOut) {
  if (Pseudo_.empty()) return 0;
  if (topOut.empty()) {
    mprinterr("Internal Error: No topology output name given.\n");
    return 1;
  }
  std::vector< FileName > OutNames = GenerateOutNames(topOut);

  for (unsigned int nt = 0; nt != Pseudo_.size(); nt++) {
    ParmFile topOut;
    Pseudo_[nt].CommonSetup(false);
    if (topOut.WriteTopology(Pseudo_[nt], OutNames[nt], ArgList(), ParmFile::AMBERPARM, 0))
      return 1;
  }
  return 0;
}
