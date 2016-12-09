#include "TopInfo.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // DigitWidth

/// DESTRUCTOR
TopInfo::~TopInfo() {
  if (toStdout_ && outfile_ != 0)
    delete outfile_;
}

/// CONSTRUCTOR - To Stdout
TopInfo::TopInfo(Topology* pIn) { SetupTopInfo( 0, pIn ); }

/// CONSTRUCTOR - File already set up
TopInfo::TopInfo(CpptrajFile* fIn, Topology* pIn) { SetupTopInfo( fIn, pIn ); }

// TopInfo::SetupTopInfo()
int TopInfo::SetupTopInfo(CpptrajFile* fIn, Topology* pIn) {
  if (pIn == 0) {
    mprinterr("Internal Error: TopInfo: Null topology\n");
    return 1;
  }
  parm_ = pIn;
  if (fIn == 0) {
    toStdout_ = true;
    outfile_ = new CpptrajFile();
    if (outfile_ == 0) {
      mprinterr("Internal Error: TopInfo: Could not allocate file.\n");
      return 1;
    }
    if (outfile_->OpenWrite("")) {
      delete outfile_;
      outfile_ = 0;
      mprinterr("Internal Error: TopInfo: Could not open file.\n");
      return 1;
    }
  } else {
    toStdout_ = false;
    outfile_ = fIn;
  }
  return 0;
}

int TopInfo::PrintAtomInfo(std::string const& maskExpression) const {
  AtomMask mask( maskExpression );
  if (parm_->SetupIntegerMask( mask )) return 1;
  if ( mask.None() )
    mprinterr("\tSelection is empty.\n");
  else {
    int width = DigitWidth(parm_->Natom());
    if (width < 5) width = 5;
    outfile_->Printf("%-*s %4s %*s %4s %*s %4s %8s %8s %8s %2s",
               width, "#Atom", "Name",
               width, "#Res",  "Name",
               width, "#Mol",  "Type", "Charge", "Mass", "GBradius", "El");
    if (parm_->Nonbond().HasNonbond())
      outfile_->Printf(" %8s %8s", "rVDW", "eVDW");
    outfile_->Printf("\n");
    for (AtomMask::const_iterator atnum = mask.begin(); atnum != mask.end(); atnum++) {
      const Atom& atom = (*parm_)[*atnum];
      int resnum = atom.ResNum();
      outfile_->Printf("%*i %4s %*i %4s %*i %4s %8.4f %8.4f %8.4f %2s",
                 width, *atnum+1, atom.c_str(),
                 width, resnum+1, parm_->Res(resnum).c_str(),
                 width, atom.MolNum()+1, *(atom.Type()), atom.Charge(),
                 atom.Mass(), atom.GBRadius(), atom.ElementName());
      if (parm_->Nonbond().HasNonbond())
        outfile_->Printf(" %8.4f %8.4f", parm_->GetVDWradius(*atnum), parm_->GetVDWdepth(*atnum));
      outfile_->Printf("\n");
    }
  }
  return 0;
}
