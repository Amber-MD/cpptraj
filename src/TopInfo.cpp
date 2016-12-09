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

// TopInfo::PrintAtomInfo()
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

// TopInfo::PrintBonds()
void TopInfo::PrintBonds(BondArray const& barray, BondParmArray const& bondparm,
                         CharMask const& mask1, CharMask const& mask2,
                         int& nb) const
{
  if (barray.empty()) return;
  int rwidth = DigitWidth(parm_->Nres()) + 7;
  for (BondArray::const_iterator batom = barray.begin();
                                 batom != barray.end(); ++batom)
  {
    int atom1 = batom->A1();
    int atom2 = batom->A2();
    bool printBond = false;
    if (mask2.MaskStringSet())
      printBond = (mask1.AtomInCharMask(atom1) && mask2.AtomInCharMask(atom2));
    else
      printBond = (mask1.AtomInCharMask(atom1) || mask1.AtomInCharMask(atom1));
    if (printBond) {
      outfile_->Printf("%8i:", nb);
      int bidx = batom->Idx();
      if ( bidx > -1 )
        outfile_->Printf(" %6.2f %6.3f", bondparm[bidx].Rk(), bondparm[bidx].Req());
      outfile_->Printf(" %-*s %-*s (%i,%i)",
              rwidth, parm_->AtomMaskName(atom1).c_str(),
              rwidth, parm_->AtomMaskName(atom2).c_str(),
              atom1+1, atom2+1);
      // Atom types
      const char* atype1 = *((*parm_)[atom1].Type());
      const char* atype2 = *((*parm_)[atom2].Type());
      outfile_->Printf(" %c%c-%c%c\n", atype1[0], atype1[1], atype2[0], atype2[1]);
    }
    nb++;
  }
  outfile_->Printf("\n");
}

// TopInfo::PrintBondInfo()
int TopInfo::PrintBondInfo(std::string const& mask1exp, std::string const& mask2exp) const {
  CharMask mask1( mask1exp );
  if (parm_->SetupCharMask( mask1 )) return 1;
  mprintf("#");
  mask1.MaskInfo();
  if (mask1.None()) return 1;
  CharMask mask2;
  if (!mask2exp.empty()) {
    mask2.SetMaskString( mask2exp );
    if (parm_->SetupCharMask( mask2 )) return 1;
    mprintf("#");
    mask2.MaskInfo();
    if (mask2.None()) return 1;
  }
  outfile_->Printf("#   Bond     Kb     Req       atom names   (numbers)\n");
  int nb = 1;
  PrintBonds( parm_->BondsH(), parm_->BondParm(), mask1, mask2, nb );
  PrintBonds( parm_->Bonds(),  parm_->BondParm(), mask1, mask2, nb );
  return 0;
}
