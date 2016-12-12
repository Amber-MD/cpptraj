#include "TopInfo.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // DigitWidth
#include "Constants.h" // RADDEG

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

// TopInfo::PrintResidueInfo()
int TopInfo::PrintResidueInfo(std::string const& maskExpression) const {
  AtomMask mask( maskExpression );
  if (parm_->SetupIntegerMask( mask )) return 1;
  if ( mask.None() )
    mprinterr("\tSelection is empty.\n");
  else {
    int awidth = DigitWidth(parm_->Natom());
    if (awidth < 5) awidth = 5;
    int rwidth = DigitWidth(parm_->Nres());
    if (rwidth < 5) rwidth = 5;
    int mwidth = DigitWidth(parm_->Nmol());
    if (mwidth < 5) mwidth = 5;
    outfile_->Printf("%-*s %4s %*s %*s %*s %*s %*s\n",
               rwidth, "#Res", "Name",
               awidth, "First", awidth, "Last",
               awidth, "Natom", rwidth, "#Orig", mwidth, "#Mol");
    int rn = -1;
    for (AtomMask::const_iterator atom = mask.begin();
                                  atom != mask.end(); ++atom)
    {
      if ((*parm_)[*atom].ResNum() > rn) {
        rn = (*parm_)[*atom].ResNum();
        Residue const& res = parm_->Res(rn);
        outfile_->Printf("%*i %4s %*i %*i %*i %*i %*i %c\n", rwidth, rn+1, res.c_str(),
                   awidth, res.FirstAtom()+1, awidth, res.LastAtom(),
                   awidth, res.NumAtoms(), rwidth, res.OriginalResNum(),
                   mwidth, (*parm_)[*atom].MolNum()+1, res.ChainID());
      }
    }
  }
  return 0;
}

/** Print residue info using single char names. */
int TopInfo::PrintShortResInfo(std::string const& maskString, int maxChar) const {
  AtomMask mask( maskString );
  if (parm_->SetupIntegerMask( mask )) return 1;
  if ( mask.None() )
    mprinterr("\tSelection is empty.\n");
  else {
    // Determine last selected residue.
    int max_res = (*parm_)[mask.back()].ResNum();
    int total = 0;
    int rn = -1, startRes = -1;
    std::string resLine;
    for (AtomMask::const_iterator atom = mask.begin();
                                  atom != mask.end(); ++atom)
    {
      int current_res = (*parm_)[*atom].ResNum();
      if (current_res > rn) {
        int n_res_skipped = 1;
        if (startRes == -1)
          startRes = current_res;
        else
          n_res_skipped = current_res - rn;
        // If we skipped any residues print last consective segment and start a new one.
        if (n_res_skipped > 1) {
          outfile_->Printf("%-8i %s\n", startRes+1, resLine.c_str());
          startRes = current_res;
          resLine = parm_->Res(current_res).SingleCharName();
          total = 1;
        } else {
          // Convert residue name.
          resLine += parm_->Res(current_res).SingleCharName();
          total++;
        }
        // Print if max line length reached or final res.
        if ((total%maxChar)==0 || current_res == max_res)
        {
          outfile_->Printf("%-8i %s\n", startRes+1, resLine.c_str());
          if (current_res == max_res) break;
          startRes = -1;
          resLine.clear();
        } else if ((total % 10) == 0)
          resLine += ' ';
        rn = current_res;
      }
    }
  }
  return 0;
}

// TopInfo::PrintMoleculeInfo()
int TopInfo::PrintMoleculeInfo(std::string const& maskString) const {
  if (parm_->Nmol() < 1)
    mprintf("\t'%s' No molecule info.\n", parm_->c_str());
  else {
    CharMask mask( maskString );
    if (parm_->SetupCharMask( mask )) return 1;
    if ( mask.None() )
      mprintf("\tSelection is empty.\n");
    else {
      int mwidth = DigitWidth(parm_->Nmol());
      if (mwidth < 5) mwidth = 5;
      int awidth = DigitWidth(parm_->Natom());
      if (awidth < 5) awidth = 5;
      int rwidth = DigitWidth(parm_->Nres());
      if (rwidth < 5) rwidth = 5;
      outfile_->Printf("%-*s %*s %*s %*s %*s %4s\n", mwidth, "#Mol", awidth, "Natom",
              rwidth, "Nres", rwidth, "Res0", rwidth, "Res1", "Name");
      unsigned int mnum = 1;
      for (Topology::mol_iterator mol = parm_->MolStart();
                                  mol != parm_->MolEnd(); ++mol)
      {
        if ( mask.AtomsInCharMask( mol->BeginAtom(), mol->EndAtom() ) ) {
          int firstres = (*parm_)[ mol->BeginAtom() ].ResNum();
          int lastres  = (*parm_)[ mol->EndAtom()-1 ].ResNum();
          outfile_->Printf("%*u %*i %*i %*i %*i %4s %c", mwidth, mnum, awidth, mol->NumAtoms(),
                  rwidth, lastres-firstres+1, rwidth, firstres+1,
                  rwidth, lastres+1, parm_->Res(firstres).c_str(), parm_->Res(firstres).ChainID());
          if ( mol->IsSolvent() ) outfile_->Printf(" SOLVENT");
          outfile_->Printf("\n");
        }
        ++mnum;
      }
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
      printBond = (mask1.AtomInCharMask(atom1) || mask1.AtomInCharMask(atom2));
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

// TopInfo::SetupMask()
int TopInfo::SetupMask(CharMask& maskIn) const {
  if (parm_->SetupCharMask( maskIn )) return 1;
  mprintf("#");
  maskIn.MaskInfo();
  if (maskIn.None()) return 1;
  return 0;
}

// TopInfo::SetupMask()
int TopInfo::SetupMask(std::string const& maskexp, CharMask& maskIn) const {
  if (maskexp.empty()) return 0;
  maskIn.SetMaskString( maskexp );
  return SetupMask( maskIn );
}

// TopInfo::PrintBondInfo()
int TopInfo::PrintBondInfo(std::string const& mask1exp, std::string const& mask2exp) const {
  CharMask mask1( mask1exp );
  if (SetupMask( mask1 )) return 1;
  CharMask mask2;
  if (!mask2exp.empty() && SetupMask( mask2exp, mask2 )) return 1;
  outfile_->Printf("#   Bond     Kb     Req       atom names   (numbers)\n");
  int nb = 1;
  PrintBonds( parm_->BondsH(), parm_->BondParm(), mask1, mask2, nb );
  PrintBonds( parm_->Bonds(),  parm_->BondParm(), mask1, mask2, nb );
  return 0;
}

// TopInfo::PrintAngles()
void TopInfo::PrintAngles(AngleArray const& aarray, AngleParmArray const& angleparm,
                          CharMask const& mask1, CharMask const& mask2, CharMask const& mask3,
                          int& na) const
{
  if (aarray.empty()) return;
  int rwidth = DigitWidth(parm_->Nres()) + 7;
  for (AngleArray::const_iterator aatom = aarray.begin();
                                  aatom != aarray.end(); ++aatom)
  {
    int atom1 = aatom->A1();
    int atom2 = aatom->A2();
    int atom3 = aatom->A3();
    bool printAngle = false;
    if (mask2.MaskStringSet() && mask3.MaskStringSet())
      printAngle = (mask1.AtomInCharMask(atom1) &&
                    mask2.AtomInCharMask(atom2) &&
                    mask3.AtomInCharMask(atom3));
    else
      printAngle = (mask1.AtomInCharMask(atom1) ||
                    mask1.AtomInCharMask(atom2) ||
                    mask1.AtomInCharMask(atom3));
    if (printAngle) {
      outfile_->Printf("%8i:", na);
      int aidx = aatom->Idx();
      if ( aidx > -1 )
        outfile_->Printf(" %6.3f %6.2f", angleparm[aidx].Tk(), 
                         angleparm[aidx].Teq() * Constants::RADDEG);
      outfile_->Printf(" %-*s %-*s %-*s (%i,%i,%i)",
              rwidth, parm_->AtomMaskName(atom1).c_str(),
              rwidth, parm_->AtomMaskName(atom2).c_str(),
              rwidth, parm_->AtomMaskName(atom3).c_str(),
              atom1+1, atom2+1, atom3+1);
      // Atom types
      const char* atype1 = *((*parm_)[atom1].Type());
      const char* atype2 = *((*parm_)[atom2].Type());
      const char* atype3 = *((*parm_)[atom3].Type());
      outfile_->Printf(" %c%c-%c%c-%c%c\n",atype1[0],atype1[1],atype2[0],atype2[1],
              atype3[0],atype3[1]);
    }
    na++;
  }
  outfile_->Printf("\n");
}

// TopInfo::PrintAngleInfo()
int TopInfo::PrintAngleInfo(std::string const& mask1exp, std::string const& mask2exp,
                            std::string const& mask3exp) const
{
  CharMask mask1( mask1exp );
  if (SetupMask( mask1 )) return 1;
  CharMask mask2;
  CharMask mask3;
  if (!mask2exp.empty() && SetupMask( mask2exp, mask2 )) return 1;
  if (!mask3exp.empty() && SetupMask( mask3exp, mask3 )) return 1;
  if (mask2exp.empty() != mask3exp.empty()) {
    mprinterr("Error: 2 masks specified. Require 1 or 3\n");
    return 1;
  }
  outfile_->Printf("# Angle   Kthet  degrees        atom names        (numbers)\n");
  int na = 1;
  PrintAngles( parm_->AnglesH(), parm_->AngleParm(), mask1, mask2, mask3, na );
  PrintAngles( parm_->Angles(),  parm_->AngleParm(), mask1, mask2, mask3, na );
  return 0;
}

// TopInfo::PrintDihedrals()
void TopInfo::PrintDihedrals(DihedralArray const& darray, DihedralParmArray const& dihedralparm,
                             CharMask const& mask1, CharMask const& mask2,
                             CharMask const& mask3, CharMask const& mask4, int& nd) const
{
  if (darray.empty()) return;
  int rwidth = DigitWidth(parm_->Nres()) + 7;
  for (DihedralArray::const_iterator datom = darray.begin();
                                     datom != darray.end(); ++datom)
  {
    int atom1 = datom->A1();
    int atom2 = datom->A2();
    int atom3 = datom->A3();
    int atom4 = datom->A4();
    bool printDihedral = false;
    if (mask2.MaskStringSet() && mask3.MaskStringSet() && mask4.MaskStringSet())
      printDihedral = (mask1.AtomInCharMask(atom1) &&
                       mask2.AtomInCharMask(atom2) &&
                       mask3.AtomInCharMask(atom3) &&
                       mask4.AtomInCharMask(atom4));
    else
      printDihedral = (mask1.AtomInCharMask(atom1) ||
                       mask1.AtomInCharMask(atom2) ||
                       mask1.AtomInCharMask(atom3) ||
                       mask1.AtomInCharMask(atom4));
    if (printDihedral) {
      // Determine dihedral type: 'E'nd, 'I'mproper, or 'B'oth
      char type = ' ';
      if      (datom->Type() == DihedralType::END     ) type = 'E';
      else if (datom->Type() == DihedralType::IMPROPER) type = 'I';
      else if (datom->Type() == DihedralType::BOTH    ) type = 'B';
      outfile_->Printf("%c %8i:", type, nd);
      int didx = datom->Idx();
      if ( didx > -1 )
        outfile_->Printf(" %6.3f %4.2f %4.1f",dihedralparm[didx].Pk(),dihedralparm[didx].Phase(),
                 dihedralparm[didx].Pn());
      outfile_->Printf(" %-*s %-*s %-*s %-*s (%i,%i,%i,%i)",
              rwidth, parm_->AtomMaskName(atom1).c_str(),
              rwidth, parm_->AtomMaskName(atom2).c_str(),
              rwidth, parm_->AtomMaskName(atom3).c_str(),
              rwidth, parm_->AtomMaskName(atom4).c_str(),
              atom1+1, atom2+1, atom3+1, atom4+1);
      // Atom types
      const char* atype1 = *((*parm_)[atom1].Type());
      const char* atype2 = *((*parm_)[atom2].Type());
      const char* atype3 = *((*parm_)[atom3].Type());
      const char* atype4 = *((*parm_)[atom4].Type());
      outfile_->Printf(" %c%c-%c%c-%c%c-%c%c\n",atype1[0],atype1[1],atype2[0],atype2[1],
              atype3[0],atype3[1],atype4[0],atype4[1]);
    }
    nd++;
  }
  outfile_->Printf("\n");
}

// TopInfo::PrintDihedralInfo()
int TopInfo::PrintDihedralInfo(std::string const& mask1exp, std::string const& mask2exp,
                               std::string const& mask3exp, std::string const& mask4exp) const
{
  CharMask mask1( mask1exp );
  if (SetupMask( mask1 )) return 1;
  CharMask mask2;
  CharMask mask3;
  CharMask mask4;
  if (!mask2exp.empty() && SetupMask( mask2exp, mask2 )) return 1;
  if (!mask3exp.empty() && SetupMask( mask3exp, mask3 )) return 1;
  if (!mask4exp.empty() && SetupMask( mask4exp, mask4 )) return 1;
  if (mask2exp.empty() != mask3exp.empty() || mask2exp.empty() != mask4exp.empty()) {
    mprinterr("Error: Require either 1 mask or 4 masks\n");
    return 1;
  }
  outfile_->Printf("#Dihedral    pk     phase pn                atoms\n");
  int nd = 1;
  PrintDihedrals( parm_->DihedralsH(), parm_->DihedralParm(), mask1, mask2, mask3, mask4, nd );
  PrintDihedrals( parm_->Dihedrals(),  parm_->DihedralParm(), mask1, mask2, mask3, mask4, nd );
  return 0;
}
