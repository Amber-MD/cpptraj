#include <algorithm> // std::min, std::max
#include "TopInfo.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // DigitWidth
#include "Constants.h" // RADDEG
#include "DistRoutines.h" // DIST_NoImage
#include "TorsionRoutines.h" // CalcAngle, Torsion

/// DESTRUCTOR
TopInfo::~TopInfo() {
  if (toStdout_ && outfile_ != 0)
    delete outfile_;
}

/// CONSTRUCTOR - To Stdout
TopInfo::TopInfo(Topology* pIn) { SetupTopInfo( 0, pIn, 0 ); }

// TopInfo::SetupTopInfo()
int TopInfo::SetupTopInfo(CpptrajFile* fIn, Topology* pIn, DataSet_Coords* cIn) {
  if (cIn == 0 && pIn == 0) {
    mprinterr("Internal Error: TopInfo: Null topology\n");
    return 1;
  }
  if (cIn != 0) {
    parm_ = cIn->TopPtr();
    coords_ = cIn->AllocateFrame();
    cIn->GetFrame(0, coords_);
  } else
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
  Awidth_ = std::max(2, DigitWidth(parm_->Natom()));
  Rwidth_ = DigitWidth(parm_->Nres()) + 6; // :, @, 4 char atom name
  max_type_len_ = 2;
  for (int i = 0; i != parm_->Natom(); i++)
    max_type_len_ = std::max( max_type_len_, (*parm_)[i].Type().len() );
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
    int awidth = std::max(5, DigitWidth(parm_->Natom()));
    int rwidth = std::max(5, DigitWidth(parm_->Nres()));
    int mwidth = std::max(5, DigitWidth(parm_->Nmol()));
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
      int awidth = std::max(5, DigitWidth(parm_->Natom()));
      int rwidth = std::max(5, DigitWidth(parm_->Nres()));
      int mwidth = std::max(5, DigitWidth(parm_->Nmol()));
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

/** \return An array of unique molecule types selected by mask. */
TopInfo::Marray TopInfo::UniqueMolCount(Topology const& top, CharMask const& mask) const
{
  Marray mols;
  typedef std::vector<int> Iarray;
  // Hold begin residue for each molecule type
  Iarray Res0;
  // Hold end residue for each molecule type
  Iarray Res1;
  for (Topology::mol_iterator CurMol = top.MolStart();
                              CurMol != top.MolEnd(); ++CurMol)
  {
    if ( mask.AtomsInCharMask( CurMol->BeginAtom(), CurMol->EndAtom() ) ) {
      // Has this molecule been seen before?
      int natom = CurMol->NumAtoms();
      int res0 = top[ CurMol->BeginAtom() ].ResNum(); // CurMol first residue
      int res1 = top[ CurMol->EndAtom()-1 ].ResNum(); // CurMol last residue
      int nres = res1 - res0 + 1;
      int matchIdx = -1;
      for (int idx = 0; idx != (int)mols.size(); idx++)
      {
        // First check number of atoms, then number residues, then residue names.
        if ( mols[idx].natom_ == natom ) {
          if ( mols[idx].nres_ == nres ) {
            matchIdx = idx;
            int cridx = res0; // current molecule residue index
            for (int rridx = Res0[idx]; rridx <= Res1[idx]; rridx++, cridx++)
            {
              if ( top.Res(rridx).Name() != top.Res(cridx).Name() ) {
                // Residue name mismatch.
                matchIdx = -1;
                break;
              }
            } // END loop over all residues
          }
        }
        if (matchIdx != -1) break;
      } // END loop over found molecule types
      if (matchIdx == -1) {
        // New molecule
        mols.push_back(
          MolType(CurMol-top.MolStart(), natom, nres, top.Res(res0).Name().Truncated()));
        Res0.push_back( res0 );
        Res1.push_back( res1 );
      } else {
        // Existing molecule. Update count.
        mols[matchIdx].UpdateCount();
      }
    } // END molecule in mask
  } // END loop over all molecules
  return mols;
} 

// TopInfo::PrintShortMolInfo()
int TopInfo::PrintShortMolInfo(std::string const& maskString) const {
  if (parm_->Nmol() < 1)
    mprintf("\t'%s' No molecule info.\n", parm_->c_str());
  else {
    CharMask mask( maskString );
    if (parm_->SetupCharMask( mask )) return 1;
    if ( mask.None() )
      mprintf("\tSelection is empty.\n");
    else {
      Marray mols = MolCount(*parm_, mask);
      // Determine max counts for nice output formatting
      int maxNatom = 0;
      int maxNres = 0;
      int maxCount = 0;
      for (Marray::const_iterator mol = mols.begin(); mol != mols.end(); ++mol) {
        maxNatom = std::max( maxNatom, mol->natom_  );
        maxNres  = std::max( maxNres,  mol->nres_   );
        maxCount = std::max( maxCount, mol->count_  );
      }
      int awidth = std::max(5, DigitWidth(maxNatom));
      int rwidth = std::max(5, DigitWidth(maxNres ));
      int mwidth = std::max(5, DigitWidth(maxCount));
      outfile_->Printf("%-4s %*s %*s %*s\n", "#Mol", mwidth, "Count", 
                       awidth, "Natom", rwidth, "Nres");
      for (Marray::const_iterator mol = mols.begin(); mol != mols.end(); ++mol)
        outfile_->Printf("%-4s %*i %*i %*i\n", mol->name_.c_str(),
                         mwidth, mol->count_,
                         awidth, mol->natom_,
                         rwidth, mol->nres_);
    } // END selection not empty
  } // END parm has molecules
  return 0;
}

// TopInfo::PrintChargeInfo()
int TopInfo::PrintChargeInfo(std::string const& maskExpression) const {
  AtomMask mask( maskExpression );
  if (parm_->SetupIntegerMask( mask )) return 1;
  double sumQ = 0.0;
  for (AtomMask::const_iterator idx = mask.begin(); idx != mask.end(); ++idx)
    sumQ += (*parm_)[*idx].Charge();
  outfile_->Printf("Sum of charges in mask [%s](%i) is %g\n",
                   mask.MaskString(), mask.Nselected(), sumQ);
  return 0;
}

// TopInfo::PrintMassInfo()
int TopInfo::PrintMassInfo(std::string const& maskExpression) const {
  AtomMask mask( maskExpression );
  if (parm_->SetupIntegerMask( mask )) return 1;
  double sumM = 0.0;
  for (AtomMask::const_iterator idx = mask.begin(); idx != mask.end(); ++idx)
    sumM += (*parm_)[*idx].Mass();
  outfile_->Printf("Sum of masses in mask [%s](%i) is %g\n",
                   mask.MaskString(), mask.Nselected(), sumM);
  return 0;
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


// TopInfo::PrintBonds()
void TopInfo::PrintBonds(BondArray const& barray, BondParmArray const& bondparm,
                         CharMask const& mask1, CharMask const& mask2,
                         int nw, int& nb) const
{
  if (barray.empty()) return;
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
      outfile_->Printf("%*i", nw, nb);
      int bidx = batom->Idx();
      if ( bidx > -1 )
        outfile_->Printf(" %6.2f %6.3f", bondparm[bidx].Rk(), bondparm[bidx].Req());
      if ( !coords_.empty() )
        outfile_->Printf(" %6.3f", DIST_NoImage(coords_.XYZ(atom1), coords_.XYZ(atom2)));
      outfile_->Printf(" %-*s %-*s %*i %*i",
              Rwidth_, parm_->AtomMaskName(atom1).c_str(),
              Rwidth_, parm_->AtomMaskName(atom2).c_str(),
              Awidth_, atom1+1, Awidth_, atom2+1);
      // Atom types
      outfile_->Printf(" %*s %*s\n",
                       max_type_len_, (*parm_)[atom1].Type().Truncated().c_str(),
                       max_type_len_, (*parm_)[atom2].Type().Truncated().c_str());
    }
    nb++;
  }
  //outfile_->Printf("\n");
}

// TopInfo::PrintBondInfo()
int TopInfo::PrintBondInfo(std::string const& mask1exp, std::string const& mask2exp) const {
  CharMask mask1( mask1exp );
  if (SetupMask( mask1 )) return 1;
  CharMask mask2;
  if (!mask2exp.empty() && SetupMask( mask2exp, mask2 )) return 1;
  int nw = std::max(4, DigitWidth(parm_->BondsH().size() + parm_->Bonds().size()));
  outfile_->Printf("%-*s", nw, "#Bnd");
  if (!parm_->BondParm().empty())
    outfile_->Printf(" %6s %6s", "RK", "REQ");
  if (!coords_.empty())
    outfile_->Printf(" %6s", "Value");
  outfile_->Printf(" %-*s %-*s %*s %*s %*s %*s\n",
                   Rwidth_, "Atom1", Rwidth_, "Atom2",
                   Awidth_, "A1", Awidth_, "A2",
                   max_type_len_, "T1", max_type_len_, "T2");
  int nb = 1;
  PrintBonds( parm_->BondsH(), parm_->BondParm(), mask1, mask2, nw, nb );
  PrintBonds( parm_->Bonds(),  parm_->BondParm(), mask1, mask2, nw, nb );
  return 0;
}

// TopInfo::PrintAngles()
void TopInfo::PrintAngles(AngleArray const& aarray, AngleParmArray const& angleparm,
                          CharMask const& mask1, CharMask const& mask2, CharMask const& mask3,
                          int nw, int& na) const
{
  if (aarray.empty()) return;
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
      outfile_->Printf("%*i", nw, na);
      int aidx = aatom->Idx();
      if ( aidx > -1 )
        outfile_->Printf(" %6.3f %6.2f", angleparm[aidx].Tk(), 
                         angleparm[aidx].Teq() * Constants::RADDEG);
      if ( !coords_.empty() )
        outfile_->Printf(" %6.2f", CalcAngle(coords_.XYZ(atom1),
                                             coords_.XYZ(atom2),
                                             coords_.XYZ(atom3)) * Constants::RADDEG);
      outfile_->Printf(" %-*s %-*s %-*s %*i %*i %*i",
              Rwidth_, parm_->AtomMaskName(atom1).c_str(),
              Rwidth_, parm_->AtomMaskName(atom2).c_str(),
              Rwidth_, parm_->AtomMaskName(atom3).c_str(),
              Awidth_, atom1+1, Awidth_, atom2+1, Awidth_, atom3+1);
      // Atom types
      outfile_->Printf(" %*s %*s %*s\n",
                       max_type_len_, (*parm_)[atom1].Type().Truncated().c_str(),
                       max_type_len_, (*parm_)[atom2].Type().Truncated().c_str(),
                       max_type_len_, (*parm_)[atom3].Type().Truncated().c_str());
    }
    na++;
  }
  //outfile_->Printf("\n");
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
    mprinterr("Error: Require either 1 mask or 3 masks.\n");
    return 1;
  }
  int nw = std::max(4, DigitWidth(parm_->AnglesH().size() + parm_->Angles().size()));
  outfile_->Printf("%-*s", nw, "#Ang");
  if (!parm_->AngleParm().empty())
    outfile_->Printf(" %6s %6s", "TK", "TEQ");
  if (!coords_.empty())
    outfile_->Printf(" %6s", "Value");
  outfile_->Printf(" %-*s %-*s %-*s %*s %*s %*s %*s %*s %*s\n",
                   Rwidth_, "Atom1", Rwidth_, "Atom2", Rwidth_, "Atom3",
                   Awidth_, "A1", Awidth_, "A2", Awidth_, "A3",
                   max_type_len_, "T1", max_type_len_, "T2", max_type_len_, "T3");
  int na = 1;
  PrintAngles( parm_->AnglesH(), parm_->AngleParm(), mask1, mask2, mask3, nw, na );
  PrintAngles( parm_->Angles(),  parm_->AngleParm(), mask1, mask2, mask3, nw, na );
  return 0;
}

// TopInfo::PrintDihedrals()
void TopInfo::PrintDihedrals(DihedralArray const& darray, DihedralParmArray const& dihedralparm,
                             CharMask const& mask1, CharMask const& mask2,
                             CharMask const& mask3, CharMask const& mask4, int nw, int& nd) const
{
  if (darray.empty()) return;
  for (DihedralArray::const_iterator dih = darray.begin();
                                     dih != darray.end(); ++dih)
  {
    int atom1 = dih->A1();
    int atom2 = dih->A2();
    int atom3 = dih->A3();
    int atom4 = dih->A4();
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
      if      (dih->Type() == DihedralType::END     ) type = 'E';
      else if (dih->Type() == DihedralType::IMPROPER) type = 'I';
      else if (dih->Type() == DihedralType::BOTH    ) type = 'B';
      outfile_->Printf("%c %*i", type, nw, nd);
      int didx = dih->Idx();
      if ( didx > -1 )
        outfile_->Printf(" %6.3f %5.2f %4.1f",dihedralparm[didx].Pk(),dihedralparm[didx].Phase(),
                 dihedralparm[didx].Pn());
      if ( !coords_.empty() )
        outfile_->Printf(" %7.2f", Torsion( coords_.XYZ(atom1),
                                            coords_.XYZ(atom2),
                                            coords_.XYZ(atom3),
                                            coords_.XYZ(atom4) ) * Constants::RADDEG );
      outfile_->Printf(" %-*s %-*s %-*s %-*s %*i %*i %*i %*i",
              Rwidth_, parm_->AtomMaskName(atom1).c_str(),
              Rwidth_, parm_->AtomMaskName(atom2).c_str(),
              Rwidth_, parm_->AtomMaskName(atom3).c_str(),
              Rwidth_, parm_->AtomMaskName(atom4).c_str(),
              Awidth_, atom1+1, Awidth_, atom2+1,
              Awidth_, atom3+1, Awidth_, atom4+1);
      // Atom types
      outfile_->Printf(" %*s %*s %*s %*s\n",
                       max_type_len_, (*parm_)[atom1].Type().Truncated().c_str(),
                       max_type_len_, (*parm_)[atom2].Type().Truncated().c_str(),
                       max_type_len_, (*parm_)[atom3].Type().Truncated().c_str(),
                       max_type_len_, (*parm_)[atom4].Type().Truncated().c_str());
    }
    nd++;
  }
  //outfile_->Printf("\n");
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
    mprinterr("Error: Require either 1 mask or 4 masks.\n");
    return 1;
  }
  int nw = std::max(3, DigitWidth(parm_->DihedralsH().size() + parm_->Dihedrals().size()));
  outfile_->Printf("# %*s", nw, "Dih");
  if (!parm_->DihedralParm().empty())
    outfile_->Printf(" %6s %5s %4s", "PK", "Phase", "PN");
  if (!coords_.empty())
    outfile_->Printf(" %7s", "Value");
  outfile_->Printf(" %-*s %-*s %-*s %-*s %*s %*s %*s %*s %*s %*s %*s %*s\n",
                   Rwidth_, "Atom1", Rwidth_, "Atom2",
                   Rwidth_, "Atom3", Rwidth_, "Atom4",
                   Awidth_, "A1", Awidth_, "A2",
                   Awidth_, "A3", Awidth_, "A4",
                   max_type_len_, "T1", max_type_len_, "T2",
                   max_type_len_, "T3", max_type_len_, "T4");
  int nd = 1;
  PrintDihedrals( parm_->DihedralsH(), parm_->DihedralParm(), mask1, mask2, mask3, mask4, nw, nd );
  PrintDihedrals( parm_->Dihedrals(),  parm_->DihedralParm(), mask1, mask2, mask3, mask4, nw, nd );
  return 0;
}
