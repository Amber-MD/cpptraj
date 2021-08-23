#include <algorithm> // std::min, std::max
#include "TopInfo.h"
#include "CpptrajFile.h"
#include "DataSet_Coords.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // DigitWidth
#include "Constants.h" // RADDEG
#include "DistRoutines.h" // DIST_NoImage
#include "TorsionRoutines.h" // CalcAngle, Torsion
#include "Mol.h"
#include "CharMask.h"

/// CONSTRUCTOR
TopInfo::TopInfo() :
  outfile_(0),
  parm_(0),
  Awidth_(0),
  amn_width_(0),
  max_aname_len_(0),
  toStdout_(false),
  noIntraRes_(false)
{}

/// CONSTRUCTOR - To Stdout
TopInfo::TopInfo(Topology const* pIn) :
  outfile_(0),
  parm_(0),
  Awidth_(0),
  amn_width_(0),
  max_type_len_(0),
  max_aname_len_(0),
  toStdout_(false),
  noIntraRes_(false)
{
  SetupTopInfo( 0, pIn, 0 );
}

/// DESTRUCTOR
TopInfo::~TopInfo() {
  if (toStdout_ && outfile_ != 0)
    delete outfile_;
}

// TopInfo::SetupTopInfo()
int TopInfo::SetupTopInfo(CpptrajFile* fIn, Topology const* pIn, DataSet_Coords* cIn) {
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
  // Determine widths for output.
  Awidth_ = std::max(2, DigitWidth(parm_->Natom()));
  max_type_len_ = 2;
  max_aname_len_ = 4;
  for (int i = 0; i != parm_->Natom(); i++) {
    max_type_len_ = std::max( max_type_len_, (*parm_)[i].Type().len() );
    max_aname_len_ = std::max( max_aname_len_, (*parm_)[i].Name().len() );
  }
  amn_width_ = DigitWidth(parm_->Nres()) + max_aname_len_ + 2; // :, @
  return 0;
}

/** \return maximum atom/atom type name length from selection. */
int TopInfo::maxAtomNamesWidth(AtomMask const& mask) const {
  // Sanity check.
  if (parm_ == 0) {
    mprinterr("Internal Error: TopInfo::maxAtomNamesWidth: parm is null.\n");
    return 0;
  }
  int nWidth = 4;
  for (AtomMask::const_iterator atnum = mask.begin(); atnum != mask.end(); ++atnum)
  {
    Atom const& at = (*parm_)[*atnum];
    int an_size = at.Name().len();
    int at_size = at.Type().len();
    if (an_size > nWidth) nWidth = an_size;
    if (at_size > nWidth) nWidth = at_size;
  }
  return nWidth;
}

// TopInfo::PrintAtomInfo()
int TopInfo::PrintAtomInfo(std::string const& maskExpression) const {
  if (maskExpression.empty()) {
    mprinterr("Error: No valid mask given to select atoms.\n");
    return 1;
  }
  AtomMask mask( maskExpression );
  if (parm_->SetupIntegerMask( mask )) return 1;
  if ( mask.None() )
    mprinterr("\tSelection is empty.\n");
  else {
    mprintf("%i atoms selected.\n", mask.Nselected());
    int width = DigitWidth(parm_->Natom());
    if (width < 5) width = 5;
    int nWidth = maxAtomNamesWidth(mask);
    outfile_->Printf("%-*s %-*s %*s %-*s %*s %-*s %8s %8s %8s %2s",
               width, "#Atom",
               nWidth, "Name",
               width, "#Res",
               nWidth, "Name",
               width, "#Mol",
               nWidth, "Type",
               "Charge", "Mass", "GBradius", "El");
    if (parm_->Nonbond().HasNonbond())
      outfile_->Printf(" %8s %8s", "rVDW", "eVDW");
    outfile_->Printf("\n");
    for (AtomMask::const_iterator atnum = mask.begin(); atnum != mask.end(); ++atnum) {
      const Atom& atom = (*parm_)[*atnum];
      int resnum = atom.ResNum();
      outfile_->Printf("%*i %-*s %*i %-*s %*i %-*s %8.4f %8.4f %8.4f %2s",
                 width, *atnum+1,
                 nWidth, atom.c_str(),
                 width, resnum+1,
                 nWidth, parm_->Res(resnum).c_str(),
                 width, atom.MolNum()+1,
                 nWidth, *(atom.Type()),
                 atom.Charge(), atom.Mass(), atom.GBRadius(), atom.ElementName());
      if (parm_->Nonbond().HasNonbond())
        outfile_->Printf(" %8.4f %8.4f", parm_->GetVDWradius(*atnum), parm_->GetVDWdepth(*atnum));
      outfile_->Printf("\n");
    }
  }
  return 0;
}

/** \return max residue name length from selection. */
int TopInfo::maxResNameWidth(std::vector<int> const& resNums) const {
  int nWidth = 4;
  for (std::vector<int>::const_iterator rnum = resNums.begin(); rnum != resNums.end(); ++rnum)
  {
    int rn_size = parm_->Res(*rnum).Name().len();
    if (rn_size > nWidth)
      nWidth = rn_size;
  }
  return nWidth;
}

// TopInfo::PrintResidueInfo()
int TopInfo::PrintResidueInfo(std::string const& maskExpression) const {
  if (maskExpression.empty()) {
    mprinterr("Error: No valid mask given to select residues.\n");
    return 1;
  }
  AtomMask mask( maskExpression );
  if (parm_->SetupIntegerMask( mask )) return 1;
  if ( mask.None() )
    mprinterr("\tSelection is empty.\n");
  else {
    std::vector<int> resNums = parm_->ResnumsSelectedBy(mask);
    mprintf("%zu residues selected.\n", resNums.size());
    int rn_width = maxResNameWidth( resNums );

    int awidth = std::max(5, DigitWidth(parm_->Natom()));
    int rwidth = std::max(5, DigitWidth(parm_->Nres()));
    int mwidth = std::max(5, DigitWidth(parm_->Nmol()));
    outfile_->Printf("%-*s %-*s %*s %*s %*s %*s %*s %c %c\n",
               rwidth, "#Res",
               rn_width, "Name",
               awidth, "First",
               awidth, "Last",
               awidth, "Natom",
               rwidth, "#Orig",
               mwidth, "#Mol",
               'C', 'I');
    for (std::vector<int>::const_iterator rnum = resNums.begin(); rnum != resNums.end(); ++rnum)
    {
      Residue const& res = parm_->Res(*rnum);
      outfile_->Printf("%*i %-*s %*i %*i %*i %*i %*i %c %c\n",
                       rwidth, *rnum + 1,
                       rn_width, res.c_str(),
                       awidth, res.FirstAtom()+1,
                       awidth, res.LastAtom(),
                       awidth, res.NumAtoms(),
                       rwidth, res.OriginalResNum(),
                       mwidth, (*parm_)[res.FirstAtom()].MolNum()+1,
                       res.ChainId(), res.Icode());
    }
  }
  return 0;
}

/** Print residue info using single char names. */ // TODO use Topology::ResnumsSelectedBy
int TopInfo::PrintShortResInfo(std::string const& maskString, int maxChar) const {
  if (maskString.empty()) {
    mprinterr("Error: No valid mask given for short residue info.\n");
    return 1;
  }
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

/** \return max molecule name length from selection. */
int TopInfo::maxMolNameWidth(std::vector<int> const& molNums) const {
  int nWidth = 4;
  for (std::vector<int>::const_iterator mnum = molNums.begin(); mnum != molNums.end(); ++mnum)
  {
    int molAtom0 = parm_->Mol( *mnum ).MolUnit().Front();
    int firstres = (*parm_)[ molAtom0 ].ResNum();
    int mn_size = parm_->Res(firstres).Name().len();
    if (mn_size > nWidth)
      nWidth = mn_size;
  }
  return nWidth;
}

// TopInfo::PrintMoleculeInfo()
int TopInfo::PrintMoleculeInfo(std::string const& maskString) const {
  if (maskString.empty()) {
    mprinterr("Error: No valid mask given to select molecules.\n");
    return 1;
  }
  if (parm_->Nmol() < 1)
    mprintf("\t'%s' No molecule info.\n", parm_->c_str());
  else {
    AtomMask mask( maskString );
    if (parm_->SetupIntegerMask( mask )) return 1;
    if ( mask.None() )
      mprintf("\tSelection is empty.\n");
    else {
      std::vector<int> molNums = parm_->MolnumsSelectedBy( mask );
      mprintf("%zu molecules.\n", molNums.size());
      // TODO determine max segments
      int mn_width = maxMolNameWidth( molNums );
      int awidth = std::max(5, DigitWidth(parm_->Natom()));
      int rwidth = std::max(5, DigitWidth(parm_->Nres()));
      int mwidth = std::max(5, DigitWidth(parm_->Nmol()));
      outfile_->Printf("%-*s %*s %*s %*s %*s %-*s %c\n",
                       mwidth, "#Mol",
                       awidth, "Natom",
                       rwidth, "Nres",
                       rwidth, "Res0",
                       rwidth, "Res1", 
                       mn_width, "Name",
                       'C');
      for (std::vector<int>::const_iterator mnum = molNums.begin();
                                            mnum != molNums.end(); ++mnum)
      {
        Molecule const& Mol = parm_->Mol(*mnum);
        outfile_->Printf("%*u %*i", mwidth, *mnum+1, awidth, Mol.NumAtoms());
        // Loop over segments
        for (Unit::const_iterator seg = Mol.MolUnit().segBegin();
                                  seg != Mol.MolUnit().segEnd(); ++seg)
        {
          int firstres = (*parm_)[ seg->Begin() ].ResNum();
          int lastres  = (*parm_)[ seg->End()-1 ].ResNum();
          outfile_->Printf(" %*i %*i %*i %-*s %c",
                           rwidth, lastres-firstres+1,
                           rwidth, firstres+1,
                           rwidth, lastres+1,
                           mn_width, parm_->Res(firstres).c_str(),
                           parm_->Res(firstres).ChainId());
        } // END loop over segments
        if ( Mol.IsSolvent() ) outfile_->Printf(" SOLVENT");
        outfile_->Printf("\n");
      }
    }
  }
  return 0;
}

// TopInfo::PrintShortMolInfo()
int TopInfo::PrintShortMolInfo(std::string const& maskString) const {
  if (maskString.empty()) {
    mprinterr("Error: No valid mask given for short molecule info.\n");
    return 1;
  }
  if (parm_->Nmol() < 1)
    mprintf("\t'%s' No molecule info.\n", parm_->c_str());
  else {
    AtomMask mask( maskString );
    if (parm_->SetupIntegerMask( mask )) return 1;
    if ( mask.None() )
      mprintf("\tSelection is empty.\n");
    else {
      std::vector<int> molNums = parm_->MolnumsSelectedBy( mask );
      Mol::Marray mols = Mol::UniqueCount(*parm_, molNums);
      // Determine max counts for nice output formatting
      int maxNatom = 0;
      int maxNres = 0;
      unsigned int maxCount = 0;
      unsigned int mn_width = 4;
      for (Mol::Marray::const_iterator mol = mols.begin(); mol != mols.end(); ++mol) {
        maxNatom = std::max( maxNatom, mol->natom_ );
        maxNres  = std::max( maxNres,  mol->nres_  );
        maxCount = std::max( maxCount, (unsigned int)mol->idxs_.size() );
        mn_width = std::max( mn_width, (unsigned int)mol->name_.size() );
      }
      int awidth = std::max(5, DigitWidth(maxNatom));
      int rwidth = std::max(5, DigitWidth(maxNres ));
      int mwidth = std::max(5, DigitWidth(maxCount));
      outfile_->Printf("%-*s %*s %*s %*s\n",
                       mn_width, "#Mol",
                       mwidth, "Count",
                       awidth, "Natom",
                       rwidth, "Nres");
      for (Mol::Marray::const_iterator mol = mols.begin(); mol != mols.end(); ++mol)
        outfile_->Printf("%-*s %*zu %*i %*i\n",
                         mn_width, mol->name_.c_str(),
                         mwidth, mol->idxs_.size(),
                         awidth, mol->natom_,
                         rwidth, mol->nres_);
    } // END selection not empty
  } // END parm has molecules
  return 0;
}

// TopInfo::PrintChargeInfo()
int TopInfo::PrintChargeInfo(std::string const& maskExpression, double& sumQ) const {
  if (maskExpression.empty()) {
    mprinterr("Error: No valid mask given to calculate charge.\n");
    return 1;
  }
  AtomMask mask( maskExpression );
  if (parm_->SetupIntegerMask( mask )) return 1;
  sumQ = 0.0;
  for (AtomMask::const_iterator idx = mask.begin(); idx != mask.end(); ++idx)
    sumQ += (*parm_)[*idx].Charge();
  outfile_->Printf("Sum of charges in mask [%s](%i) is %g\n",
                   mask.MaskString(), mask.Nselected(), sumQ);
  return 0;
}

// TopInfo::PrintMassInfo()
int TopInfo::PrintMassInfo(std::string const& maskExpression, double& sumM) const {
  if (maskExpression.empty()) {
    mprinterr("Error: No valid mask given to calculate mass.\n");
    return 1;
  }
  AtomMask mask( maskExpression );
  if (parm_->SetupIntegerMask( mask )) return 1;
  sumM = 0.0;
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
  if (maskIn.SetMaskString( maskexp )) return 1;
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
    if (noIntraRes_ && printBond) {
      if ( (*parm_)[atom1].ResNum() == (*parm_)[atom2].ResNum() )
        printBond = false;
    }
    if (printBond) {
      outfile_->Printf("%*i", nw, nb);
      int bidx = batom->Idx();
      if ( bidx > -1 )
        outfile_->Printf(" %6.2f %6.3f", bondparm[bidx].Rk(), bondparm[bidx].Req());
      if ( !coords_.empty() )
        outfile_->Printf(" %6.3f", DIST_NoImage(coords_.XYZ(atom1), coords_.XYZ(atom2)));
      outfile_->Printf(" %-*s %-*s %*i %*i",
              amn_width_, parm_->AtomMaskName(atom1).c_str(),
              amn_width_, parm_->AtomMaskName(atom2).c_str(),
              Awidth_, atom1+1, Awidth_, atom2+1);
      // Atom types
      outfile_->Printf(" %*s %*s\n",
                       max_type_len_, *((*parm_)[atom1].Type()),
                       max_type_len_, *((*parm_)[atom2].Type()));
    }
    nb++;
  }
  //outfile_->Printf("\n");
}

// TopInfo::PrintBondInfo()
int TopInfo::PrintBondInfo(std::string const& mask1exp, std::string const& mask2exp,
                           bool printUB) const
{
  if (printUB && !parm_->Chamber().HasChamber()) {
    mprintf("Warning: '%s' does not have any CHARMM parameters.\n", parm_->c_str());
    return 0;
  }
  CharMask mask1( mask1exp );
  if (SetupMask( mask1 )) return 1;
  CharMask mask2;
  if (!mask2exp.empty() && SetupMask( mask2exp, mask2 )) return 1;
  size_t n_bonds;
  const char* label = "#Bnd";
  bool hasParams;
  if (printUB) {
    label = "#UB";
    n_bonds = parm_->Chamber().UB().size();
    hasParams = !parm_->Chamber().UBparm().empty();
  } else {
    n_bonds = parm_->BondsH().size() + parm_->Bonds().size();
    hasParams = !parm_->BondParm().empty();
  }
  int nw = std::max(4, DigitWidth(n_bonds));
  outfile_->Printf("%-*s", nw, label);
  if (hasParams)
    outfile_->Printf(" %6s %6s", "RK", "REQ");
  if (!coords_.empty())
    outfile_->Printf(" %6s", "Value");
  outfile_->Printf(" %-*s %-*s %*s %*s %*s %*s\n",
                   amn_width_, "Atom1", amn_width_, "Atom2",
                   Awidth_, "A1", Awidth_, "A2",
                   max_type_len_, "T1", max_type_len_, "T2");
  int nb = 1;
  if (printUB)
    PrintBonds( parm_->Chamber().UB(), parm_->Chamber().UBparm(), mask1, mask2, nw, nb );
  else {
    PrintBonds( parm_->BondsH(), parm_->BondParm(), mask1, mask2, nw, nb );
    PrintBonds( parm_->Bonds(),  parm_->BondParm(), mask1, mask2, nw, nb );
  }
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
              amn_width_, parm_->AtomMaskName(atom1).c_str(),
              amn_width_, parm_->AtomMaskName(atom2).c_str(),
              amn_width_, parm_->AtomMaskName(atom3).c_str(),
              Awidth_, atom1+1, Awidth_, atom2+1, Awidth_, atom3+1);
      // Atom types
      outfile_->Printf(" %*s %*s %*s\n",
                       max_type_len_, *((*parm_)[atom1].Type()),
                       max_type_len_, *((*parm_)[atom2].Type()),
                       max_type_len_, *((*parm_)[atom3].Type()));
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
                   amn_width_, "Atom1", amn_width_, "Atom2", amn_width_, "Atom3",
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
        outfile_->Printf(" %7.3f %5.2f %4.1f",dihedralparm[didx].Pk(),dihedralparm[didx].Phase(),
                 dihedralparm[didx].Pn());
      if ( !coords_.empty() )
        outfile_->Printf(" %7.2f", Torsion( coords_.XYZ(atom1),
                                            coords_.XYZ(atom2),
                                            coords_.XYZ(atom3),
                                            coords_.XYZ(atom4) ) * Constants::RADDEG );
      outfile_->Printf(" %-*s %-*s %-*s %-*s %*i %*i %*i %*i",
              amn_width_, parm_->AtomMaskName(atom1).c_str(),
              amn_width_, parm_->AtomMaskName(atom2).c_str(),
              amn_width_, parm_->AtomMaskName(atom3).c_str(),
              amn_width_, parm_->AtomMaskName(atom4).c_str(),
              Awidth_, atom1+1, Awidth_, atom2+1,
              Awidth_, atom3+1, Awidth_, atom4+1);
      // Atom types
      outfile_->Printf(" %*s %*s %*s %*s\n",
                       max_type_len_, *((*parm_)[atom1].Type()),
                       max_type_len_, *((*parm_)[atom2].Type()),
                       max_type_len_, *((*parm_)[atom3].Type()),
                       max_type_len_, *((*parm_)[atom4].Type()));
    }
    nd++;
  }
  //outfile_->Printf("\n");
}

// TopInfo::PrintDihedralInfo()
int TopInfo::PrintDihedralInfo(std::string const& mask1exp, std::string const& mask2exp,
                               std::string const& mask3exp, std::string const& mask4exp,
                               bool printImpropers) const
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
  size_t n_torsions;
  const char* label = "Dih";
  bool hasParams;
  if (printImpropers) {
    label = "Imp";
    n_torsions = parm_->Chamber().Impropers().size();
    hasParams = !parm_->Chamber().ImproperParm().empty();
  } else {
    n_torsions = parm_->DihedralsH().size() + parm_->Dihedrals().size();
    hasParams = !parm_->DihedralParm().empty();
  }
  int nw = std::max(3, DigitWidth(n_torsions));
  outfile_->Printf("# %*s", nw, label);
  if (hasParams)
    outfile_->Printf(" %7s %5s %4s", "PK", "Phase", "PN");
  if (!coords_.empty())
    outfile_->Printf(" %7s", "Value");
  outfile_->Printf(" %-*s %-*s %-*s %-*s %*s %*s %*s %*s %*s %*s %*s %*s\n",
                   amn_width_, "Atom1", amn_width_, "Atom2",
                   amn_width_, "Atom3", amn_width_, "Atom4",
                   Awidth_, "A1", Awidth_, "A2",
                   Awidth_, "A3", Awidth_, "A4",
                   max_type_len_, "T1", max_type_len_, "T2",
                   max_type_len_, "T3", max_type_len_, "T4");
  int nd = 1;
  if (printImpropers) {
    PrintDihedrals( parm_->Chamber().Impropers(), parm_->Chamber().ImproperParm(), mask1, mask2, mask3, mask4, nw, nd );
  } else {
    PrintDihedrals( parm_->DihedralsH(), parm_->DihedralParm(), mask1, mask2, mask3, mask4, nw, nd );
    PrintDihedrals( parm_->Dihedrals(),  parm_->DihedralParm(), mask1, mask2, mask3, mask4, nw, nd );
  }
  return 0;
}
