// Action_Pairwise
#include <cmath> //sqrt
#include "Action_Pairwise.h"
#include "CpptrajStdio.h"
#include "Trajout_Single.h"
#include "Constants.h" // ELECTOAMBER
#include "StringRoutines.h" // ByteString()

// CONSTRUCTOR
Action_Pairwise::Action_Pairwise() :
  printMode_(ONLY_CUT),
  nb_calcType_(NORMAL),
  CurrentParm_(0),
  N_ref_interactions_(0),
  nframes_(0),
  ds_vdw_(0),
  ds_elec_(0),
  vdwMat_(0),
  eleMat_(0),
  ELJ_(0),
  Eelec_(0),
  cut_evdw_(1.0),
  cut_eelec_(1.0),
  Eout_(0)
{} 

void Action_Pairwise::Help() const {
  mprintf("\t[<name>] [<mask>] [out <filename>] [cuteelec <ecut>] [cutevdw <vcut>]\n"
          "\t[ %s ] [cutout <cut mol2 prefix>]\n", DataSetList::RefArgs);
  mprintf("\t[vmapout <vdw map>] [emapout <elec map>] [avgout <avg file>]\n"
          "\t[eout <eout file>] [pdbout <pdb file>] [printmode {only|or|and}]\n"
          "  Calculate pairwise (non-bonded) energy for atoms in <mask>.\n"
          "  If 'eout' is specified individual interaction energies will be written to\n"
          "  <eout file>. If a reference structure is given the energies will be\n"
          "  Eref - Eframe. Only energies with absolute value greater than <ecut> and\n"
          "  <vcut> (by default 1.0 kcal/mol) will be printed.\n"
          "  printmode only : Only print energy cutoff is satisfied.\n"
          "            or   : Print both energies if either cutoff is satisfied.\n"
          "            and  : Print both energies if both cutoffs are satisfied.\n");
}

const double Action_Pairwise::QFAC = Constants::ELECTOAMBER * Constants::ELECTOAMBER;

// Action_Pairwise::Init()
Action::RetType Action_Pairwise::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
# ifdef MPI
  if (init.TrajComm().Size() > 1) {
    mprinterr("Error: 'pairwise' action does not work with > 1 thread (%i threads currently).\n",
              init.TrajComm().Size());
    return Action::ERR;
  }
# endif
  // Get Keywords
  DataFile* dataout = init.DFL().AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  DataFile* vmapout = init.DFL().AddDataFile( actionArgs.GetStringKey("vmapout"), actionArgs );
  DataFile* emapout = init.DFL().AddDataFile( actionArgs.GetStringKey("emapout"), actionArgs );
  avgout_ = actionArgs.GetStringKey("avgout");
  std::string eout = actionArgs.GetStringKey("eout");
  cut_eelec_ = fabs(actionArgs.getKeyDouble("cuteelec",1.0));
  cut_evdw_ = fabs(actionArgs.getKeyDouble("cutevdw",1.0));
  mol2Prefix_ = actionArgs.GetStringKey("cutout");
  std::string pdbout = actionArgs.GetStringKey("pdbout");
  printMode_ = ONLY_CUT;
  std::string pmode = actionArgs.GetStringKey("printmode");
  if (!pmode.empty()) {
    if (pmode == "only") printMode_ = ONLY_CUT;
    else if (pmode == "or") printMode_ = OR_CUT;
    else if (pmode == "and") printMode_ = AND_CUT;
    else {
      mprinterr("Error: Unrecognized print mode: %s\n", pmode.c_str());
      return Action::ERR;
    }
  }
  ReferenceFrame REF = init.DSL().GetReferenceFrame( actionArgs );
  
  // Get Masks
  Mask0_.SetMaskString( actionArgs.GetMaskNext() );
  std::string refmask = actionArgs.GetMaskNext();
  if (!refmask.empty())
    RefMask_.SetMaskString(refmask);
  else
    RefMask_.SetMaskString( Mask0_.MaskString() );

  // Datasets
  std::string ds_name = actionArgs.GetStringNext();
  if (ds_name.empty())
    ds_name = init.DSL().GenerateDefaultName("PW");
  ds_vdw_  = init.DSL().AddSet(DataSet::DOUBLE, MetaData(ds_name, "EVDW"));
  ds_elec_ = init.DSL().AddSet(DataSet::DOUBLE, MetaData(ds_name, "EELEC"));
  if (ds_vdw_ == 0 || ds_elec_ == 0) return Action::ERR;
  // Add DataSets to data file list
  if (dataout != 0) {
    dataout->AddDataSet(ds_vdw_);
    dataout->AddDataSet(ds_elec_);
  }
  vdwMat_ = (DataSet_MatrixDbl*)init.DSL().AddSet(DataSet::MATRIX_DBL, MetaData(ds_name, "VMAP"));
  eleMat_ = (DataSet_MatrixDbl*)init.DSL().AddSet(DataSet::MATRIX_DBL, MetaData(ds_name, "EMAP"));
  if (vdwMat_ == 0 || eleMat_ == 0) return Action::ERR;
  if (vmapout != 0) vmapout->AddDataSet(vdwMat_);
  if (emapout != 0) emapout->AddDataSet(eleMat_);

  // Get reference structure
  if (REF.error()) return Action::ERR;
  if (!REF.empty()) { 
    // Set up reference mask
    if ( REF.Parm().SetupIntegerMask(RefMask_) ) return Action::ERR;
    if (RefMask_.None()) {
      mprinterr("Error: No atoms selected in reference mask.\n");
      return Action::ERR;
    }
    // Set up nonbonded params for reference
    if ( (N_ref_interactions_=SetupNonbondParm( RefMask_, REF.Parm() )) == -1 ) 
      return Action::ERR;
    // Calculate energy for reference
    nb_calcType_ = SET_REF;
    NonbondEnergy(REF.Coord(), REF.Parm(), RefMask_);
    nb_calcType_ = COMPARE_REF;
  }

  // Output for individual atom energy | dEnergy
  if (!eout.empty()) {
    Eout_ = init.DFL().AddCpptrajFile(eout, "Atom Energies");
    if (Eout_ == 0) {
      mprinterr("Error: Could not set up file %s for eout.\n",eout.c_str());
      return Action::ERR;
    }
  }

  // Set up output pdb
  if (!pdbout.empty()) {
    if (PdbOut_.OpenWrite( pdbout )) return Action::ERR;
  }

  // Action Info
  mprintf("    PAIRWISE: Atoms in mask [%s].\n",Mask0_.MaskString());
  if (!eout.empty())
    mprintf("\tEnergy info for each atom will be written to %s\n",eout.c_str());
  if (nb_calcType_ == COMPARE_REF) { 
    mprintf("\tReference %s, mask [%s]\n", REF.refName(), RefMask_.MaskString());
    mprintf("\tReference energy (kcal/mol): EVDW= %12.5e  EELEC= %12.5e\n", ELJ_, Eelec_);
    mprintf("\tSize of reference energy array is %zu elements (%s)\n",
            ref_nonbondEnergy_.size(),
            ByteString(ref_nonbondEnergy_.size() * 2 * sizeof(double), BYTE_DECIMAL).c_str());
  }
  mprintf("\tEelec print absolute cutoff (kcal/mol): %.4f\n", cut_eelec_);
  mprintf("\tEvdw print absolute cutoff (kcal/mol) : %.4f\n", cut_evdw_);
  if (!mol2Prefix_.empty())
    mprintf("\tAtoms satisfying cutoff will be printed to %s.e<type>.mol2\n",
            mol2Prefix_.c_str());
  if (PdbOut_.IsOpen())
    mprintf("\tPDB with evdw/eelec in occ/b-fac columns will be written to %s\n",
            PdbOut_.Filename().full());
  
  return Action::OK;
}

// Action_Pairwise::SetupNonbondParm()
/** Set up the exclusion list based on the given mask and parm.
  * \return the total number of interactions, -1 on error.
  */
int Action_Pairwise::SetupNonbondParm(AtomMask const& maskIn, Topology const& ParmIn)
{
  // Check if LJ parameters present - need at least 2 atoms for it to matter.
  if (ParmIn.Natom() > 1 && !ParmIn.Nonbond().HasNonbond()) {
    mprinterr("Error: Topology does not have LJ information.\n");
    return -1;
  }

  // Determine the actual number of pairwise interactions that will be calcd.
  unsigned int n_interactions = 0;
  for (AtomMask::const_iterator at0 = maskIn.begin(); at0 != maskIn.end(); ++at0) {
    Atom::excluded_iterator ex = ParmIn[*at0].excludedbegin();
    for (AtomMask::const_iterator at1 = at0 + 1; at1 != maskIn.end(); ++at1) {
      if (ex != ParmIn[*at0].excludedend() && *at1 == *ex)
        // Atom 1 is excluded from Atom0; just increment to next excluded atom.
        ++ex;
      else
        ++n_interactions;
    }
  }

  // Print total number of interactions for this parm
  mprintf("\t%u interactions for topology '%s'.\n", n_interactions, ParmIn.c_str());

  // DEBUG - Print exclusion list for each atom
  /*for (unsigned int atom = 0; atom < exclusionList.size(); atom++) {
    mprintf("\t%8u:",atom + 1);
    for (std::vector<int>::iterator eat = exclusionList[atom].begin();
                                    eat != exclusionList[atom].end();
                                    eat++)
    {
      mprintf(" %i",*eat + 1);
    }
    mprintf("\n");
  }*/
  return (int)n_interactions;
}

// Action_Pairwise::Setup()
/** Set up mask, allocate memory for exclusion list.
  */
Action::RetType Action_Pairwise::Setup(ActionSetup& setup) {
  // Set up mask
  if ( setup.Top().SetupIntegerMask( Mask0_ ) ) return Action::ERR;
  if (Mask0_.None()) {
    mprintf("Warning: Mask has no atoms.\n");
    return Action::SKIP;
  }

  // Set up exclusion list and determine total # interactions.
  int N_interactions = SetupNonbondParm(Mask0_, setup.Top());
  if (N_interactions < 0) return Action::ERR;
  if (N_interactions < 1) {
    mprintf("Warning: No pairwise interactions to calculate for mask '%s'\n", Mask0_.MaskString());
    return Action::SKIP;
  }

  // Allocate/check matrix memory
  if (vdwMat_->Size() == 0) {
    vdwMat_->AllocateTriangle( Mask0_.Nselected() );
    eleMat_->AllocateTriangle( Mask0_.Nselected() );
  } else {
    size_t nselected = (size_t)Mask0_.Nselected();
    size_t newSize = (nselected * (nselected-1)) / 2;
    if (vdwMat_->Size() != newSize) {
      mprinterr("Error: Attempting to reallocate matrix with different size.\n"
                "Error:   Original size= %zu, new size= %zu\n"
                "Error:   This can occur when different #s of atoms are selected in\n"
                "Error:   different topology files.\n", vdwMat_->Size(), newSize);
      return Action::ERR;
    }
  }

  // If comparing to a reference frame for atom-by-atom comparison make sure
  // the number of interactions is the same in reference and parm.
  if (nb_calcType_ == COMPARE_REF) {
    if (N_interactions != N_ref_interactions_) {
      mprinterr("Error: # reference interactions (%i) != # interactions for this parm (%i)\n",
                N_ref_interactions_, N_interactions);
      return Action::ERR;
    }
  }
  // Set up cumulative energy arrays
  atom_eelec_.clear();
  atom_eelec_.resize(Mask0_.Nselected(), 0.0);
  atom_evdw_.clear();
  atom_evdw_.resize(Mask0_.Nselected(), 0.0);
  // Print pairwise info for this parm
  Mask0_.MaskInfo();
  CurrentParm_ = setup.TopAddress();
  return Action::OK;  
}

// Action_Pairwise::WriteEnergies()
void Action_Pairwise::WriteEnergies(Topology const& parmIn, int atom1, int atom2,
                                    double evdw, double eelec, const char* etype)
{
  if (fabs(evdw) > cut_evdw_) {
    Eout_->Printf("\tAtom %6i@%4s-%6i@%4s %sEvdw= %12.4f\n",
                    atom1+1, parmIn[atom1].c_str(),
                    atom2+1, parmIn[atom2].c_str(),
                    etype,  evdw);
  }
  if (fabs(eelec) > cut_eelec_) {
    Eout_->Printf("\tAtom %6i@%4s-%6i@%4s %sEelec= %12.4f\n",
                    atom1+1, parmIn[atom1].c_str(),
                    atom2+1, parmIn[atom2].c_str(),
                    etype, eelec);
  }
}

// Action_Pairwise::NonbondEnergy()
/** Calculate non-bonded energy using the nonbondParm array. The total
  * LJ (vdw) energy is put in ELJ, and the total Coulomb (elec) energy
  * is put in Eelec. Depending on the value of nb_calcType, each pair
  * energy is either compared to a reference, distributed over both atoms
  * evenly in the cumulative array, or reference values are set. If comparing
  * to a reference structure, pairs for which the energy difference exceeds
  * the cutoffs are printed.
  */
void Action_Pairwise::NonbondEnergy(Frame const& frameIn, Topology const& parmIn, 
                                    AtomMask const& maskIn)
{
  double delta2;
  NonbondEnergyType refE;

  ELJ_ = 0.0;
  Eelec_ = 0.0;
  std::vector<NonbondEnergyType>::const_iterator refpair = ref_nonbondEnergy_.begin();
  // Loop over all atom pairs and set information
  // Outer loop
  for (int idx1 = 0; idx1 != maskIn.Nselected(); idx1++)
  {
    int maskatom1 = maskIn[idx1];
    // Get coordinates for first atom.
    Vec3 coord1 = frameIn.XYZ( maskatom1 );
    // Set up exclusion list for this atom
    Atom::excluded_iterator excluded_atom = parmIn[maskatom1].excludedbegin();
    // Inner loop
    for (int idx2 = idx1 + 1; idx2 != maskIn.Nselected(); idx2++)
    {
      int maskatom2 = maskIn[idx2];
      // If atom is excluded, just increment to next excluded atom;
      // otherwise perform energy calc.
      if ( excluded_atom != parmIn[maskatom1].excludedend() && maskatom2 == *excluded_atom )
        ++excluded_atom;
      else {
        // Calculate the vector pointing from atom2 to atom1
        Vec3 JI = coord1 - Vec3(frameIn.XYZ( maskatom2 ));
        double rij2 = JI.Magnitude2();
        // Normalize
        double rij = sqrt(rij2);
        JI /= rij;
        // LJ energy
        NonbondType const& LJ = parmIn.GetLJparam(maskatom1, maskatom2);
        double r2    = 1.0 / rij2;
        double r6    = r2 * r2 * r2;
        double r12   = r6 * r6;
        double f12   = LJ.A() * r12;  // A/r^12
        double f6    = LJ.B() * r6;   // B/r^6
        double e_vdw = f12 - f6;     // (A/r^12)-(B/r^6)
        ELJ_ += e_vdw;
        // LJ Force 
        //force=((12*f12)-(6*f6))*r2; // (12A/r^13)-(6B/r^7)
        //scalarmult(f,JI,F);
        // Coulomb energy 
        double qiqj = QFAC * parmIn[maskatom1].Charge() * parmIn[maskatom2].Charge();
        double e_elec = qiqj / rij;
        Eelec_ += e_elec;
        // Coulomb Force
        //force=e_elec/rij; // kes_*(qiqj/r)*(1/r)
        //scalarmult(f,JI,F);

        // ----------------------------------------
        if (nb_calcType_ == COMPARE_REF) {
          // 1 - Comparison to reference, cumulative dEnergy on atoms
          // dEvdw
          double delta_vdw = refpair->evdw - e_vdw;
          // dEelec
          double delta_eelec = refpair->eelec - e_elec;
          // Output
          if (Eout_ != 0)
            WriteEnergies(parmIn, maskatom1, maskatom2, delta_vdw, delta_eelec, "d");
          vdwMat_->Element(idx1, idx2) += delta_vdw;
          eleMat_->Element(idx1, idx2) += delta_eelec;
          // Divide the total pair dEvdw between both atoms.
          delta2 = delta_vdw * 0.5;
          atom_evdw_[idx1] += delta2;
          atom_evdw_[idx2] += delta2;
          // Divide the total pair dEelec between both atoms.
          delta2 = delta_eelec * 0.5;
          atom_eelec_[idx1] += delta2;
          atom_eelec_[idx2] += delta2;
        } else if (nb_calcType_ == NORMAL) {
          // 2 - No reference, just cumulative Energy on atoms
          if (Eout_ != 0)
            WriteEnergies(parmIn, maskatom1, maskatom2, e_vdw, e_elec, "");
          vdwMat_->Element(idx1, idx2) += e_vdw;
          eleMat_->Element(idx1, idx2) += e_elec;
          // Cumulative evdw - divide between both atoms
          delta2 = e_vdw * 0.5;
          atom_evdw_[idx1] += delta2;
          atom_evdw_[idx2] += delta2;
          // Cumulative eelec - divide between both atoms
          delta2 = e_elec * 0.5;
          atom_eelec_[idx1] += delta2;
          atom_eelec_[idx2] += delta2;
        } else { // if nb_calcType_ == SET_REF
          // 3 - Store the reference nonbond energy for this pair
          refE.evdw = e_vdw;
          refE.eelec = e_elec;
          ref_nonbondEnergy_.push_back( refE );
        }
        ++refpair;
        // ----------------------------------------
      } // END pair not excluded
    } // END Inner loop
  } // END Outer loop
}

// Action_Pairwise::WriteCutFrame()
/** Write file containing only cut atoms and energies as charges. */
int Action_Pairwise::WriteCutFrame(int frameNum, Topology const& Parm, AtomMask const& CutMask, 
                                   Darray const& CutCharges,
                                   Frame const& frame, std::string const& outfilename) 
{
  if (CutMask.Nselected() != (int)CutCharges.size()) {
    mprinterr("Error: WriteCutFrame: # of charges (%u) != # mask atoms (%i)\n",
              CutCharges.size(), CutMask.Nselected());
    return 1;
  }
  Frame CutFrame(frame, CutMask);
  Topology* CutParm = Parm.modifyStateByMask( CutMask );
  if (CutParm == 0) return 1;
  // Set new charges
  for (int i = 0; i != CutParm->Natom(); i++)
    CutParm->SetAtom(i).SetCharge( CutCharges[i] );
  int err = 0;
  Trajout_Single tout;
  if (tout.PrepareTrajWrite(outfilename, "multi", CutParm, CoordinateInfo(), 1,
                            TrajectoryFile::MOL2FILE))
  {
    mprinterr("Error: Could not set up cut mol2 file %s\n", outfilename.c_str());
    err = 1;
  } else {
    tout.WriteSingle(frameNum, CutFrame);
    tout.EndTraj();
  }
  delete CutParm;
  return err;
}

static const char* CalcString[] = { "Evdw", "Eelec" };
static const char* CutName[] = { ".evdw.mol2", ".eelec.mol2" };

// Action_Pairwise::PrintCutAtoms()
/** Print atoms for which the cumulative energy satisfies the given
  * cutoffs. Also create MOL2 files containing those atoms.
  */
int Action_Pairwise::PrintCutAtoms(Frame const& frame, int frameNum, EoutType ctype,
                                   Darray const& Earray, double cutIn)
{
  AtomMask CutMask;  // Hold atoms that satisfy the cutoff
  Darray CutCharges; // Hold evdw/eelec corresponding to CutMask atoms.

  if (Eout_ != 0) {
    if (nb_calcType_==COMPARE_REF)
      Eout_->Printf("\tPAIRWISE: Cumulative d%s:", CalcString[ctype]);
    else
      Eout_->Printf("\tPAIRWISE: Cumulative %s:", CalcString[ctype]);
    Eout_->Printf(" %s < %.4f, %s > %.4f\n", CalcString[ctype], -cutIn,
                 CalcString[ctype], cutIn);
  }
  for (int idx = 0; idx != Mask0_.Nselected(); idx++)
  {
    int atom = Mask0_[idx];
    if (fabs(Earray[idx]) > cutIn)
    {
      if (Eout_ != 0) 
        Eout_->Printf("\t\t%6i@%s: %12.4f\n", atom+1,
                    (*CurrentParm_)[atom].c_str(), Earray[idx]);
      CutMask.AddAtom(atom);
      CutCharges.push_back(Earray[idx]);
    }
  }
  // Write mol2 with atoms satisfying cutoff
  if (!mol2Prefix_.empty() && !CutMask.None()) {
    if (WriteCutFrame(frameNum, *CurrentParm_, CutMask, CutCharges, 
                      frame, mol2Prefix_ + CutName[ctype])) 
      return 1;
  }

  return 0;
}

// Action_Pairwise::DoAction()
Action::RetType Action_Pairwise::DoAction(int frameNum, ActionFrame& frm) {
  // Reset cumulative energy arrays
  atom_eelec_.assign(Mask0_.Nselected(), 0.0);
  atom_evdw_.assign(Mask0_.Nselected(), 0.0);
  if (Eout_ != 0) Eout_->Printf("PAIRWISE: Frame %i\n",frm.TrajoutNum());
  NonbondEnergy( frm.Frm(), *CurrentParm_, Mask0_ );
  nframes_++;
  // Write cumulative energy arrays
  if (PrintCutAtoms( frm.Frm(), frm.TrajoutNum(), VDWOUT, atom_evdw_, cut_evdw_ ))
    return Action::ERR;
  if (PrintCutAtoms( frm.Frm(), frm.TrajoutNum(), ELECOUT, atom_eelec_, cut_eelec_ ))
    return Action::ERR;
  // Write PDB with atoms that satisfy cutoff colored in.
  if (PdbOut_.IsOpen()) {
    PdbOut_.WriteMODEL(frm.TrajoutNum() + 1); // FIXME in parallel this needs to be separate files
    for (int idx = 0; idx != Mask0_.Nselected(); idx++)
    {
      int atom = Mask0_[idx];
      float occ = 0.0;
      float bfac = 0.0;
      if (fabs(atom_evdw_[idx]) > cut_evdw_)
        occ = (float)atom_evdw_[idx];
      if (fabs(atom_eelec_[idx]) > cut_eelec_)
        bfac = (float)atom_eelec_[idx];
      const double* XYZ = frm.Frm().XYZ( atom );
      Atom const& AT = (*CurrentParm_)[atom];
      int rn = AT.ResNum();
      PdbOut_.WriteCoord(PDBfile::ATOM, atom+1, AT.c_str(), CurrentParm_->Res(rn).c_str(),
                         rn + 1, XYZ[0], XYZ[1], XYZ[2], occ, bfac, AT.ElementName(),
                         (int)AT.Charge());
    }
    PdbOut_.WriteENDMDL();
  }
  ds_vdw_->Add(frameNum, &ELJ_);
  ds_elec_->Add(frameNum, &Eelec_);

  return Action::OK;
}

// Action_Pairwise::Print()
void Action_Pairwise::Print() {
  if (nframes_ < 1) return;
  // Divide matrices by # of frames
  double norm = 1.0 / (double)nframes_;
  for (unsigned int i = 0; i != vdwMat_->Size(); i++)
  {
    (*vdwMat_)[i] *= norm;
    (*eleMat_)[i] *= norm;
  }
  // Write out final results
  CpptrajFile AvgOut;
  if (AvgOut.OpenWrite( avgout_ )) return;
  if (nb_calcType_ == NORMAL)
    mprintf("  PAIRWISE: Writing all pairs with |<evdw>| > %.4f, |<eelec>| > %.4f\n",
            cut_evdw_, cut_eelec_);
  else if (nb_calcType_ == COMPARE_REF)
    mprintf("  PAIRWISE: Writing all pairs with |<dEvdw>| > %.4f, |<dEelec>| > %.4f\n",
            cut_evdw_, cut_eelec_);
  AvgOut.Printf("%-16s %5s -- %16s %5s : ENE\n","#Name1", "At1", "Name2", "At2");
  for (int idx1 = 0; idx1 != Mask0_.Nselected(); idx1++)
  {
    int m1 = Mask0_[idx1];
    for (int idx2 = idx1 + 1; idx2 != Mask0_.Nselected(); idx2++)
    {
      int m2 = Mask0_[idx2];
      double EV = vdwMat_->GetElement(idx1, idx2);
      double EE = eleMat_->GetElement(idx1, idx2);
      bool outputv = ( fabs(EV) > cut_evdw_ );
      bool outpute = ( fabs(EE) > cut_eelec_ );
      if (printMode_ == ONLY_CUT) {
        if (outputv || outpute) {
          AvgOut.Printf("%16s %5i -- %16s %5i :",
                  CurrentParm_->TruncResAtomName(m1).c_str(), m1 + 1,
                  CurrentParm_->TruncResAtomName(m2).c_str(), m2 + 1);
          if (outputv) AvgOut.Printf("  EVDW= %12.5e", EV);
          if (outpute) AvgOut.Printf(" EELEC= %12.5e", EE);
          AvgOut.Printf("\n");
        }
      } else if (printMode_ == OR_CUT) {
        if (outputv || outpute)
          AvgOut.Printf("%16s %5i -- %16s %5i :  EVDW= %12.5e EELEC= %12.5e\n",
                        CurrentParm_->TruncResAtomName(m1).c_str(), m1 + 1,
                        CurrentParm_->TruncResAtomName(m2).c_str(), m2 + 1,
                        EV, EE);
      } else if (printMode_ == AND_CUT) {
        if (outputv && outpute)
          AvgOut.Printf("%16s %5i -- %16s %5i :  EVDW= %12.5e EELEC= %12.5e\n",
                        CurrentParm_->TruncResAtomName(m1).c_str(), m1 + 1,
                        CurrentParm_->TruncResAtomName(m2).c_str(), m2 + 1,
                        EV, EE);
      }
    }
  }
}
