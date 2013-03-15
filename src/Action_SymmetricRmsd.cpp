#include "Hungarian.h" 
#include "Action_SymmetricRmsd.h"
#include "CpptrajStdio.h"
#include "AtomMap.h"
#include "DistRoutines.h"

// CONSTRUCTOR
Action_SymmetricRmsd::Action_SymmetricRmsd() : debug_(0) {}

void Action_SymmetricRmsd::Help() {
  mprintf("\t[<name>] <mask> [<refmask>] [out <filename>] [nofit | norotate] [mass]\n");
  mprintf("\t[ first | ref <filename> | refindex <#> |\n");
  mprintf("\t  reftraj <trajname> [parm <parmname> | parmindex <#>] ]\n");
  mprintf("\tPerform symmetry-corrected RMSD calculation.\n");
}

// Action_SymmetricRmsd::Init()
Action::RetType Action_SymmetricRmsd::Init(ArgList& actionArgs, TopologyList* PFL, 
                          FrameList* FL, DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  debug_ = debugIn;
  // Check for keywords
  GetRmsKeywords( actionArgs );
  DataFile* outfile = DFL->AddDataFile(actionArgs.GetStringKey("out"), actionArgs);
  // Reference keywords
  bool previous = actionArgs.hasKey("previous");
  bool first = actionArgs.hasKey("first");
  ReferenceFrame  REF = FL->GetFrameFromArgs( actionArgs );
  std::string reftrajname = actionArgs.GetStringKey("reftraj");
  Topology* RefParm = PFL->GetParm( actionArgs );
  // Get the RMS mask string for target
  std::string mask1 = GetRmsMasks(actionArgs);
  // Initialize reference
  if (InitRef(previous, first, UseMass(), Fit(), reftrajname, REF, RefParm, mask1,
              actionArgs, "symmrmsd"))
    return Action::ERR;
  // Set up the RMSD data set. 
  rmsd_ = DSL->AddSet(DataSet::DOUBLE, actionArgs.GetStringNext(),"RMSD");
  if (rmsd_==0) return Action::ERR;
  rmsd_->SetScalar( DataSet::M_RMS );
  // Add dataset to data file list
  if (outfile != 0) outfile->AddSet( rmsd_ );
  
  mprintf("    SYMMRMSD: (%s), reference is %s", TgtMask().MaskString(),
          RefModeString());
  PrintRmsStatus();
  return Action::OK;
}

// Action_SymmetricRmsd::Setup()
Action::RetType Action_SymmetricRmsd::Setup(Topology* currentParm, Topology** parmAddress) {
  // Target setup
  if (SetupRmsMask(*currentParm, "symmrmsd")) return Action::ERR;
  // Target remap setup
  remapFrame_.SetupFrameM( currentParm->Atoms() );
  // Reference setup
  if (SetupRef(*currentParm, TgtMask().Nselected(), "symmrmsd"))
    return Action::ERR;
  //InitialFitMask.ResetMask();
  // Create char mask to see what symmetric atoms are selected. 
  AtomMask cMask = TgtMask();
  cMask.ConvertToCharMask();
  // Create initial 1 to 1 atom map for all atoms; indices in 
  // SymmetricAtomIndices will correspond to positions in AMap.
  AMap_.clear();
  for (int atom = 0; atom < currentParm->Natom(); atom++)
    AMap_.push_back(atom);
  // Determine last selected residue
  int last_res = (*currentParm)[TgtMask().back()].ResNum() + 1;
  // Create atom maps for each selected atom in residues
  SymmetricAtomIndices_.clear();
  AtomMap resmap;
  if (debug_ > 1) resmap.SetDebug(1);
  for (int residue = 0; residue < last_res; ++residue) {
    int res_first_atom = currentParm->Res(residue).FirstAtom();
    if (debug_ > 0)
      mprintf("DEBUG: Residue %s\n", currentParm->TruncResNameNum(residue).c_str());
    if (resmap.SetupResidue(*currentParm, residue) != 0) return Action::ERR;
    if (resmap.CheckBonds() != 0) return Action::ERR;
    resmap.DetermineAtomIDs();
    // NOTE: Indices for resmap start at 0.
    // Generate maps for symmetric atoms
    // AtomStatus: 0=Unselected, 1=Selected/Non-symm., 2=Selected/Symm
    std::vector<int> AtomStatus(resmap.Natom(), 0);
    for (int atom1 = 0; atom1 < resmap.Natom(); atom1++) {
      int actual_atom1 = atom1 + res_first_atom; // Actual atom index in currentParm
      if (cMask.AtomInCharMask(actual_atom1) && AtomStatus[atom1] == 0) {
        AtomStatus[atom1] = 1; // Initially select as non-symmetric
        // Check if atom is duplicated and not bound to a chiral center. 
        // If so, find all selected duplicates.
        if (!resmap[atom1].BoundToChiral() && resmap[atom1].Nduplicated() > 0) {
          AtomStatus[atom1] = 2; // Select as symmetric
          Iarray symmatoms(1, actual_atom1);
          for (int atom2 = atom1 + 1; atom2 < resmap.Natom(); atom2++) {
            int actual_atom2 = atom2 + res_first_atom;
            if (cMask.AtomInCharMask(actual_atom2) &&
                resmap[atom1].Unique() == resmap[atom2].Unique() &&
                !resmap[atom2].BoundToChiral()) 
            {
              AtomStatus[atom2] = 2; // Select as symmetric
              symmatoms.push_back(actual_atom2);
            }
          } // END loop over atom2
          if (debug_ > 0)
            mprintf("DEBUG:\t\tAtom %s ID %s is duplicated %u times:", 
                    currentParm->TruncResAtomName(symmatoms.front()).c_str(), 
                    resmap[atom1].Unique().c_str(), symmatoms.size());
          if (symmatoms.size() > 1) {
            SymmetricAtomIndices_.push_back( symmatoms );
            if (debug_ > 0)
              for (Iarray::const_iterator sa = symmatoms.begin(); sa != symmatoms.end(); ++sa)
                mprintf(" %i", *sa + 1);
          } else {
            // Only one atom selected, no symmetry. Change to non-symmetric.
            AtomStatus[symmatoms.front() - res_first_atom] = 1; // Select as non-symmetric
          }
          if (debug_ > 0) mprintf("\n");
        } // END if atom is duplicated
      }
    } // END loop over atom1
    // TODO: If fitting, set up mask to perform initial fit with selected nonsymmetric atoms
    if (debug_ > 0 && Fit()) {
      mprintf("DEBUG:\tSelected Non-symmetric atoms:");
      for (int atom1 = 0; atom1 < resmap.Natom(); atom1++)
        if (AtomStatus[atom1] == 1) { // If selected/non-symmetric
          mprintf(" %i", atom1 + res_first_atom + 1);
          //InitialFitMask.AddAtom(atom1 + res_first_atom);
        }
      mprintf("\n");
    }
  } // End loop over residue
  //if (Fit()) 
  //  mprintf("DEBUG: Initial fit mask has %i atoms.\n", InitialFitMask.Nselected());

  return Action::OK;
}

// Action_SymmetricRmsd::DoAction()
Action::RetType Action_SymmetricRmsd::DoAction(int frameNum, Frame* currentFrame, 
                                               Frame** frameAddress) 
{
  Matrix_2D cost_matrix;
  // Perform any needed reference actions
  ActionRef( *currentFrame, Fit(), UseMass() );
  // Calculate initial best-fit RMSD
  if (Fit())
    CalcRmsd( *currentFrame, SelectedRef(), RefTrans() );
  // Correct RMSD for symmetry
  // FIXME: currently calcs all atom RMSD
  for (AtomIndexArray::iterator symmatoms = SymmetricAtomIndices_.begin();
                                symmatoms != SymmetricAtomIndices_.end(); ++symmatoms)
  {
    // For each array of symmetric atoms, determine the lowest distance score
    cost_matrix.Setup((*symmatoms).size(), (*symmatoms).size());
    for (Iarray::iterator tgtatom = (*symmatoms).begin(); 
                          tgtatom != (*symmatoms).end(); ++tgtatom)
    {
      for (Iarray::iterator refatom = (*symmatoms).begin();
                            refatom != (*symmatoms).end(); ++refatom)
      {
        double dist2 = DIST2_NoImage( RefFrame().XYZ(*refatom), currentFrame->XYZ(*tgtatom) );
//        mprintf("\t\t%i to %i: %f\n", *tgtatom + 1, *refatom + 1, dist2);
        cost_matrix.AddElement( dist2 ); 
      }
    }
//    mprintf("\tSymmetric atoms starting with %i", (*symmatoms).front() + 1);
    Hungarian HA(cost_matrix);
    Iarray resMap = HA.Optimize();
    // Fill in overall map
    Iarray::iterator rmap = resMap.begin();
    for (Iarray::iterator atmidx = (*symmatoms).begin();
                          atmidx != (*symmatoms).end(); ++atmidx, ++rmap)
    {
      AMap_[*atmidx] = (*symmatoms)[*rmap];
//      mprintf("\tAssigned atom %i to atom %i\n", *atmidx + 1, (*symmatoms)[*rmap] + 1);
    }
  }
/*  int ref = 0;
  for (Iarray::iterator map = AMap_.begin(); map != AMap_.end(); ++map)
    mprintf("\t%i -> %i\n", ref++, *map);*/
  remapFrame_.SetCoordinatesByMap(*currentFrame, AMap_);
  double rmsdval = CalcRmsd(remapFrame_, SelectedRef(), RefTrans());
  rmsd_->Add(frameNum, &rmsdval);
  return Action::OK;
}
