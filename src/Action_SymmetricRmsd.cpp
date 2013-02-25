#include "Hungarian.h" 
#include "Action_SymmetricRmsd.h"
#include "CpptrajStdio.h"
#include "AtomMap.h"
#include "DistRoutines.h"

// CONSTRUCTOR
Action_SymmetricRmsd::Action_SymmetricRmsd() {}

void Action_SymmetricRmsd::Help() {
  mprintf("\t[<name>] <mask> [<refmask>] [out filename] [nofit | norotate] [mass]\n");
  mprintf("\t[ first | ref <filename> | refindex <#> |\n");
  mprintf("\treftraj <filename> [parm <parmname> | parmindex <#>] ]\n");
}

// Action_SymmetricRmsd::Init()
Action::RetType Action_SymmetricRmsd::Init(ArgList& actionArgs, TopologyList* PFL, 
                          FrameList* FL, DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  // Check for keywords
  GetRmsKeywords( actionArgs );
  DataFile* outfile = DFL->AddDataFile(actionArgs.GetStringKey("out"), actionArgs);
  // Reference keywords
  bool previous = actionArgs.hasKey("previous");
  bool first = actionArgs.hasKey("first");
  ReferenceFrame  REF = FL->GetFrame( actionArgs );
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
  // FIXME: No mass information yet
  rmsTgtFrame_.SetupFrame( TgtMask().Nselected() );
  // Reference setup
  if (SetupRef(*currentParm, TgtMask().Nselected(), "symmrmsd"))
    return Action::ERR;
  // Check for symmetric atoms
  //AtomMask cMask = TgtMask();
  //cMask.ConvertMaskType(); // Convert to char mask
  // Create initial 1 to 1 atom map
  AMap_.clear();
  for (int atom = 0; atom < currentParm->Natom(); atom++)
    AMap_.push_back(atom);
  // Create atom maps for each residue
  SymmetricAtomIndices_.clear();
  AtomMap resmap;
  resmap.SetDebug(0); // DEBUG
  for (int residue = 0; residue < currentParm->Nres(); ++residue) {
    int res_first_atom = currentParm->Res(residue).FirstAtom();
    mprintf("DEBUG: Residue %s\n", currentParm->TruncResNameNum(residue).c_str());
    if (resmap.SetupResidue(*currentParm, residue) != 0) return Action::ERR;
    if (resmap.CheckBonds() != 0) return Action::ERR;
    resmap.DetermineAtomIDs();
    // NOTE: Indices for resmap start at 0.
    // Generate maps for symmetric atoms
    std::vector<bool> selected(resmap.Natom(), false);
    for (int atom1 = 0; atom1 < resmap.Natom(); atom1++) {
      if (!selected[atom1]) {
        if (resmap[atom1].Nduplicated() > 0) { // This atom is duplicated
          selected[atom1] = true;
          Iarray symmatoms(1, atom1 + res_first_atom);
          for (int atom2 = atom1 + 1; atom2 < resmap.Natom(); atom2++) {
            if (resmap[atom1].Unique() == resmap[atom2].Unique()) {
              selected[atom2] = true;
              symmatoms.push_back(atom2 + res_first_atom);
            }
          } // End loop over atom2
          SymmetricAtomIndices_.push_back( symmatoms );
          mprintf("DEBUG:\t\tAtom %s ID %s is duplicated %u times:", 
                  currentParm->TruncResAtomName(symmatoms.front()).c_str(), 
                  resmap[atom1].Unique().c_str(), symmatoms.size());
          for (Iarray::const_iterator sa = symmatoms.begin(); sa != symmatoms.end(); ++sa)
            mprintf(" %i", *sa + 1);
          mprintf("\n");
        }
      }
    } // End loop over atom1
    mprintf("DEBUG:\tNon-symmetric atoms:");
    for (int atom1 = 0; atom1 < resmap.Natom(); atom1++)
      if (!selected[atom1]) mprintf(" %i", atom1 + res_first_atom + 1);
    mprintf("\n");
  } // End loop over residue

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
        mprintf("\t\t%i to %i: %f\n", *tgtatom + 1, *refatom + 1, dist2);
        cost_matrix.AddElement( dist2 ); 
      }
    }
    mprintf("\tSymmetric atoms starting with %i", (*symmatoms).front() + 1);
    cost_matrix.Print("Cost Matrix:"); // DEBUG
    Hungarian HA(cost_matrix);
    Iarray resMap = HA.Optimize();
    // Fill in overall map
    Iarray::iterator rmap = resMap.begin();
    for (Iarray::iterator atmidx = (*symmatoms).begin();
                          atmidx != (*symmatoms).end(); ++atmidx, ++rmap)
    {
      AMap_[*atmidx] = (*symmatoms)[*rmap];
      mprintf("\tAssigned atom %i to atom %i\n", *atmidx + 1, (*symmatoms)[*rmap] + 1);
    }
  }
/*  int ref = 0;
  for (Iarray::iterator map = AMap_.begin(); map != AMap_.end(); ++map)
    mprintf("\t%i -> %i\n", ref++, *map);*/
  rmsTgtFrame_.SetTargetByMap(*currentFrame, AMap_);
  double rmsdval = CalcRmsd(rmsTgtFrame_, SelectedRef(), RefTrans());
  rmsd_->Add(frameNum, &rmsdval);
  return Action::OK;
}
