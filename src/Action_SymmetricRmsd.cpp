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
  AtomMap resmap;
  resmap.SetDebug(0); // DEBUG
  for (int residue = 0; residue < currentParm->Nres(); ++residue) {
    mprintf("DEBUG: Residue %i\n", residue+1);
    if (resmap.SetupResidue(*currentParm, residue) != 0) return Action::ERR;
    if (resmap.CheckBonds() != 0) return Action::ERR;
    resmap.DetermineAtomIDs();
    // Generate maps for symmetric atoms
    PairArray pairs;
    std::vector<bool> selected(resmap.Natom(), false);
    for (int atom1 = 0; atom1 < resmap.Natom(); atom1++) {
      if (!selected[atom1]) {
        if (resmap[atom1].Nduplicated() > 0) {
          selected[atom1] = true;
          Iarray symmatoms(1, atom1);
          for (int atom2 = atom1 + 1; atom2 < resmap.Natom(); atom2++) {
            if (resmap[atom1].Unique() == resmap[atom2].Unique()) {
              selected[atom2] = true;
              symmatoms.push_back(atom2);
            }
          } // End loop over atom2
          pairs.push_back( AtomPair(symmatoms) );
          mprintf("DEBUG:\t\tAtom ID %s is duplicated %u times:", resmap[atom1].Unique().c_str(),
                  symmatoms.size());
          for (Iarray::const_iterator sa = symmatoms.begin(); sa != symmatoms.end(); ++sa)
            mprintf(" %i", *sa + 1);
          mprintf("\n");
        }
      }
    } // End loop over atom1
    residues_.push_back( pairs );
    mprintf("DEBUG:\tNon-symmetric atoms:");
    for (int atom1 = 0; atom1 < resmap.Natom(); atom1++)
      if (!selected[atom1]) mprintf(" %i", atom1+1);
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
  for (std::vector<PairArray>::iterator res = residues_.begin(); res != residues_.end(); res++) 
  {
    // For each array of symmetric atoms, determine the lowest distance score
    for (PairArray::iterator pair = (*res).begin(); pair != (*res).end(); ++pair)
    {
      cost_matrix.Setup((*pair).AtomIndexes_.size(), (*pair).AtomIndexes_.size());
      for (Iarray::iterator tgtatom = (*pair).AtomIndexes_.begin(); 
                            tgtatom != (*pair).AtomIndexes_.end(); ++tgtatom)
      {
        for (Iarray::iterator refatom = (*pair).AtomIndexes_.begin();
                              refatom != (*pair).AtomIndexes_.end(); ++refatom)
        {
          double dist2 = DIST2_NoImage( RefFrame().XYZ(*refatom), currentFrame->XYZ(*tgtatom) );
          //mprintf("\t\t%i to %i: %f\n", *tgtatom + 1, *refatom + 1, dist2);
          cost_matrix.AddElement( dist2 ); 
        }
        //AMap_[ minrefatom ] = *tgtatom;
      }
      mprintf("\tRes %i", res - residues_.begin() + 1);
      cost_matrix.Print("Cost Matrix:");
      Hungarian HA(cost_matrix);
      Iarray resMap = HA.Optimize();
      // Fill in overall map
      Iarray::iterator rmap = resMap.begin();
      for (Iarray::iterator atmidx = (*pair).AtomIndexes_.begin();
                            atmidx != (*pair).AtomIndexes_.end(); ++atmidx, ++rmap)
      {
        AMap_[*atmidx] = (*pair).AtomIndexes_[*rmap];
        mprintf("\tAssigned atom %i to atom %i\n", *atmidx + 1, (*pair).AtomIndexes_[*rmap] + 1);
      }
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
