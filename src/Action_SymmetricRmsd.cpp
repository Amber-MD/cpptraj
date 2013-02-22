#include "Action_SymmetricRmsd.h"
#include "CpptrajStdio.h"
#include "AtomMap.h"

// CONSTRUCTOR
Action_SymmetricRmsd::Action_SymmetricRmsd() {}

void Action_SymmetricRmsd::Help() {
  mprintf("\t[<name>] <mask> [<refmask>] [out filename] [nofit | norotate] [mass]\n");
  mprintf("\t[ first | ref <filename> | refindex <#> |\n");
  mprintf("\treftraj <filename> [parm <parmname> | parmindex <#>] ]\n");
}

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
  // Per-res keywords
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

Action::RetType Action_SymmetricRmsd::Setup(Topology* currentParm, Topology** parmAddress) {
  // Target setup
  if (SetupRmsMask(*currentParm, "symmrmsd")) return Action::ERR;
  // Reference setup
  if (SetupRef(*currentParm, TgtMask().Nselected(), "symmrmsd"))
    return Action::ERR;
  // Check for symmetric atoms
  //AtomMask cMask = TgtMask();
  //cMask.ConvertMaskType(); // Convert to char mask
  AtomMap resmap;
  resmap.SetDebug(1); // DEBUG
  for (int residue = 0; residue < currentParm->Nres(); ++residue) {
    if (resmap.SetupResidue(*currentParm, residue) != 0) return Action::ERR;
    if (resmap.CheckBonds() != 0) return Action::ERR;
    resmap.DetermineAtomIDs();
  }

  return Action::OK;
}

Action::RetType Action_SymmetricRmsd::DoAction(int frameNum, Frame* currentFrame, 
                                               Frame** frameAddress) 
{
  // Perform any needed reference actions
  ActionRef( *currentFrame, Fit(), UseMass() );
  // Calculate RMSD
  double rmsdval = CalcRmsd( *currentFrame, SelectedRef(), RefTrans() );
  rmsd_->Add(frameNum, &rmsdval);
  return Action::OK;
}
