#include "Action_MultiPucker.h"
#include "CpptrajStdio.h"

using namespace Cpptraj;

Action_MultiPucker::Action_MultiPucker() :
  outfile_(0),
  masterDSL_(0)
{}

// Action_MultiPucker::Help()
void Action_MultiPucker::Help() const {

}

// Action_MultiPucker::Init()
Action::RetType Action_MultiPucker::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  // Get Keywords
  outfile_ = init.DFL().AddDataFile( actionArgs.GetStringKey("out"), actionArgs);
  std::string resrange_arg = actionArgs.GetStringKey("resrange");
  if (!resrange_arg.empty())
    if (resRange_.SetRange( resrange_arg )) return Action::ERR;
  // Search for known pucker keywords
  if (puckerSearch_.SearchForArgs( actionArgs )) return Action::ERR;
  // Get custom pucker args
  if (puckerSearch_.SearchForNewTypeArgs( actionArgs )) return Action::ERR;
  // If no pucker types are yet selected, this will select all.
  puckerSearch_.SearchForAll();

  // Setup DataSet(s) name
  dsetname_ = actionArgs.GetStringNext();

  mprintf("    MULTIPUCKER: Calculating");
  puckerSearch_.PrintTypes();
  if (!resRange_.Empty())
    mprintf(" puckers for residues in range %s\n", resRange_.RangeArg());
  else
    mprintf(" puckers for all solute residues.\n");
  if (!dsetname_.empty())
    mprintf("\tDataSet name: %s\n", dsetname_.c_str());
  if (outfile_ != 0) mprintf("\tOutput to %s\n", outfile_->DataFilename().base());
  //if (minTorsion_ > -180.0) 
  //  mprintf("\tOutput range is 0 to 360 degrees.\n");
  //else
  //  mprintf("\tOutput range is -180 to 180 degrees.\n");
  init.DSL().SetDataSetsPending(true);
  masterDSL_ = init.DslPtr();
  return Action::OK;
}

// Action_MultiPucker::Setup()
Action::RetType Action_MultiPucker::Setup(ActionSetup& setup)
{
  Range actualRange;
  // If range is empty (i.e. no resrange arg given) look through all 
  // solute residues.
  if (resRange_.Empty())
    actualRange = setup.Top().SoluteResidues();
  else {
    // If user range specified, create new range shifted by -1 since internal
    // resnums start from 0.
    actualRange = resRange_;
    actualRange.ShiftBy(-1);
  }
  // Exit if no residues specified
  if (actualRange.Empty()) {
    mprinterr("Error: No residues specified for %s\n",setup.Top().c_str());
    return Action::ERR;
  }
  // Search for specified puckers in each residue in the range
  if (puckerSearch_.FindPuckers(setup.Top(), actualRange))
    return Action::SKIP;
  mprintf("\tResRange=[%s]", resRange_.RangeArg());
  puckerSearch_.PrintTypes();
  mprintf(", %u puckers.\n", puckerSearch_.Npuckers());

  // Print selected puckers, set up DataSets
  data_.clear();
  if (dsetname_.empty())
    dsetname_ = masterDSL_->GenerateDefaultName("MPUCKER");
  for (Pucker::PuckerSearch::mask_it pucker = puckerSearch_.begin();
                                     pucker != puckerSearch_.end(); ++pucker)
  {
    int resNum = pucker->ResNum() + 1;
    // See if Dataset already present. FIXME should AddSet do this?
    MetaData md( dsetname_, pucker->Name(), resNum );
    DataSet* ds = masterDSL_->CheckForSet(md);
    if (ds == 0) {
      // Create new DataSet
      md.SetScalarMode( MetaData::M_PUCKER );
      md.SetScalarType( MetaData::PUCKER   ); // TODO pucker types
      ds = masterDSL_->AddSet( DataSet::DOUBLE, md );
      if (ds == 0) return Action::ERR;
      // Add to outfile
      if (outfile_ != 0)
        outfile_->AddDataSet( ds );
    }
    data_.push_back( ds ); 
    //if (debug_ > 0) {
      mprintf("\tPUCKER [%s]: %s", ds->legend(), pucker->PuckerMaskString(setup.Top()).c_str());
    //}
  }
  return Action::OK;
}

// Action_MultiPucker::DoAction()
Action::RetType Action_MultiPucker::DoAction(int frameNum, ActionFrame& frm)
{

}
