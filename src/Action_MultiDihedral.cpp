#include "Action_MultiDihedral.h"
#include "DataSet_double.h"
#include "CpptrajStdio.h"
#include "Constants.h"
#include "TorsionRoutines.h"
#include "StringRoutines.h" // convertToInteger

Action_MultiDihedral::Action_MultiDihedral() :
  minTorsion_(-180.0),
  debug_(0),
  outfile_(0),
  masterDSL_(0)
{}

void Action_MultiDihedral::Help() const {
  mprintf("\t[<name>] <dihedral types> [resrange <range>] [out <filename>] [range360]\n");
  mprintf("\t[dihtype <name>:<a0>:<a1>:<a2>:<a3>[:<offset>] ...]\n");
  DihedralSearch::OffsetHelp();
  //mprintf("\t[range360]\n");
  mprintf("\t<dihedral types> = ");
  DihedralSearch::ListKnownTypes();
  mprintf("  Calculate specified dihedral angle types for residues in given <range>.\n");
}

Action::RetType Action_MultiDihedral::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
  debug_ = debugIn;
  // Get keywords
  outfile_ = init.DFL().AddDataFile( actionArgs.GetStringKey("out"), actionArgs);
  if (actionArgs.hasKey("range360"))
    minTorsion_ = 0.0;
  else
    minTorsion_ = -180.0;
  std::string resrange_arg = actionArgs.GetStringKey("resrange");
  if (!resrange_arg.empty())
    if (resRange_.SetRange( resrange_arg )) return Action::ERR;
  // Search for known dihedral keywords
  dihSearch_.SearchForArgs(actionArgs);
  // Get custom dihedral arguments: dihtype <name>:<a0>:<a1>:<a2>:<a3>[:<offset>]
  std::string dihtype_arg = actionArgs.GetStringKey("dihtype");
  while (!dihtype_arg.empty()) {
    ArgList dihtype(dihtype_arg, ":");
    if (dihtype.Nargs() < 5) {
      mprinterr("Error: Malformed dihtype arg.\n");
      return Action::ERR;
    }
    int offset = 0;
    if (dihtype.Nargs() == 6) offset = convertToInteger(dihtype[5]);
    dihSearch_.SearchForNewType(offset,dihtype[1],dihtype[2],dihtype[3],dihtype[4], dihtype[0]);
    dihtype_arg = actionArgs.GetStringKey("dihtype");
  }
  // If no dihedral types yet selected, this will select all.
  dihSearch_.SearchForAll();

  // Setup DataSet(s) name
  dsetname_ = actionArgs.GetStringNext();

  mprintf("    MULTIDIHEDRAL: Calculating");
  dihSearch_.PrintTypes();
  if (!resRange_.Empty())
    mprintf(" dihedrals for residues in range %s\n", resRange_.RangeArg());
  else
    mprintf(" dihedrals for all solute residues.\n");
  if (!dsetname_.empty())
    mprintf("\tDataSet name: %s\n", dsetname_.c_str());
  if (outfile_ != 0) mprintf("\tOutput to %s\n", outfile_->DataFilename().base());
  if (minTorsion_ > -180.0) 
    mprintf("\tOutput range is 0 to 360 degrees.\n");
  else
    mprintf("\tOutput range is -180 to 180 degrees.\n");
  init.DSL().SetDataSetsPending(true);
  masterDSL_ = init.DslPtr();
  return Action::OK;
}

// Action_MultiDihedral::Setup();
Action::RetType Action_MultiDihedral::Setup(ActionSetup& setup) {
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
  // Search for specified dihedrals in each residue in the range
  if (dihSearch_.FindDihedrals(setup.Top(), actualRange))
    return Action::SKIP;
  mprintf("\tResRange=[%s]", resRange_.RangeArg());
  dihSearch_.PrintTypes();
  mprintf(", %i dihedrals.\n", dihSearch_.Ndihedrals());

  // Print selected dihedrals, set up DataSets
  data_.clear();
  if (dsetname_.empty())
    dsetname_ = masterDSL_->GenerateDefaultName("MDIH");
  for (DihedralSearch::mask_it dih = dihSearch_.begin();
                               dih != dihSearch_.end(); ++dih)
  {
    int resNum = dih->ResNum() + 1;
    // See if Dataset already present. FIXME should AddSet do this?
    MetaData md( dsetname_, dih->Name(), resNum );
    DataSet* ds = masterDSL_->CheckForSet(md);
    if (ds == 0) {
      // Create new DataSet
      md.SetScalarMode( MetaData::M_TORSION );
      md.SetScalarType( dih->Type() );
      ds = masterDSL_->AddSet( DataSet::DOUBLE, md );
      if (ds == 0) return Action::ERR;
      // Add to outfile
      if (outfile_ != 0)
        outfile_->AddDataSet( ds );
    }
    data_.push_back( ds ); 
    if (debug_ > 0) {
      mprintf("\tDIH [%s]:", ds->legend());
      mprintf(" :%i@%i",   setup.Top()[dih->A0()].ResNum()+1, dih->A0() + 1);
      mprintf(" :%i@%i",   setup.Top()[dih->A1()].ResNum()+1, dih->A1() + 1);
      mprintf(" :%i@%i",   setup.Top()[dih->A2()].ResNum()+1, dih->A2() + 1);
      mprintf(" :%i@%i\n", setup.Top()[dih->A3()].ResNum()+1, dih->A3() + 1);
    }
  }
  return Action::OK;
}

// Action_MultiDihedral::DoAction()
Action::RetType Action_MultiDihedral::DoAction(int frameNum, ActionFrame& frm) {
  std::vector<DataSet*>::const_iterator ds = data_.begin();
  for (DihedralSearch::mask_it dih = dihSearch_.begin();
                               dih != dihSearch_.end(); ++dih, ++ds)
  {
    double torsion = Torsion( frm.Frm().XYZ(dih->A0()),
                              frm.Frm().XYZ(dih->A1()),
                              frm.Frm().XYZ(dih->A2()),
                              frm.Frm().XYZ(dih->A3()) );
    torsion *= Constants::RADDEG;
    if (torsion < minTorsion_)
      torsion += 360.0;
    (*ds)->Add(frameNum, &torsion);
  }
  return Action::OK;
}
