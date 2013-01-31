#include "Action_MultiDihedral.h"
#include "DataSet_double.h"
#include "CpptrajStdio.h"
#include "Constants.h"
#include "TorsionRoutines.h"

Action_MultiDihedral::Action_MultiDihedral() :
  useMass_(false),
  range360_(false),
  findPHI_(false),
  findPSI_(false),
  outfile_(0)
{}

void Action_MultiDihedral::Help() {
  mprintf("\t[<name>] [phi] [psi] [resrange <range>] [mass] [out <filename>]\n");
  mprintf("\t[range360]\n");
}

Action::RetType Action_MultiDihedral::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  // Get keywords
  outfile_ = DFL->AddDataFile( actionArgs.GetStringKey("out"), actionArgs);
  useMass_ = actionArgs.hasKey("mass");
  range360_ = actionArgs.hasKey("range360");
  std::string resrange_arg = actionArgs.GetStringKey("resrange");
  if (!resrange_arg.empty())
    if (resRange_.SetRange( resrange_arg )) return Action::ERR;
  findPHI_ = actionArgs.hasKey("phi");
  findPSI_ = actionArgs.hasKey("psi");
  if (!findPHI_ && !findPSI_) { 
    findPHI_ = true;
    findPSI_ = true;
  }

  // Setup DataSet(s) name
  dsetname_ = actionArgs.GetStringNext();
  if (dsetname_.empty())
    dsetname_ = DSL->GenerateDefaultName("MDIH");

  mprintf("    MULTIDIHEDRAL: Calculating");
  if (findPHI_) mprintf(" phi");
  if (findPSI_) mprintf(" psi");
  if (!resRange_.Empty())
    mprintf(" dihedrals for residues in range %s\n", resRange_.RangeArg());
  else
    mprintf(" dihedrals for all residues.\n");
  mprintf("\tDataSet name: %s\n", dsetname_.c_str());
  if (outfile_ != 0) mprintf("\tOutput to %s\n", outfile_->Filename());
  if (useMass_) mprintf("\tMass-weighted\n");
  if (range360_) 
    mprintf("\tRange 0-360 deg.\n");
  else
    mprintf("\tRange -180-180 deg.\n");
  masterDSL_ = DSL;
  return Action::OK;
}

void Action_MultiDihedral::FindDihedralAtoms(Topology* currentParm, int resIn, int offset,
                                             NameType const& a1, NameType const& a2,
                                             NameType const& a3, NameType const& a4,
                                             std::string const& aspectIn)
{
  int firstAtomRes = resIn;
  int lastAtomRes = resIn;
  if (offset == -1)
    --firstAtomRes;
  else if (offset == 1)
    ++lastAtomRes;
  int atom1 = currentParm->FindAtomInResidue(firstAtomRes, a1);
  if (atom1 == -1) return;
  int atom2 = currentParm->FindAtomInResidue(resIn, a2);
  if (atom2 == -1) return;
  int atom3 = currentParm->FindAtomInResidue(resIn, a3);
  if (atom3 == -1) return;
  int atom4 = currentParm->FindAtomInResidue(lastAtomRes, a4);
  if (atom4 == -1) return;
  // All atoms found at this point.
  maskAtoms_.push_back(atom1);
  maskAtoms_.push_back(atom2);
  maskAtoms_.push_back(atom3);
  maskAtoms_.push_back(atom4);
  // See if Dataset already present
  DataSet* ds = masterDSL_->GetSet(dsetname_, resIn+1, aspectIn);
  if (ds == 0) { 
    // Create new DataSet
    ds = masterDSL_->AddSetIdxAspect( DataSet::DOUBLE, dsetname_, resIn+1, aspectIn);
    // Add to outfile
    if (outfile_ != 0)
      outfile_->AddSet( ds );
  }
  if (ds != 0)
    data_.push_back( ds );
}

Action::RetType Action_MultiDihedral::Setup(Topology* currentParm, Topology** parmAddress) {
  Range actualRange;
  // If range is empty (i.e. no resrange arg given) look through all 
  // solute residues.
  if (resRange_.Empty())
    actualRange.SetRange(0, currentParm->FinalSoluteRes());
  else {
    // If user range specified, create new range shifted by -1 since internal
    // resnums start from 0.
    actualRange = resRange_;
    actualRange.ShiftBy(-1);
  }
  // Exit if no residues specified
  if (actualRange.Empty()) {
    mprinterr("Error: multidihedral: No residues specified for %s\n",currentParm->c_str());
    return Action::ERR;
  }
  // Search for specified dihedrals in each residue in the range
  maskAtoms_.clear();
  data_.clear();
  for (Range::const_iterator res = actualRange.begin(); res != actualRange.end(); ++res)
  {
    // PHI: C0-N1-CA1-C1
    if (findPHI_)
      FindDihedralAtoms(currentParm, *res, -1, "C", "N", "CA", "C", "phi");
    // PSI: N0-CA0-C0-N1
    if (findPSI_)
      FindDihedralAtoms(currentParm, *res, 1, "N", "CA", "C", "N", "psi");
  }
  if (maskAtoms_.empty()) {
    mprintf("Warning: No dihedrals selected.\n");
    return Action::ERR;
  }
  // Print selected dihedrals
  for (std::vector<int>::iterator atom = maskAtoms_.begin();
                                  atom != maskAtoms_.end(); ++atom)
  {
    mprintf("\tDIH: %s", currentParm->TruncResAtomName(*(atom++)).c_str());
    mprintf("-%s", currentParm->TruncResAtomName(*(atom++)).c_str());
    mprintf("-%s", currentParm->TruncResAtomName(*(atom++)).c_str());
    mprintf("-%s\n", currentParm->TruncResAtomName(*atom).c_str());
  }
  // Print created DataSets
  for (std::vector<DataSet*>::iterator ds = data_.begin(); ds != data_.end(); ++ds)
    mprintf("\tDS: %s\n", (*ds)->Legend().c_str());
  return Action::OK;
}

Action::RetType Action_MultiDihedral::DoAction(int frameNum, Frame* currentFrame, 
                                               Frame** frameAddress)
{
  std::vector<DataSet*>::iterator ds = data_.begin();
  for (std::vector<int>::iterator atom = maskAtoms_.begin();
                                  atom != maskAtoms_.end(); atom += 4, ++ds)
  {
    double torsion = Torsion( currentFrame->XYZ(*atom),
                              currentFrame->XYZ(*(atom+1)),
                              currentFrame->XYZ(*(atom+2)),
                              currentFrame->XYZ(*(atom+3)) );
    torsion *= RADDEG;
    (*ds)->Add(frameNum, &torsion);
  }
  return Action::OK;
}

void Action_MultiDihedral::Print() {

}
