//#include <algorithm> // next_permuation
//#include <list>
#include <cmath>
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
  resmap.SetDebug(0); // DEBUG
  for (int residue = 0; residue < currentParm->Nres(); ++residue) {
    mprintf("DEBUG: Residue %i\n", residue+1);
    if (resmap.SetupResidue(*currentParm, residue) != 0) return Action::ERR;
    if (resmap.CheckBonds() != 0) return Action::ERR;
    resmap.DetermineAtomIDs();
    // Generate maps for symmetric atoms
    residues_.push_back( SymRes() ); // DEBUG
    std::vector<bool> selected(resmap.Natom(), false);
    for (int atom1 = 0; atom1 < resmap.Natom(); atom1++) {
      // DEBUG
      residues_.back().NonSymAtoms_.push_back( atom1 + currentParm->Res(residue).FirstAtom() );
      // DEBUG
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
          residues_.back().SymMasks_.push_back( symmatoms );
          mprintf("DEBUG:\t\tAtom ID %s is duplicated %u times:", resmap[atom1].Unique().c_str(),
                  symmatoms.size());
          //do {
          //for (unsigned int iteration = 0; iteration < symmatoms.size(); ++iteration) {
            for (Iarray::const_iterator sa = symmatoms.begin(); sa != symmatoms.end(); ++sa)
              mprintf(" %i", *sa + 1);
            mprintf("\n");
          //  symmatoms.push_front( symmatoms.back() );
          //  symmatoms.pop_back();
          //}
          //} while (std::next_permutation(symmatoms.begin(), symmatoms.end()));
        }
      }
    } // End loop over atom1
    mprintf("DEBUG:\tNon-symmetric atoms:");
    for (int atom1 = 0; atom1 < resmap.Natom(); atom1++)
      if (!selected[atom1]) mprintf(" %i", atom1+1);
    mprintf("\n");
  } // End loop over residue

  return Action::OK;
}

Action::RetType Action_SymmetricRmsd::DoAction(int frameNum, Frame* currentFrame, 
                                               Frame** frameAddress) 
{
  // Perform any needed reference actions
  ActionRef( *currentFrame, Fit(), UseMass() );
  // Calculate initial best-fit RMSD
  if (Fit())
    CalcRmsd( *currentFrame, SelectedRef(), RefTrans() );
  // Correct RMSD for symmetry
  double rms_return = 0.0;
  int total = 0;
  // FIXME: currently calcs all atom RMSD
  for (std::vector<SymRes>::iterator res = residues_.begin(); res != residues_.end(); res++) 
  {
    // First get RMSD of non-symmetric atoms.
    for (Iarray::const_iterator atom = (*res).NonSymAtoms_.begin();
                                atom != (*res).NonSymAtoms_.end(); ++atom)
    { 
      Vec3 diff( RefFrame().XYZ(*atom) );
      diff -= Vec3( currentFrame->XYZ(*atom) );
      rms_return += diff.Magnitude2();
      ++total;
    }
    // For each array of symmetric atoms, determine the lowest distance score
    for (std::vector<Iarray>::iterator symmask = (*res).SymMasks_.begin();
                                       symmask != (*res).SymMasks_.end(); ++symmask)
    {
      for (Iarray::iterator tgtatom = (*symmask).begin(); 
                            tgtatom != (*symmask).end(); ++tgtatom)
      {
        for (Iarray::iterator refatom = (*symmask).begin();
                              refatom != (*symmask).end(); ++refatom)
        {
          Vec3 diff( RefFrame().XYZ(*refatom) );
          diff -= Vec3( currentFrame->XYZ(*tgtatom) );
          mprintf("\t\t%i to %i: %f\n", *tgtatom + 1, *refatom + 1, diff.Magnitude2());
        }
      }
    }
  }
  rms_return = sqrt(rms_return / (double)total);
  //rmsd_->Add(frameNum, &rmsdval);
  rmsd_->Add(frameNum, &rms_return);
  return Action::OK;
}
