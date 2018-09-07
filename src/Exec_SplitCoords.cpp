#include "Exec_SplitCoords.h"
#include "CpptrajStdio.h"
#include "DataSet_Coords_CRD.h"
#include "ProgressBar.h"

// Exec_SplitCoords::Help()
void Exec_SplitCoords::Help() const
{
  mprintf("\t<crd set> name <output set name>\n");
}

// Exec_SplitCoords::Execute()
Exec::RetType Exec_SplitCoords::Execute(CpptrajState& State, ArgList& argIn)
{
  std::string setname = argIn.GetStringNext();
  if (setname.empty()) {
    mprinterr("Error: %s: Specify COORDS dataset name.\n", argIn.Command());
    return CpptrajState::ERR;
  }
  DataSet_Coords* CRD = (DataSet_Coords*)State.DSL().FindCoordsSet( setname );
  if (CRD == 0) {
    mprinterr("Error: %s: No COORDS set with name %s found.\n", argIn.Command(), setname.c_str());
    return CpptrajState::ERR;
  }
  mprintf("\tUsing set '%s'\n", CRD->legend());
  mprintf("\tSplitting by molecule.\n");
  // Output COORDS set
  std::string dsname = argIn.GetStringKey("name");
  if (dsname.empty()) {
    mprinterr("Error: Must specify output COORDS name.\n");
    return CpptrajState::ERR;
  }
  DataSet_Coords_CRD* OUT = (DataSet_Coords_CRD*)
                            State.DSL().AddSet(DataSet::COORDS, MetaData(dsname));
  if (OUT == 0) return CpptrajState::ERR;
  mprintf("\tOutput to COORDS set '%s'\n", OUT->legend());
  // In order for this to work, currently must have all molecules be the same
  // size. TODO check that residue names match as well?
  Topology const& topIn = CRD->Top();
  if (topIn.Nmol() < 2) {
    mprinterr("Error: Topology for '%s' has less than 2 molecules.\n", CRD->legend());
    return CpptrajState::ERR;
  }
  int molNatoms = -1;
  int molNres = -1;
  for (Topology::mol_iterator mol = topIn.MolStart(); mol != topIn.MolEnd(); ++mol)
  {
    int nres = topIn[mol->BeginAtom()].ResNum() - topIn[mol->EndAtom()-1].ResNum() + 1;
    if (molNatoms == -1) {
      molNatoms = mol->NumAtoms();
      molNres = nres;
    } else if (molNatoms != mol->NumAtoms()) {
      mprinterr("Error: Molecule %u has different number of atoms (%i) than first molecule (%i)\n",
                mol - topIn.MolStart() + 1, mol->NumAtoms(), molNatoms);
      return CpptrajState::ERR;
    } else if (molNres != nres) {
      mprinterr("Error: Molecule %u has different number of residues (%i) than first molecule (%i)\n",
                mol - topIn.MolStart() + 1, nres, molNres);
      return CpptrajState::ERR;
    } 
  }
  // Create a mask for each input molecule.
  std::vector<AtomMask> Masks;
  Masks.reserve( topIn.Nmol() );
  Topology* topOut = 0;
  for (Topology::mol_iterator mol = topIn.MolStart(); mol != topIn.MolEnd(); ++mol)
  {
    Masks.push_back( AtomMask(mol->BeginAtom(), mol->EndAtom()) );
    // Set total number of atoms
    Masks.back().SetNatoms( topIn.Natom() );
    // First time around set up the output topology.
    if (topOut == 0) {
      topOut = topIn.modifyStateByMask( Masks.back() );
      if (topOut == 0) return CpptrajState::ERR;
    }
  }
  topOut->Brief("Split topology");
  // Set up output COORDS
  if (OUT->CoordsSetup( *topOut, CRD->CoordsInfo() )) return CpptrajState::ERR;
  OUT->Allocate( DataSet::SizeArray(1, CRD->Size() * topIn.Nmol()) );
  // OUT now has a copy of topOut, so it is no longer needed here.
  delete topOut;
  // Set up output frame
  Frame frameOut;
  frameOut.SetupFrameV(OUT->Top().Atoms(), OUT->CoordsInfo());
  // Loop over all input frames.
  Frame frameIn = CRD->AllocateFrame();
  ProgressBar progress( CRD->Size() * topIn.Nmol() );
  int idx = 0;
  for (unsigned int frm = 0; frm != CRD->Size(); frm++)
  {
    CRD->GetFrame(frm, frameIn);
    for (std::vector<AtomMask>::const_iterator mask = Masks.begin(); mask != Masks.end(); ++mask)
    {
      progress.Update( idx++ );
      frameOut.SetFrame( frameIn, *mask );
      OUT->AddFrame( frameOut );
    }
  }
  mprintf("\t'%s' : %zu frames.\n", OUT->legend(), OUT->Size());
  return CpptrajState::OK;
}
