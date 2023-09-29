#include "Exec_Zmatrix.h"
#include "CpptrajStdio.h"
#include "DataSet_Zmatrix.h"
#include "Structure/Zmatrix.h"

using namespace Cpptraj::Structure;

// Exec_Zmatrix::Help()
void Exec_Zmatrix::Help() const
{

}

// Exec_Zmatrix::Execute()
Exec::RetType Exec_Zmatrix::Execute(CpptrajState& State, ArgList& argIn)
{
  int molnum = argIn.getKeyInt("molnum", 1) - 1;
  int frmidx = argIn.getKeyInt("frame", 1) - 1;
  std::string dsname = argIn.GetStringKey("name");
  DataFile* outfile = State.DFL().AddDataFile( argIn.GetStringKey("out"), argIn );

  std::string setname = argIn.GetStringNext();
  if (setname.empty()) {
    mprinterr("Error: %s: Specify COORDS dataset name.\n", argIn.Command());
    return CpptrajState::ERR;
  }
  DataSet_Coords* CRD = (DataSet_Coords*)State.DSL().FindSetOfGroup( setname, DataSet::COORDINATES );
  if (CRD == 0) {
    mprinterr("Error: %s: No COORDS set with name %s found.\n", argIn.Command(), setname.c_str());
    return CpptrajState::ERR;
  }
  mprintf("\tUsing set '%s'\n", CRD->legend());
  if (CRD->Size() < 1) {
    mprinterr("Error: Set '%s' is empty.\n");
    return CpptrajState::ERR;
  }

  DataSet_Zmatrix* zset = (DataSet_Zmatrix*)State.DSL().AddSet( DataSet::ZMATRIX, dsname, "ZMATRIX" );
  if (zset == 0) {
    mprinterr("Error: Could not allocate zmatrix set.\n");
    return CpptrajState::ERR;
  }
  zset->SetDim(Dimension::X, Dimension(1, 1, "IC"));
  if (outfile != 0)
    outfile->AddDataSet( zset );

  mprintf("\tOutput set  : %s\n", zset->legend());
  mprintf("\tMolecule    : %i\n", molnum + 1 );
  mprintf("\tFrame       : %i\n", frmidx + 1 );
  if (outfile != 0)
    mprintf("\tOutput file : %s\n", outfile->DataFilename().full());

  Frame frmIn = CRD->AllocateFrame();
  CRD->GetFrame( frmidx, frmIn );

  Zmatrix& zmatrix = *(zset->Zptr());

  zmatrix.SetDebug( State.Debug() );
  int errStat = zmatrix.SetFromFrame( frmIn, CRD->Top(), molnum );
  zmatrix.print(); // DEBUG
  if (errStat != 0) return CpptrajState::ERR;
  return CpptrajState::OK;
}
