#include "Exec_Zmatrix.h"
#include "CpptrajStdio.h"
#include "DataSet_Zmatrix.h"
#include "Structure/Zmatrix.h"

using namespace Cpptraj::Structure;

/** Get Zmatrix for specified molecule at specified frame number. */
int Exec_Zmatrix::getZmatrix(DataSet_Coords* CRD, int molnum, int frmidx,
                             std::string const& dsname, DataFile* outfile, CpptrajState& State)
const
{
 if (CRD == 0) {
    mprinterr("Error: No COORDS set specified.\n");
    return CpptrajState::ERR;
  }
  mprintf("\tGetting Z-matrix from set '%s'\n", CRD->legend());
  if (CRD->Size() < 1) {
    mprinterr("Error: Set '%s' is empty.\n", CRD->legend());
    return 1; 
  }

  DataSet_Zmatrix* zset = (DataSet_Zmatrix*)State.DSL().AddSet( DataSet::ZMATRIX, dsname, "ZMATRIX" );
  if (zset == 0) {
    mprinterr("Error: Could not allocate zmatrix set.\n");
    return 1; 
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
  if (debug_ > 0) zmatrix.print(); // DEBUG
  return errStat;
}

/** Create COORDS using Zmatrix. */
int Exec_Zmatrix::putZmatrix(DataSet_Coords* CRD, Topology* topIn,
                             std::string const& zsetname,
                             std::string const& dsname, CpptrajState& State)
const
{
  // Get zmatrix set
  DataSet_Zmatrix* zmatrix = (DataSet_Zmatrix*)State.DSL().FindSetOfType( zsetname, DataSet::ZMATRIX );
  if (zmatrix == 0) {
    mprinterr("Error: No zmatrix set with name '%s' found.\n", zsetname.c_str());
    return 1;
  }
  mprintf("\tCreating COORDS from Zmatrix set : %s\n", zmatrix->legend());

  // Create COORDS set
  DataSet_Coords* out = (DataSet_Coords*)State.DSL().AddSet(DataSet::COORDS, dsname, "ZMCRD");
  if (out == 0) {
    mprinterr("Error: Could not create COORDS set for assigning from Zmatrix.\n");
    return 1;
  }
  mprintf("\tOutput COORDS set: %s\n", out->legend());

  Topology* myTop = 0;
  if (topIn != 0) {
    mprintf("\tTopology: %s\n", topIn->c_str());
    myTop = topIn;
  } else if (CRD != 0) {
    mprintf("\tTopology from COORDS set: %s\n", CRD->legend());
    myTop = CRD->TopPtr();
  } else {
    mprinterr("Error: No COORDS set or Topology specified for 'zset'.\n");
    return 1;
  }
  CoordinateInfo myInfo;
  if (CRD != 0)
    myInfo = CRD->CoordsInfo();
  out->CoordsSetup( *myTop, myInfo);

  // Assign
  Frame frm = out->AllocateFrame();
  frm.ZeroCoords();
  int err = 0;
  if (debug_ > 0) {
    Zmatrix tmpZ = *(zmatrix->Zptr()); // TODO just pass debug level into Zmatrix?
    tmpZ.SetDebug( debug_ );
    err = tmpZ.SetToFrame( frm );
  } else {
    err = zmatrix->Zptr()->SetToFrame( frm );
  }
  if (err != 0) {
    mprinterr("Error: Zmatrix to Cartesian coords failed.\n");
    return 1;
  }
  out->SetCRD(0, frm);

  return 0;
}

// Exec_Zmatrix::Help()
void Exec_Zmatrix::Help() const
{
  mprintf("\t<COORDS set name> [name <output set name>]\n"
          "\t{ zset <input zmatrix set> [parm <top>|parmindex <#>] |\n"
          "\t  [molnum <mol#>] [frame <frame#>] [out <zmatrix file>] }\n"
          "  If 'zset' is specified, use Z-matrix to generate coordinates using\n"
          "  a specified topology or topology from specified COORDS set;\n"
          "  output is a new COORDS set.\n"
          "  Otherwise calculate Zmatrix for specified molecule/frame of\n"
          "  specified COORDS set; output is a Z-matrix set.\n");
}

// Exec_Zmatrix::Execute()
Exec::RetType Exec_Zmatrix::Execute(CpptrajState& State, ArgList& argIn)
{
  debug_ = State.Debug();
  std::string dsname = argIn.GetStringKey("name");
  std::string zsetname = argIn.GetStringKey("zset");
  int molnum = -1;
  int frmidx = -1;
  DataFile* outfile = 0;
  Topology* topIn = 0;
  if (zsetname.empty()) {
    // Calculating Z-matrix
    molnum = argIn.getKeyInt("molnum", 1) - 1;
    frmidx = argIn.getKeyInt("frame", 1) - 1;
    outfile = State.DFL().AddDataFile( argIn.GetStringKey("out"), argIn );
  } else {
    // Applying Zmatrix
    if (argIn.Contains("parm") || argIn.Contains("parmindex"))
      topIn = (Topology*)State.DSL().GetTopology(argIn);
  }
  // Get COORDS set name
  std::string setname = argIn.GetStringNext();
  // ----- No more args below here -----
  if (setname.empty()) {
    mprinterr("Error: %s: Specify COORDS dataset name.\n", argIn.Command());
    return CpptrajState::ERR;
  }
  DataSet_Coords* CRD = (DataSet_Coords*)State.DSL().FindSetOfGroup( setname, DataSet::COORDINATES );


  int errStat = 0;
  if (!zsetname.empty()) {
    errStat = putZmatrix(CRD, topIn, zsetname, dsname, State);
  } else {
    errStat = getZmatrix(CRD, molnum, frmidx, dsname, outfile, State);
  }

  if (errStat != 0) return CpptrajState::ERR;
  return CpptrajState::OK;
}
