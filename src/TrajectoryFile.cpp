#include "TrajectoryFile.h"
#include "CpptrajStdio.h"
// All TrajectoryIO classes go here
#include "Traj_AmberCoord.h"
#ifdef BINTRAJ
  #include "Traj_AmberNetcdf.h"
  #include "Traj_AmberRestartNC.h"
#endif
#include "Traj_PDBfile.h"
#include "Traj_AmberRestart.h"
#include "Traj_Mol2File.h"
#include "Traj_Conflib.h"
#include "Traj_CharmmDcd.h"
#include "Traj_Binpos.h"

// ----- STATIC VARS / ROUTINES ------------------------------------------------ 
const char TrajectoryFile::FORMAT_STRINGS[10][17] = {
"Unknown", "Amber NetCDF", "Amber NC Restart", "PDB", "Mol2", "Charmm DCD",
"BINPOS", "Amber Restart", "Amber Trajectory", "LMOD conflib" "\0" };

// TrajectoryFile::GetFormatFromArg()
/** Given an arglist, search for one of the file format keywords.
  * Default to AmberTraj if no arglist given or no keywords present. 
  */
TrajectoryFile::TrajFormatType TrajectoryFile::GetFormatFromArg(ArgList& argIn)
{
  TrajFormatType writeFormat = AMBERTRAJ;
  if      ( argIn.hasKey("pdb")      ) writeFormat=PDBFILE;
  else if ( argIn.hasKey("netcdf")   ) writeFormat=AMBERNETCDF;
  else if ( argIn.hasKey("restart")  ) writeFormat=AMBERRESTART;
  else if ( argIn.hasKey("ncrestart")) writeFormat=AMBERRESTARTNC;
  else if ( argIn.hasKey("restartnc")) writeFormat=AMBERRESTARTNC;
  else if ( argIn.hasKey("mol2")     ) writeFormat=MOL2FILE;
  else if ( argIn.hasKey("dcd")      ) writeFormat=CHARMMDCD;
  else if ( argIn.hasKey("charmm")   ) writeFormat=CHARMMDCD;
  else if ( argIn.hasKey("binpos")   ) writeFormat=BINPOS;
  return writeFormat;
}

// TrajectoryFile::GetExtensionForType()
std::string TrajectoryFile::GetExtensionForType(TrajFormatType typeIn) {
  std::string ext;
  switch (typeIn) {
    case PDBFILE : ext=".pdb"; break;
    case AMBERTRAJ: ext=".crd"; break;
    case AMBERNETCDF: ext=".nc"; break;
    case AMBERRESTART: ext=".rst7"; break;
    case AMBERRESTARTNC: ext=".ncrst"; break;
    case MOL2FILE: ext=".mol2"; break;
    case CHARMMDCD: ext=".dcd"; break;
    case BINPOS: ext=".binpos"; break;
    default: ext="";
  }
  return ext;
}

// TrajectoryFile::GetTypeFromExtension()
TrajectoryFile::TrajFormatType TrajectoryFile::GetTypeFromExtension( std::string const& extIn)
{
  if      ( extIn == ".nc" ) return AMBERNETCDF;
  else if ( extIn == ".ncrst" ) return AMBERRESTARTNC;
  else if ( extIn == ".pdb" ) return PDBFILE;
  else if ( extIn == ".mol2" ) return MOL2FILE;
  else if ( extIn == ".dcd" ) return CHARMMDCD;
  else if ( extIn == ".rst7" ) return AMBERRESTART;
  else if ( extIn == ".crd") return AMBERTRAJ;
  else if ( extIn == ".binpos") return BINPOS;
  // No entry for CONFLIB
  return UNKNOWN_TRAJ;
}

// TrajectoryFile::FormatString()
const char* TrajectoryFile::FormatString( TrajectoryFile::TrajFormatType tIn ) {
  return TrajectoryFile::FORMAT_STRINGS[ tIn ];
}
// -----------------------------------------------------------------------------

// CONSTRUCTOR
TrajectoryFile::TrajectoryFile() :
  debug_(0),
  trajParm_(0)
{}

// TrajectoryFile::SetDebug()
/** Set debug level. */
void TrajectoryFile::SetDebug(int debugIn) {
  debug_ = debugIn;
  if (debug_>0) mprintf("\tTrajectoryFile debug level set to %i\n",debug_);
}

// TrajectoryFile::SetFileNames()
void TrajectoryFile::SetFileNames(std::string const& full, std::string const& base) {
  trajName_ = full;
  baseName_ = base;
} 

int TrajectoryFile::SetTrajParm( Topology* tparmIn ) {
  if (tparmIn==0) {
    mprinterr("Error: TrajectoryFile: Parm file is NULL.\n");
    return 1;
  }
  trajParm_ = tparmIn;
  return 0;
}

// TrajectoryFile::SetupTrajectoryIO()
/** \param tformat Trajectory format to set up for.
  * \return TrajectoryIO class based on tformat. 0 on error.
  */
TrajectoryIO* TrajectoryFile::SetupTrajectoryIO(TrajFormatType tformat) {
  TrajectoryIO *tio = 0;
  switch ( tformat ) {
    case AMBERRESTART: tio = new Traj_AmberRestart(); break;
    case AMBERTRAJ   : tio = new Traj_AmberCoord();    break;
    case AMBERNETCDF :
#     ifdef BINTRAJ
      tio = new Traj_AmberNetcdf();
#     else
      mprinterr("    Error: Can not set up trajectory (%s):\n",trajName_.c_str());
      mprinterr("           Compiled without NETCDF support. Recompile with -DBINTRAJ\n");
#     endif
      break;
    case AMBERRESTARTNC :
#     ifdef BINTRAJ
      tio = new Traj_AmberRestartNC();
#     else
      mprinterr("    Error: Can not set up trajectory (%s):\n",trajName_.c_str());
      mprinterr("           Compiled without NETCDF support. Recompile with -DBINTRAJ\n");
#     endif
      break;
    case PDBFILE     : tio = new Traj_PDBfile();      break;
    case CONFLIB     : tio = new Traj_Conflib();      break;
    case MOL2FILE    : tio = new Traj_Mol2File();     break;
    case CHARMMDCD   : tio = new Traj_CharmmDcd();    break;
    case BINPOS      : tio = new Traj_Binpos();    break;
    default:
      return 0;
  }
  return tio;
}

// TrajectoryFile::DetectFormat()
TrajectoryIO* TrajectoryFile::DetectFormat(CpptrajFile& fileIn) {
  TrajectoryIO* tio = 0;
  for (int fmtidx = (int)UNKNOWN_TRAJ + 1; fmtidx != (int)NTRAJ; ++fmtidx) {
    tio = SetupTrajectoryIO( (TrajFormatType)fmtidx );
    if (tio != 0) {
      // Set TrajectoryIO file information from input file. 
      tio->TrajectoryIO::operator=( fileIn ); // NOTE: Should also set debug
      // Check if this file matches current TrajectoryIO format
      if (tio->ID_TrajFormat())
        break;
      // Not the right format; delete TrajectoryIO to try next format.
      delete tio;
      tio = 0;
    }
  }
  return tio;
}

