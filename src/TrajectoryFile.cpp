#include "TrajectoryFile.h"
#include "CpptrajStdio.h"
// All TrajectoryIO classes go here
#include "Traj_AmberCoord.h"
#ifdef BINTRAJ
  #include "Traj_AmberNetcdf.h"
  #include "Traj_AmberRestartNC.h"
  #include "Traj_NcEnsemble.h"
#endif
#include "Traj_PDBfile.h"
#include "Traj_AmberRestart.h"
#include "Traj_Mol2File.h"
#include "Traj_Conflib.h"
#include "Traj_CharmmDcd.h"
#include "Traj_Binpos.h"
#include "Traj_GmxTrX.h"
#include "Traj_SQM.h"
#include "Traj_CIF.h"
#include "Traj_SDF.h"
#include "Traj_Tinker.h"
#include "Traj_CharmmCor.h"
#include "Traj_Gro.h"
#include "Traj_GmxXtc.h"
#include "Traj_CharmmRestart.h"
#include "Traj_XYZ.h"

// ----- STATIC VARS / ROUTINES ------------------------------------------------ 
// NOTE: Must be in same order as TrajFormatType
/** Static array containing traj allocators, optionally read/write help functions. */
const FileTypes::AllocToken TrajectoryFile::TF_AllocArray[] = {
# ifdef BINTRAJ
  { "Amber NetCDF",       Traj_AmberNetcdf::ReadHelp, Traj_AmberNetcdf::WriteHelp, Traj_AmberNetcdf::Alloc    },
  { "Amber NC Restart",   Traj_AmberRestartNC::ReadHelp, Traj_AmberRestartNC::WriteHelp, Traj_AmberRestartNC::Alloc },
# else
  { "Amber NetCDF",       0, 0, 0                          },
  { "Amber NC Restart",   0, 0, 0                          },
# endif
# if defined (ENABLE_SINGLE_ENSEMBLE) && defined (BINTRAJ)
  { "Amber NC Ensemble",  Traj_NcEnsemble::ReadHelp, Traj_NcEnsemble::WriteHelp, Traj_NcEnsemble::Alloc },
# else
  { "Amber NC Ensemble",  0, 0, 0                          },
# endif
  { "PDB",                0, Traj_PDBfile::WriteHelp, Traj_PDBfile::Alloc        },
  { "Mol2",               0, Traj_Mol2File::WriteHelp, Traj_Mol2File::Alloc       },
  { "CIF",                0, 0, Traj_CIF::Alloc            },
  { "Charmm DCD",         Traj_CharmmDcd::ReadHelp, Traj_CharmmDcd::WriteHelp, Traj_CharmmDcd::Alloc      },
  { "Gromacs TRX",        0, Traj_GmxTrX::WriteHelp, Traj_GmxTrX::Alloc         },
# ifdef NO_XDRFILE
  { "Gromacs XTC", 0, 0, 0                  },
# else
  { "Gromacs XTC", 0, Traj_GmxXtc::WriteHelp, Traj_GmxXtc::Alloc },
# endif
  { "BINPOS",             0, 0, Traj_Binpos::Alloc         },
  { "Amber Restart",      Traj_AmberRestart::ReadHelp, Traj_AmberRestart::WriteHelp, Traj_AmberRestart::Alloc   },
  { "GRO file",           0, 0, Traj_Gro::Alloc            },
  { "Tinker file",        0, 0, Traj_Tinker::Alloc         },
  { "Charmm COR",         0, Traj_CharmmCor::WriteHelp, Traj_CharmmCor::Alloc      },
  { "Charmm Restart",     Traj_CharmmRestart::ReadHelp, Traj_CharmmRestart::WriteHelp, Traj_CharmmRestart::Alloc },
  { "Amber Trajectory",   0, Traj_AmberCoord::WriteHelp, Traj_AmberCoord::Alloc     },
  { "SQM Input",          0, Traj_SQM::WriteHelp, Traj_SQM::Alloc            },
  { "SDF",                0, 0, Traj_SDF::Alloc            },
  { "XYZ",                0, Traj_XYZ::WriteHelp, Traj_XYZ::Alloc            },
  { "LMOD conflib",       0, 0, Traj_Conflib::Alloc        },
  { "Unknown trajectory", 0, 0, 0                          }
};

/** Input file formats. */
const FileTypes::KeyToken TrajectoryFile::TF_KeyArray[] = {
  { AMBERNETCDF,    "netcdf",    ".nc"      },
  { AMBERNETCDF,    "cdf",       ".nc"      },
  { AMBERRESTARTNC, "ncrestart", ".ncrst"   },
  { AMBERRESTARTNC, "restartnc", ".ncrst"   },
# ifdef ENABLE_SINGLE_ENSEMBLE
  { AMBERNCENSEMBLE,"ncensemble",".ncens"   },
# endif
  { PDBFILE,        "pdb",       ".pdb"     },
  { MOL2FILE,       "mol2",      ".mol2"    },
  { CIF,            "cif",       ".cif"     },
  { CHARMMDCD,      "dcd",       ".dcd"     },
  { CHARMMDCD,      "charmm",    ".dcd"     },
  { GMXTRX,         "trr",       ".trr"     },
  { GMXXTC,         "xtc",       ".xtc"     },
  { BINPOS,         "binpos",    ".binpos"  },
  { AMBERRESTART,   "restart",   ".rst7"    },
  { AMBERRESTART,   "restrt",    ".rst7"    },
  { AMBERRESTART,   "rest",      ".rst7"    },
  { GRO,            "gro",       ".gro"     },
  { TINKER,         "arc",       ".arc"     },
  { CHARMMCOR,      "cor",       ".cor"     },
  { CHARMMREST,     "charmmres", ".res"     },
  { AMBERTRAJ,      "crd",       ".crd"     },
  { CONFLIB,        "conflib",   ".conflib" },
  { SQM,            "sqm",       ".sqm"     },
  { SDF,            "sdf",       ".sdf"     },
  { XYZ,            "xyz",       ".xyz"     },
  { UNKNOWN_TRAJ,   0,           0          }
};

/** Output file formats. */
const FileTypes::KeyToken TrajectoryFile::TF_WriteKeyArray[] = {
  { AMBERNETCDF,    "netcdf",    ".nc"      },
  { AMBERNETCDF,    "cdf",       ".nc"      },
  { AMBERRESTARTNC, "ncrestart", ".ncrst"   },
  { AMBERRESTARTNC, "restartnc", ".ncrst"   },
# ifdef ENABLE_SINGLE_ENSEMBLE
  { AMBERNCENSEMBLE,"ncensemble",".ncens"   },
# endif
  { PDBFILE,        "pdb",       ".pdb"     },
  { MOL2FILE,       "mol2",      ".mol2"    },
  { CHARMMDCD,      "dcd",       ".dcd"     },
  { CHARMMDCD,      "charmm",    ".dcd"     },
  { GMXTRX,         "trr",       ".trr"     },
  { GMXXTC,         "xtc",       ".xtc"     },
  { BINPOS,         "binpos",    ".binpos"  },
  { AMBERRESTART,   "restart",   ".rst7"    },
  { AMBERRESTART,   "restrt",    ".rst7"    },
  { AMBERRESTART,   "rest",      ".rst7"    },
  { AMBERRESTART,   "rest",      ".rst"     },
  { AMBERTRAJ,      "crd",       ".crd"     },
  { SQM,            "sqm",       ".sqm"     },
  { XYZ,            "xyz",       ".xyz"     },
  { UNKNOWN_TRAJ,   0,           0          }
};

// -----------------------------------------------------------------------------

// TrajectoryFile::DetectFormat()
TrajectoryIO* TrajectoryFile::DetectFormat(FileName const& fname, TrajFormatType& ttype) {
  CpptrajFile file;
  if (file.SetupRead(fname, 0) == 0) {
    for (int i = 0; i < (int)UNKNOWN_TRAJ; i++) {
      ttype = (TrajFormatType)i;
      TrajectoryIO* IO = (TrajectoryIO*)FileTypes::AllocIO(TF_AllocArray, ttype, true );
      if (IO != 0 && IO->ID_TrajFormat( file ))
        return IO;
      delete IO;
    }
  }
  ttype = UNKNOWN_TRAJ;
  return 0;
}

// TrajectoryFile::DetectFormat()
TrajectoryFile::TrajFormatType TrajectoryFile::DetectFormat(FileName const& fname) {
  TrajFormatType ttype;
  TrajectoryIO* tio = DetectFormat(fname, ttype);
  delete tio;
  return ttype;
}

/** This version of DetectFormat() optionally takes a format keyword in
  * place of automatic type determination.
  */
TrajectoryIO* TrajectoryFile::DetectFormat(FileName const& fname, std::string const& fmtarg,
                                           TrajFormatType& ttype)
{
  TrajectoryIO* tio = 0;
  if (!fmtarg.empty()) {
    ttype = (TrajFormatType)FileTypes::GetFormatFromString( TF_KeyArray, fmtarg, UNKNOWN_TRAJ );
    if (ttype == UNKNOWN_TRAJ) {
      mprinterr("Error: Trajectory format '%s' is not recognized.\n", fmtarg.c_str());
      return 0;
    }
    tio = (TrajectoryIO*)FileTypes::AllocIO( TF_AllocArray, ttype, false );
  } else {
    if ( (tio = DetectFormat( fname, ttype )) == 0 ) {
      mprinterr("Error: Could not determine trajectory '%s' format.\n", fname.full());
      return 0;
    }
  }
  return tio;
}
