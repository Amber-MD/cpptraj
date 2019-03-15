#include "ParmFile.h"
#include "CpptrajStdio.h"
// All ParmIO classes go here
#include "Parm_Amber.h"
#include "Parm_PDB.h"
#include "Parm_Mol2.h"
#include "Parm_CharmmPsf.h"
#include "Parm_CIF.h"
#include "Parm_SDF.h"
#include "Parm_Tinker.h"
#include "Parm_Gromacs.h"
#include "BondSearch.h"

// ----- STATIC VARS / ROUTINES ------------------------------------------------
// NOTE: Must be in same order as ParmFormatType
const FileTypes::AllocToken ParmFile::PF_AllocArray[] = {
  { "Amber Topology",   0,                  Parm_Amber::WriteHelp, Parm_Amber::Alloc     },
  { "PDB File",         Parm_PDB::ReadHelp, 0,                     Parm_PDB::Alloc       },
  { "Mol2 File",        0,                  0,                     Parm_Mol2::Alloc      },
  { "Charmm PSF",       0,                  0,                     Parm_CharmmPsf::Alloc },
  { "CIF File",         0,                  0,                     Parm_CIF::Alloc       },
  { "Gromacs Topology", 0,                  0,                     Parm_Gromacs::Alloc   },
  { "SDF File",         0,                  0,                     Parm_SDF::Alloc       },
  { "Tinker File",      0,                  0,                     Parm_Tinker::Alloc    },
  { "Unknown Topology", 0,                  0,                     0                     }
};

const FileTypes::KeyToken ParmFile::PF_KeyArray[] = {
  { AMBERPARM,    "amber",   ".parm7" },
  { PDBFILE,      "pdb",     ".pdb"   },
  { MOL2FILE,     "mol2",    ".mol2"  },
  { CHARMMPSF,    "psf",     ".psf"   },
  { CIFFILE,      "cif",     ".cif"   },
  { GMXTOP,       "gromacs", ".top"   },
  { SDFFILE,      "sdf",     ".sdf"   },
  { TINKER,       "tinker",  ".arc"   },
  { TINKER,       "arc",     ".arc"   },
  { UNKNOWN_PARM, 0,         0        }
};

const FileTypes::KeyToken ParmFile::PF_WriteKeyArray[] = {
  { AMBERPARM,    "amber",   ".parm7" },
  { CHARMMPSF,    "psf",     ".psf"   },
  { UNKNOWN_PARM, 0,         0        }
};

// ParmFile::DetectFormat()
ParmIO* ParmFile::DetectFormat(FileName const& fname, ParmFormatType& ptype) {
  CpptrajFile file;
  if (file.SetupRead(fname, 0) == 0) {
    for (int i = 0; i < (int)UNKNOWN_PARM; i++) {
      ptype = (ParmFormatType)i;
      ParmIO* IO = (ParmIO*)FileTypes::AllocIO(PF_AllocArray, ptype, true );
      if (IO != 0 && IO->ID_ParmFormat( file ))
        return IO;
      delete IO;
    }
  }
  ptype = UNKNOWN_PARM;
  return 0;
}

// ParmFile::DetectFormat()
ParmFile::ParmFormatType ParmFile::DetectFormat(FileName const& fname) {
  ParmFormatType ptype;
  ParmIO* pio = DetectFormat(fname, ptype);
  delete pio;
  return ptype;
}

// ParmFile::ReadTopology()
int ParmFile::ReadTopology(Topology& Top, FileName const& fnameIn, 
                           ArgList const& argListIn, int debugIn) 
{
  if (fnameIn.empty()) {
    mprinterr("Error: No input topology name given.\n");
    return 1;
  }
  if (!File::Exists( fnameIn )) {
    File::ErrorMsg( fnameIn.full() );
    return 1;
  }
  parmName_ = fnameIn;
  ArgList argIn = argListIn;
  ParmFormatType pfType;
  ParmIO* parmio = 0;
  Top.SetDebug( debugIn );
  BondSearchType bstype;
  if (argIn.hasKey("nobondsearch"))
    bstype = SEARCH_NONE;
  else {
    bstype = SEARCH_REGULAR;
    std::string bsarg = argIn.GetStringKey("searchtype");
    if (!bsarg.empty()) {
      if (bsarg == "pairlist") {
        mprintf("\tWill use pair list to search for bonds between atoms.\n");
        mprintf("Warning: Searching for bonds via pair list is still experimental.\n");
        bstype = SEARCH_PAIRLIST;
      } else if (bsarg == "grid") {
        mprintf("\tWill use grid to search for bonds between residues.\n");
        mprintf("Warning: Searching for bonds via grid is still experimental.\n");
        bstype = SEARCH_GRID;
      } else
        mprintf("Warning: Unrecognized search type '%s'. Ignoring.\n", bsarg.c_str());
    }
  }
  double bondoffset = argIn.getKeyDouble("bondsearch", -1.0);
  bool molsearch = !argIn.hasKey("nomolsearch");
  if (!molsearch)
    mprintf("\tDisabling molecule search. Topology will have no molecule info.\n");
  // Only force bond search when 'bondsearch' is specified.
//  bool bondsearch = false;
//  if (argIn.Contains("bondsearch")) {
//    Top.SetOffset( argIn.getKeyDouble("bondsearch", -1.0) );
//    bondsearch = true;
//  }
  // 'as' keyword specifies a format
  std::string as_arg = argIn.GetStringKey("as");
  if (!as_arg.empty()) {
    pfType = (ParmFormatType)FileTypes::GetFormatFromString( PF_KeyArray, as_arg, UNKNOWN_PARM );
    if (pfType == UNKNOWN_PARM) {
      mprinterr("Error: Topology format '%s' not recognized.\n", as_arg.c_str());
      return 1;
    }
    parmio = (ParmIO*)FileTypes::AllocIO( PF_AllocArray, pfType, false );
  } else
    parmio = DetectFormat( parmName_, pfType );
  if (parmio == 0) {
    mprinterr("Error: Could not determine format of topology '%s'\n", parmName_.full());
    return 1;
  }
  mprintf("\tReading '%s' as %s\n", parmName_.full(),
          FileTypes::FormatDescription(PF_AllocArray, pfType) );
  parmio->SetDebug( debugIn );
  parmio->SetOffset( bondoffset );
  parmio->SetBondSearchType( bstype );
  if (parmio->processReadArgs(argIn)) return 1;
  int err = parmio->ReadParm( parmName_.Full(), Top);
  // Perform setup common to all parm files.
  if (err == 0) 
    err = Top.CommonSetup( molsearch );
  else
    mprinterr("Error reading topology file '%s'\n", parmName_.full());
  delete parmio;
  if (err > 0) return 1;
  return 0;
}

// ParmFile::WritePrefixTopology()
int ParmFile::WritePrefixTopology(Topology const& Top, std::string const& prefix,
                                  ParmFormatType fmtIn, int debugIn)
{
  if (prefix.empty()) return 1;
  FileName prefixName;
  if (Top.OriginalFilename().empty())
    prefixName.SetFileName_NoExpansion( prefix + ".parm7" );
  else
    prefixName.SetFileName_NoExpansion( prefix + "." + Top.OriginalFilename().Base() );
  return WriteTopology(Top, prefixName, ArgList(), fmtIn, debugIn);
}

// ParmFile::WriteTopology()
int ParmFile::WriteTopology(Topology const& Top, FileName const& fnameIn, 
                            ArgList const& argListIn, ParmFormatType fmtIn, int debugIn)
{
  parmName_ = fnameIn;
  ArgList argIn = argListIn;
  ParmFormatType fmt = fmtIn;
  if (fmt == UNKNOWN_PARM) {
    // Check arg list to see if format specified.
    fmt = (ParmFormatType)FileTypes::GetFormatFromArg(PF_WriteKeyArray, argIn, UNKNOWN_PARM);
    // If still UNKNOWN check file extension. Default to AMBERPARM
    if (fmt == UNKNOWN_PARM)
      fmt = (ParmFormatType)FileTypes::GetTypeFromExtension(PF_WriteKeyArray, parmName_.Ext(),
                                                            AMBERPARM);
  }
  ParmIO* parmio = (ParmIO*)FileTypes::AllocIO(PF_AllocArray, fmt, true);
  if (parmio == 0) return 1;
  parmio->SetDebug( debugIn );
  parmio->processWriteArgs( argIn );
  mprintf("\tWriting topology %i (%s) to '%s' with format %s\n", Top.Pindex(),
          Top.c_str(), parmName_.full(), FileTypes::FormatDescription(PF_AllocArray, fmt));
  int err = parmio->WriteParm( parmName_.Full(), Top );
  delete parmio;
  if (err != 0 ) {
    mprinterr("Error: writing topology file '%s'\n", parmName_.full());
    return 1;
  }
  return 0;
}
