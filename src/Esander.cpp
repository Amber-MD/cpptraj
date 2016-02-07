#ifdef USE_SANDERLIB
#include <locale>
#include "Esander.h"
#include "CpptrajStdio.h"

Energy_Sander::~Energy_Sander() { if (is_setup()) sander_cleanup(); }

const char* Energy_Sander::Estring_[] = {
  "Total", "VDW", "Elec", "GB", "Bond", "Angle", "Dihedral", "VDW14", "Elec14",
  "Constraint", "Polar", "Hbond", "Surf", "SCF", "Dispersion", "DVDL", "Angle_UB",
  "Improper", "CMap", "EMap", "LES", "NOE", "PB", "RISM", "CT", "aMD_Boost", 0
};

std::string Energy_Sander::Easpect(Etype typeIn) {
  if (typeIn == N_ENERGYTYPES) return std::string();
  std::locale loc;
  std::string aspect( Estring_[typeIn] );
  for (std::string::iterator it = aspect.begin(); it != aspect.end(); ++it)
    *it = std::tolower( *it, loc );
  return aspect;
}

double Energy_Sander::Energy(Etype typeIn) const {
  switch (typeIn) {
    case TOTAL: return energy_.tot;
    case VDW: return energy_.vdw;
    case ELEC: return energy_.elec;
    case GB: return energy_.gb;
    case BOND: return energy_.bond;
    case ANGLE: return energy_.angle;
    case DIHEDRAL: return energy_.dihedral;
    case VDW14: return energy_.vdw_14;
    case ELEC14: return energy_.elec_14;
    case CONSTRAINT: return energy_.constraint;
    case POLAR: return energy_.polar;
    case HBOND: return energy_.hbond;
    case SURF: return energy_.surf;
    case SCF: return energy_.scf;
    case DISP: return energy_.disp;
    case DVDL: return energy_.dvdl;
    case ANGLE_UB: return energy_.angle_ub;
    case IMP: return energy_.imp;
    case CMAP: return energy_.cmap;
    case EMAP: return energy_.emap;
    case LES: return energy_.les;
    case NOE: return energy_.noe;
    case PB: return energy_.pb;
    case RISM: return energy_.rism;
    case CT: return energy_.ct;
    case AMD_BOOST: return energy_.amd_boost;
    case N_ENERGYTYPES: break;
  }
  return 0.0; 
}

const double* Energy_Sander::Eptr(Etype typeIn) const {
  switch (typeIn) {
    case TOTAL: return &energy_.tot;
    case VDW: return &energy_.vdw;
    case ELEC: return &energy_.elec;
    case GB: return &energy_.gb;
    case BOND: return &energy_.bond;
    case ANGLE: return &energy_.angle;
    case DIHEDRAL: return &energy_.dihedral;
    case VDW14: return &energy_.vdw_14;
    case ELEC14: return &energy_.elec_14;
    case CONSTRAINT: return &energy_.constraint;
    case POLAR: return &energy_.polar;
    case HBOND: return &energy_.hbond;
    case SURF: return &energy_.surf;
    case SCF: return &energy_.scf;
    case DISP: return &energy_.disp;
    case DVDL: return &energy_.dvdl;
    case ANGLE_UB: return &energy_.angle_ub;
    case IMP: return &energy_.imp;
    case CMAP: return &energy_.cmap;
    case EMAP: return &energy_.emap;
    case LES: return &energy_.les;
    case NOE: return &energy_.noe;
    case PB: return &energy_.pb;
    case RISM: return &energy_.rism;
    case CT: return &energy_.ct;
    case AMD_BOOST: return &energy_.amd_boost;
    case N_ENERGYTYPES: break;
  }
  return 0;
}

// TODO const frame and top?
int Energy_Sander::Initialize(Topology const& topIn, Frame& fIn) {
  if (fIn.Natom() != topIn.Natom()) return 3;
  if (is_setup()) sander_cleanup();
  // FIXME: requires file name be set for now
  top_filename_ = topIn.OriginalFilename();
  top_pindex_ = topIn.Pindex();
  if (top_filename_.empty()) return 2;
  mprintf("DEBUG: Topology filename= '%s'\n", top_filename_.full());
  input_.extdiel = 80.0;
  input_.intdiel = 1.0;
  input_.rgbmax = 25.0;
  input_.saltcon = 0.0;
  // cut set below
  input_.dielc = 1.0;
  input_.rdt = 0.0;
  input_.fswitch = -1.0;

  // igb set below
  input_.alpb = 0;
  input_.gbsa = 0;
  input_.lj1264 = 0;
  input_.ipb = 0;
  input_.inp = 0;
  input_.vdwmeth = 1;
  input_.ew_type = 0;
  // ntb set below
  input_.ifqnt = 0;
  input_.jfastw = 1;
  input_.ntf = 2;
  input_.ntc = 2;
  if (!fIn.BoxCrd().HasBox()) {
    input_.cut = 9999.0;
    input_.igb = 1;
    input_.ntb = 0;
  } else {
    input_.cut = 8.0;
    input_.igb = 0;
    input_.ntb = 1;
  }

  mprintf("\textdiel= %g  intdiel= %g  rgbmax= %g  saltcon= %g  cut= %g\n"
          "\tdielc= %g  rdt= %g  fswitch= %g  igb= %i  alpb= %i  gbsa= %i  ntb= %i\n"
          "\tlj1264= %i  ipb= %i  inp= %i  vdwmeth= %i  ew_type= %i\n"
          "\tifqnt= %i  jfastw= %i  ntf= %i  ntc= %i\n", input_.extdiel, input_.intdiel,
          input_.rgbmax, input_.saltcon, input_.cut, input_.dielc,
          input_.rdt, input_.fswitch, input_.igb, input_.alpb, input_.gbsa, 
          input_.ntb, input_.lj1264, input_.ipb, input_.inp, input_.vdwmeth,
          input_.ew_type, input_.ifqnt, input_.jfastw, input_.ntf, input_.ntc);

  forces_.resize( fIn.size(), 0.0 );
  
  return sander_setup_mm(top_filename_.full(), fIn.xAddress(), fIn.bAddress(), &input_);
}

int Energy_Sander::CalcEnergy(Frame& fIn) {
  if (!is_setup()) return 1;

  set_positions( fIn.xAddress() );
  set_box( fIn.BoxCrd().BoxX(), fIn.BoxCrd().BoxY(), fIn.BoxCrd().BoxZ(),
           fIn.BoxCrd().Alpha(), fIn.BoxCrd().Beta(), fIn.BoxCrd().Gamma() );
  energy_forces( &energy_, &(forces_[0]) );
  return 0;
};
#endif
