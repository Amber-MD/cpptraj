#ifdef USE_SANDERLIB
#include <locale>
#include "Energy_Sander.h"
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

/** Set default input. */
void Energy_Sander::SetDefaultInput() {
  input_.extdiel = 80.0;
  input_.intdiel = 1.0;
  input_.rgbmax = 25.0;
  input_.saltcon = 0.0;
  input_.cut = 8.0;
  input_.dielc = 1.0;
  input_.rdt = 0.0;
  input_.fswitch = -1.0;

  input_.igb = 0;
  input_.alpb = 0;
  input_.gbsa = 0;
  input_.lj1264 = 0;
  input_.ipb = 0;
  input_.inp = 0;
  input_.vdwmeth = 1;
  input_.ew_type = 0;
  input_.ntb = 0;
  input_.ifqnt = 0;
  input_.jfastw = 1;
  input_.ntf = 2;
  input_.ntc = 2;
}
 
/** Set input for Sander. Determine which energy terms will be active. */
int Energy_Sander::SetInput(ArgList& argIn) {
  SetDefaultInput();
  input_.extdiel = argIn.getKeyDouble("extdiel", input_.extdiel);
  input_.intdiel = argIn.getKeyDouble("intdiel", input_.intdiel);
  input_.rgbmax = argIn.getKeyDouble("rgbmax", input_.rgbmax);
  input_.saltcon = argIn.getKeyDouble("saltcon", input_.saltcon);
  specified_cut_ = argIn.Contains("cut");
  input_.cut = argIn.getKeyDouble("cut", input_.cut);
  input_.dielc = argIn.getKeyDouble("dielc", input_.dielc);
  input_.rdt = argIn.getKeyDouble("rdt", input_.rdt);
  input_.fswitch = argIn.getKeyDouble("fswitch", input_.fswitch);

  specified_igb_ = argIn.Contains("igb");
  input_.igb = argIn.getKeyInt("igb", input_.igb);
  input_.alpb = argIn.getKeyInt("alpb", input_.alpb);
  input_.gbsa = argIn.getKeyInt("gbsa", input_.gbsa);
  input_.lj1264 = argIn.getKeyInt("lj1264", input_.lj1264);
  input_.ipb = argIn.getKeyInt("ipb", input_.ipb);
  input_.inp = argIn.getKeyInt("inp", input_.inp);
  input_.vdwmeth = argIn.getKeyInt("vdwmeth", input_.vdwmeth);
  input_.ew_type = argIn.getKeyInt("ew_type", input_.ew_type);
  specified_ntb_ = argIn.Contains("ntb");
  input_.ntb = argIn.getKeyInt("ntb", input_.ntb);
  input_.ifqnt = argIn.getKeyInt("ifqnt", input_.ifqnt);
  input_.jfastw = argIn.getKeyInt("jfastw", input_.jfastw);
  input_.ntf = argIn.getKeyInt("ntf", input_.ntf);
  input_.ntc = argIn.getKeyInt("ntc", input_.ntc);
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
  // Set some input options if not already specified.
  if (!fIn.BoxCrd().HasBox()) {
    if (!specified_cut_) input_.cut = 9999.0;
    if (!specified_igb_ && !specified_ntb_) {
      mprintf("Warning: No box info and igb/ntb not specified. Defaulting to igb=1, ntb=0\n");
      input_.igb = 1;
      input_.ntb = 0;
    }
  } else {
    if (!specified_cut_) input_.cut = 8.0;
    if (!specified_igb_ && !specified_ntb_) {
      mprintf("Warning: Box info and igb/ntb not specified. Defaulting to igb=0, ntb=1\n");
      input_.igb = 0;
      input_.ntb = 1;
    }
  }

  mprintf("    SANDER INPUT OPTIONS:\n"
          "\textdiel= %6.2g  intdiel= %6.2g  rgbmax=  %6.2g  saltcon= %6.2g  cut=  %6.2g\n"
          "\tdielc=   %6.2g  rdt=     %6.2g  fswitch= %6.2g  igb=     %6i  alpb= %6i\n"
          "\tgbsa=    %6i  ntb=     %6i  lj1264=  %6i  ipb=     %6i  inp=  %6i\n"
          "\tvdwmeth= %6i  ew_type= %6i  ifqnt=   %6i  jfastw=  %6i  ntf=  %6i\n"
          "\tntc=     %6i\n", input_.extdiel, input_.intdiel,
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
