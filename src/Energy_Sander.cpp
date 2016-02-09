#ifdef USE_SANDERLIB
#include <locale>
#include "Energy_Sander.h"
#include "CpptrajStdio.h"
#include "ParmFile.h" // For writing temporary top

Energy_Sander::Energy_Sander() :
  debug_(0),
  top_pindex_(-1),
  specified_cut_(false),
  specified_igb_(false),
  specified_ntb_(false)
{}

Energy_Sander::~Energy_Sander() { if (is_setup()) sander_cleanup(); }

const char* Energy_Sander::Estring_[] = {
  "Total", "VDW", "Elec", "GB", "Bond", "Angle", "Dihedral", "VDW14", "Elec14",
  "Constraint", "Polar", "Hbond", "Surf", "Cavity", "SCF", "Dispersion", "DVDL", "Angle_UB",
  "Improper", "CMap", "EMap", "LES", "NOE", "PB", "RISM", "CT", "aMD_Boost", 0
};

// Energy_Sander::Easpect()
std::string Energy_Sander::Easpect(Etype typeIn) {
  if (typeIn == N_ENERGYTYPES) return std::string();
  std::locale loc;
  std::string aspect( Estring_[typeIn] );
  for (std::string::iterator it = aspect.begin(); it != aspect.end(); ++it)
    *it = std::tolower( *it, loc );
  return aspect;
}

// Energy_Sander::Energy()
double Energy_Sander::Energy(Etype typeIn) const {
  switch (typeIn) {
    case TOTAL:      return energy_.tot; // Total PE
    case VDW:        return energy_.vdw; // van der Waals
    case ELEC:       return energy_.elec; // Electrostatic
    case GB:         return energy_.gb; // Generalized Born
    case BOND:       return energy_.bond; // Bond
    case ANGLE:      return energy_.angle; // Angle
    case DIHEDRAL:   return energy_.dihedral; // Torsion
    case VDW14:      return energy_.vdw_14; // 1-4 non-bonded
    case ELEC14:     return energy_.elec_14; // 1-4 electrostatic
    case CONSTRAINT: return energy_.constraint; // Constraint
    case POLAR:      return energy_.polar; // Polarization
    case HBOND:      return energy_.hbond; // LJ 10-12 (hydrogen bond)
    case SURF: 
    case CAVITY:     return energy_.surf; // Surface area or cavity energy
    case SCF:        return energy_.scf; // QM/MM SCF energy
    case DISP:       return energy_.disp; // implicit solvation dispersion energy
    case DVDL:       return energy_.dvdl; // TI charging free energy
    case ANGLE_UB:   return energy_.angle_ub; // CHARMM Urey-Bradley
    case IMP:        return energy_.imp; // CHARMM improper
    case CMAP:       return energy_.cmap; // CHARMM CMAP
    case EMAP:       return energy_.emap; // EMAP constraint energy
    case LES:        return energy_.les; // LES
    case NOE:        return energy_.noe; // NOE
    case PB:         return energy_.pb; // Poisson-Boltzmann
    case RISM:       return energy_.rism; // 3D-RISM total solvation free energy
    case CT:         return energy_.ct; // Charge transfer
    case AMD_BOOST:  return energy_.amd_boost; // AMD boost energy
    case N_ENERGYTYPES: break;
  }
  return 0.0; 
}

// Energy_Sander::Eptr()
const double* Energy_Sander::Eptr(Etype typeIn) const {
  switch (typeIn) {
    case TOTAL:      return &energy_.tot;
    case VDW:        return &energy_.vdw;
    case ELEC:       return &energy_.elec;
    case GB:         return &energy_.gb;
    case BOND:       return &energy_.bond;
    case ANGLE:      return &energy_.angle;
    case DIHEDRAL:   return &energy_.dihedral;
    case VDW14:      return &energy_.vdw_14;
    case ELEC14:     return &energy_.elec_14;
    case CONSTRAINT: return &energy_.constraint;
    case POLAR:      return &energy_.polar;
    case HBOND:      return &energy_.hbond;
    case SURF:
    case CAVITY:     return &energy_.surf;
    case SCF:        return &energy_.scf;
    case DISP:       return &energy_.disp;
    case DVDL:       return &energy_.dvdl;
    case ANGLE_UB:   return &energy_.angle_ub;
    case IMP:        return &energy_.imp;
    case CMAP:       return &energy_.cmap;
    case EMAP:       return &energy_.emap;
    case LES:        return &energy_.les;
    case NOE:        return &energy_.noe;
    case PB:         return &energy_.pb;
    case RISM:       return &energy_.rism;
    case CT:         return &energy_.ct;
    case AMD_BOOST:  return &energy_.amd_boost;
    case N_ENERGYTYPES: break;
  }
  return 0;
}

/** Set default input. */
void Energy_Sander::SetDefaultInput() {
  input_.extdiel = 78.5;
  input_.intdiel = 1.0;
  input_.rgbmax = 25.0;
  input_.saltcon = 0.0;
  input_.cut = 8.0;
  input_.dielc = 1.0;
  input_.rdt = 0.0;
  input_.fswitch = -1.0;
  input_.restraint_wt = 0.0;

  input_.igb = 0;
  input_.alpb = 0;
  input_.gbsa = 0;
  input_.lj1264 = -1;
  input_.ipb = 0;
  input_.inp = 2;
  input_.vdwmeth = 1;
  input_.ew_type = 0;
  input_.ntb = 0;
  input_.ifqnt = 0;
  input_.jfastw = 0;
  input_.ntf = 2;
  input_.ntc = 2;
  input_.ntr = 0;

  input_.restraintmask[0] = '\0';
}
 
/** Check and set input for Sander.*/
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
  input_.restraint_wt = argIn.getKeyDouble("restraint_wt", input_.restraint_wt);

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
  // TODO QM/MM not yet supported
  if (input_.ifqnt != 0) {
    mprinterr("Error: 'ifqnt' > 0 not yet supported.\n");
    return 1;
  }
  input_.jfastw = argIn.getKeyInt("jfastw", input_.jfastw);
  input_.ntf = argIn.getKeyInt("ntf", input_.ntf);
  input_.ntc = argIn.getKeyInt("ntc", input_.ntc);
  input_.ntr = argIn.getKeyInt("ntr", input_.ntr);
  // FIXME Currently no way with API to pass in reference coords
  if (input_.ntr > 0) {
    mprinterr("Error: ntr > 0 not yet supported.\n");
    return 1;
  }

  std::string restraintmask = argIn.GetStringKey("restraintmask");
  if (restraintmask.size() > 255) {
    mprinterr("Error: 'restraintmask' is too big.\n");
    return 1;
  }
  restraintmask.copy( input_.restraintmask, restraintmask.size(), 0 );
  // Temporary parm file name
  top_filename_ = argIn.GetStringKey("parmname", "CpptrajEsander.parm7");
  return 0;
}

/** Currently the SANDER API requires a topology file to work. However, 
  * it is not guaranteed that there will be a valid Amber topology file
  * present (for example, after stripping a Topology et cetera). Write
  * a temporary Topology file that can be used by the API.
  * FIXME Should probably use mkstemp etc.
  */
int Energy_Sander::WriteTop( Topology const& topIn ) {
  ParmFile pfile;
  return pfile.WriteTopology( topIn, top_filename_, ArgList(), ParmFile::AMBERPARM, debug_ );
}

#ifdef MPI
int Energy_Sander::Initialize(Topology const& topIn, Frame& fIn, Parallel::Comm const& commIn)
{
  int err = 0;
  if (commIn.Master()) { // TODO MasterBcast?
    // Master writes temporary top file.
    err = WriteTop( topIn );
    if (commIn.CheckError( err )) return 1;
    err = CommonInit( topIn, fIn );
    // Master sends reference frame.
    for (int rank = 1; rank < commIn.Size(); rank++)
      fIn.SendFrame( rank, commIn ); // FIXME make SendFrame const
  } else {
    // Check if master was able to write temporary top file.
    if (commIn.CheckError( err )) return 1;
    // Receive reference frame
    fIn.RecvFrame( 0, commIn );
    err = CommonInit( topIn, fIn );
  }
  // Check that everyone initialized.
  if (commIn.CheckError( err )) return 1;
  return 0;
}
#else
int Energy_Sander::Initialize(Topology const& topIn, Frame& fIn) {
  if (WriteTop( topIn )) return 1;
  return CommonInit( topIn, fIn );
}
#endif

/** Initialize or re-initialize for given Topology with given Frame. Also
  * determine which energy terms will be active.
  */
int Energy_Sander::CommonInit(Topology const& topIn, Frame& fIn) { // TODO const Frame?
  if (fIn.Natom() != topIn.Natom()) {
    mprinterr("Internal Error: Energy_Sander: Input top # atoms (%i) != frame (%i)\n",
              topIn.Natom());
    return 3;
  }
  if (top_filename_.empty()) {
    mprinterr("Internal Error: Energy_Sander: No file name set.\n");
    return 2; // SANITY CHECK
  }
  if (is_setup()) sander_cleanup();
  if (debug_ > 0)
    mprintf("DEBUG: Topology filename= '%s'\n", top_filename_.full());
  // Set some input options if not already specified.
  if (!specified_ntb_) {
    if (fIn.BoxCrd().HasBox())
      input_.ntb = 1;
    else
      input_.ntb = 0;
    mprintf("Warning: 'ntb' not specified; setting to %i based on box type '%s'\n",
            input_.ntb, fIn.BoxCrd().TypeName());
  }
  if (!specified_cut_) {
    if (input_.ntb == 0)
      input_.cut = 9999.0;
    else
      input_.cut = 8.0;
    mprintf("Warning: 'cut' not specified; setting to %.1f based on 'ntb'\n", input_.cut);
  }
  if (!specified_igb_) {
    if (input_.ntb > 0)
      input_.igb = 0; // No warning for ntb > 1 and no igb set
    else if (input_.ntb == 0) {
      input_.igb = 1;
      mprintf("Warning: 'igb' not specified and ntb==0; setting to %i\n", input_.igb);
    }
  }

  mprintf("    SANDER INPUT OPTIONS:\n"
          "      extdiel= %6.2f  intdiel= %6.2f  rgbmax=  %6.2f  saltcon= %6.2f  cut=  %8.2f\n"
          "      dielc=   %6.2f  rdt=     %6.2f  fswitch= %6.2f  igb=     %6i  alpb= %8i\n"
          "      gbsa=    %6i  ntb=     %6i  lj1264=  %6i  ipb=     %6i  inp=  %8i\n"
          "      vdwmeth= %6i  ew_type= %6i  ifqnt=   %6i  jfastw=  %6i  ntf=  %8i\n"
          "      ntc=     %6i\n", input_.extdiel, input_.intdiel,
          input_.rgbmax, input_.saltcon, input_.cut, input_.dielc,
          input_.rdt, input_.fswitch, input_.igb, input_.alpb, input_.gbsa, 
          input_.ntb, input_.lj1264, input_.ipb, input_.inp, input_.vdwmeth,
          input_.ew_type, input_.ifqnt, input_.jfastw, input_.ntf, input_.ntc);

  forces_.resize( fIn.size(), 0.0 );
  // Determine which energy terms will be active based on input
  isActive_.assign( (int)N_ENERGYTYPES, false );
  // Always active
  isActive_[TOTAL] = true;
  isActive_[VDW] = true;
  isActive_[ELEC] = true;
  isActive_[BOND] = true;
  isActive_[ANGLE] = true;
  isActive_[DIHEDRAL] = true;
  isActive_[VDW14] = true;
  isActive_[ELEC14] = true;
  // Implicit solvent
  if (input_.igb > 0) {
    if (input_.igb == 10 || input_.ipb != 0) {
      isActive_[PB] = true;
      isActive_[CAVITY] = true;
      isActive_[DISP] = true;
    } else {
      isActive_[GB] = true;
      if (input_.gbsa > 0) isActive_[SURF] = true;
    }
  }
  // Restraints?
  if (input_.ntr > 0) isActive_[CONSTRAINT] = true;
  // Polar?
  if (topIn.Ipol() != 0) isActive_[POLAR] = true;
  // Are hbond terms present?
  if (!topIn.Nonbond().HBarray().empty()) isActive_[HBOND] = true;
  // Are CHARMM terms present?
  if (topIn.Chamber().HasChamber()) {
    isActive_[ANGLE_UB] = true;
    isActive_[IMP] = true;
    isActive_[CMAP] = true;
  }
  if (debug_ > 0) {
    mprintf("DEBUG: The following terms are active:");
    for (unsigned int i = 0; i != isActive_.size(); i++)
      if (isActive_[i]) mprintf(" %s", Estring_[i]);
    mprintf("\n");
  }
  return sander_setup_mm(top_filename_.full(), fIn.xAddress(), fIn.bAddress(), &input_);
}

// Energy_Sander::CalcEnergy()
int Energy_Sander::CalcEnergy(Frame& fIn) {
  if (!is_setup()) return 1;

  set_positions( fIn.xAddress() );
  set_box( fIn.BoxCrd().BoxX(),  fIn.BoxCrd().BoxY(), fIn.BoxCrd().BoxZ(),
           fIn.BoxCrd().Alpha(), fIn.BoxCrd().Beta(), fIn.BoxCrd().Gamma() );
  energy_forces( &energy_, &(forces_[0]) );
  return 0;
};

// Energy_Sander::CalcEnergyForces()
int Energy_Sander::CalcEnergyForces(Frame& fIn) {
  if (!is_setup()) return 1;

  set_positions( fIn.xAddress() );
  set_box( fIn.BoxCrd().BoxX(),  fIn.BoxCrd().BoxY(), fIn.BoxCrd().BoxZ(),
           fIn.BoxCrd().Alpha(), fIn.BoxCrd().Beta(), fIn.BoxCrd().Gamma() );
  energy_forces( &energy_, fIn.fAddress() );
  return 0;
}
#endif
