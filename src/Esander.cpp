#ifdef USE_SANDERLIB
#include "Esander.h"
#include "CpptrajStdio.h"

Energy_Sander::~Energy_Sander() { if (is_setup()) sander_cleanup(); }

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
