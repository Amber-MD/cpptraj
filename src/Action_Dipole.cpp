#include <cmath> // sqrt
#include "Action_Dipole.h"
#include "CpptrajStdio.h"
#include "Constants.h" // SMALL

// CONSTRUCTOR
Action_Dipole::Action_Dipole() :
  max_(0),
  CurrentParm_(0)
{}

void Action_Dipole::Help() {
  mprintf("\t<filename> %s\n", Grid::HelpText);
  mprintf("\t<mask1> {origin | box} [max <max_percent>]\n");
}

// Action_Dipole::init()
Action::RetType Action_Dipole::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  // Get output filename
  filename_ = actionArgs.GetStringNext();
  if (filename_.empty()) {
    mprinterr("Error: Dipole: no filename specified.\n");
    return Action::ERR;
  }
  // 'negative' means something different here than for other grid actions,
  // so get it here. Done this way to be consistent with PTRAJ behavior.
  if (actionArgs.hasKey("negative"))
    max_ = 1;
  else
    max_ = actionArgs.getKeyDouble("max", 0);
  // Get grid options
  if (grid_.GridInit( "Dipole", actionArgs ))
    return Action::ERR;
  // Setup dipole x, y, and z grids
  dipolex_.resize( grid_.GridSize(), 0 );
  dipoley_.resize( grid_.GridSize(), 0 );
  dipolez_.resize( grid_.GridSize(), 0 );

  // Get mask
  std::string maskexpr = actionArgs.GetMaskNext();
  if (maskexpr.empty()) {
    mprinterr("Error: Dipole: No mask specified.\n");
    return Action::ERR;
  }
  mask_.SetMaskString(maskexpr);

  // Info
  grid_.GridInfo();
  mprintf("\tGrid will be printed to file %s\n",filename_.c_str());
  mprintf("\tMask expression: [%s]\n",mask_.MaskString());
  if (max_ > 0)
    mprintf("\tOnly keeping density >= to %.0lf%% of the maximum density\n", max_);

  return Action::OK;
}

// Action_Dipole::setup()
Action::RetType Action_Dipole::Setup(Topology* currentParm, Topology** parmAddress) {
  if (currentParm->Nsolvent() < 1) {
    mprinterr("Error: Dipole: no solvent present in %s.\n", currentParm->c_str());
    return Action::ERR;
  }
  // Traverse over solvent molecules to find out the 
  // "largest" solvent molecule; allocate space for this
  // many coordinates.
  int NsolventAtoms = 0;
  for (Topology::mol_iterator Mol = currentParm->MolStart();
                              Mol != currentParm->MolEnd(); ++Mol)
  {
    if ( (*Mol).IsSolvent() ) {
      if ( (*Mol).NumAtoms() > NsolventAtoms )
        NsolventAtoms = (*Mol).NumAtoms();
    }
  }
  //sol_.resize( NsolventAtoms );
  mprintf("\tLargest solvent mol is %i atoms.\n", NsolventAtoms);

  // Setup grid, checks box info.
  if (grid_.GridSetup( *currentParm )) return Action::ERR;

  // Setup mask
  if (currentParm->SetupCharMask( mask_ ))
    return Action::ERR;
  mprintf("\t[%s] %i atoms selected.\n", mask_.MaskString(), mask_.Nselected());
  if (mask_.None()) {
    mprinterr("Error: Dipole: No atoms selected for parm %s\n", currentParm->c_str());
    return Action::ERR;
  }
  CurrentParm_ = currentParm;
  return Action::OK;
}

// Action_Dipole::action()
Action::RetType Action_Dipole::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  double sol[3];
  Vec3 cXYZ;
  double dipolar_vector[3], COM[3];

  // Set up center to origin or box center
  if (grid_.GridMode() == Grid::BOX) 
    cXYZ = currentFrame->BoxCrd().Center();
  else if (grid_.GridMode() == Grid::CENTER)
    cXYZ = currentFrame->VGeometricCenter( grid_.CenterMask() );
  else
    cXYZ.Zero();

  // Traverse over solvent molecules.
  //int i_solvent = 0; // DEBUG
  for (Topology::mol_iterator solvmol = CurrentParm_->MolStart();
                              solvmol != CurrentParm_->MolEnd(); ++solvmol)
  {
    if (!(*solvmol).IsSolvent()) continue;
    //++i_solvent; // DEBUG
    dipolar_vector[0] = 0.0;
    dipolar_vector[1] = 0.0;
    dipolar_vector[2] = 0.0;
    COM[0] = 0.0;
    COM[1] = 0.0;
    COM[2] = 0.0;
    double total_mass = 0;
    // Loop over solvent atoms
    for (int satom = (*solvmol).BeginAtom(); satom < (*solvmol).EndAtom(); ++satom)
    {
      if ( mask_.AtomInCharMask(satom) ) {
        // Get coordinates and shift to origin and then to appropriate spacing
        const double* XYZ = currentFrame->XYZ( satom );
        sol[0] = XYZ[0] + grid_.SX() - cXYZ[0];
        sol[1] = XYZ[1] + grid_.SY() - cXYZ[1];
        sol[2] = XYZ[2] + grid_.SZ() - cXYZ[2];
        // Calculate dipole vector. The center of mass of the solvent is used 
        // as the "origin" for the vector.
        // NOTE: the total charge on the solvent should be neutral for this 
        //       to have any meaning.
        double mass = (*CurrentParm_)[satom].Mass();
        total_mass += mass;
        COM[0] += (mass * sol[0]);
        COM[1] += (mass * sol[1]);
        COM[2] += (mass * sol[2]);

        double charge = (*CurrentParm_)[satom].Charge();
        dipolar_vector[0] += (charge * sol[0]);
        dipolar_vector[1] += (charge * sol[1]);
        dipolar_vector[2] += (charge * sol[2]);
      }
    }
    // If no atoms selected for this solvent molecule, skip.
    if (total_mass < SMALL) continue;

    // Grid COM
    COM[0] /= total_mass;
    COM[1] /= total_mass;
    COM[2] /= total_mass;
    int bin = grid_.BinPoint( COM[0], COM[1], COM[2] );
    //mprintf("CDBG: Solvent %i XYZ %8.3lf %8.3lf %8.3lf\n",i_solvent,COM[0],COM[1],COM[2]);
    //mprintf("CDBG: Bin = %i\n", bin); 

    // Grid dipole if COM was binned
    if (bin != -1) {
      dipolex_[bin] += dipolar_vector[0] ;
      dipoley_[bin] += dipolar_vector[1] ;
      dipolez_[bin] += dipolar_vector[2] ;
    }
  } // END loop over solvent molecules

  return Action::OK;
}

// Action_Dipole::print()
/** Print dipole data in format for Chris Bayly's discern delegate that 
  * comes with Midas/Plus.
  */
void Action_Dipole::Print() {
  double max_density;
  CpptrajFile outfile;

  if (outfile.OpenWrite(filename_)) {
    mprinterr("Error: Dipole: Cannot open output file.\n");
    return;
  }

  // Write header
  outfile.Printf("field 8\nsize 1\nnside 3\nnlayer 1\ndirectional\nvector\ndata\n");

  // Determine max density
  float maxF = 0;
  for (Grid::iterator gval = grid_.begin(); gval != grid_.end(); ++gval) {
    if ( maxF < *gval )
      maxF = *gval;
  }
  mprintf("\tDipole: maximum density is %f\n", maxF);

  if ( max_ > 0) {
    max_density = max_ * (double)maxF / 100.0;
    mprintf("\tWriting density if >= to %lf\n", max_density);
  } else
    max_density = 1.0;

  // Write data
  for (int k = 0; k < grid_.NZ(); ++k) {
    for (int j = 0; j < grid_.NY(); ++j) {
      for (int i = 0; i < grid_.NX(); ++i) {
        double density = grid_.GridVal(i, j, k);
        //mprintf("CDBG: %5i %5i %5i %lf\n",i,j,k,density);
        if ( density >= max_density ) {
          // Print Bin Coords
          outfile.Printf("%8.3f %8.3f %8.3f", grid_.Xbin(i), grid_.Ybin(j), grid_.Zbin(k));
          // Normalize dipoles by density
          int idx = i*grid_.NY()*grid_.NZ() + j*grid_.NZ() + k;
          double dx = dipolex_[idx] / density; 
          double dy = dipoley_[idx] / density; 
          double dz = dipolez_[idx] / density;
          // Write dipole components and length
          outfile.Printf(" %8.3f %8.3f %8.3f", dx, dy, dz);
          outfile.Printf(" %8.3f %8.3f\n", sqrt(dx*dx+dy*dy+dz*dz), density);
        }
      }
    }
  }
  outfile.CloseFile();
} 
