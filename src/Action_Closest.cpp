#include <cmath>
#include <algorithm> // sort
#include <cfloat> // DBL_MAX
#ifdef _OPENMP
#  include <omp.h>
#endif
#include "Action_Closest.h"
#include "CpptrajStdio.h"
#include "ParmFile.h"
#ifdef CUDA
#  include <cuda_runtime_api.h>
#  include <cuda.h>

// CUDA kernels
extern void Action_NoImage_Center(double*,double*,double[3],double,int,int,float&);
extern void Action_NoImage_no_Center(double*,double*,double*,double,int,int,int,float&);
#endif

// CONSTRUCTOR
Action_Closest::Action_Closest() :
  outFile_(0),
  framedata_(0),
  moldata_(0),
  distdata_(0),
  atomdata_(0),
  Nclosest_(0),
  closestWaters_(0),
  firstAtom_(false),
  useMaskCenter_(false),
  newParm_(0),
  NsolventMolecules_(0),
  debug_(0)
{} 

void Action_Closest::Help() const {
  mprintf("\t<# to keep> <mask> [noimage] [first | oxygen] [center]\n"
          "\t[closestout <filename> [name <setname>]] [outprefix <parmprefix>]\n"
          "\t[parmout <file>]\n"
          "  Keep only the closest <# to keep> solvent molecules to atoms in <mask>.\n"
          "  Molecules can be marked as solvent with the 'solvent' command.\n"
          "  If 'center' specified use geometric center of atoms in <mask>.\n");
}

// DESTRUCTOR
Action_Closest::~Action_Closest() {
  if (newParm_!=0) delete newParm_;
}

// Action_Closest::Init()
Action::RetType Action_Closest::Init(ArgList& actionArgs, ActionInit& init, int debugIn)
{
# ifdef MPI
  trajComm_ = init.TrajComm();
# endif
  debug_ = debugIn;
  // Get Keywords
  closestWaters_ = actionArgs.getNextInteger(-1);
  if (closestWaters_ < 0) {
    mprinterr("Error: Invalid # solvent molecules to keep (%i).\n",
              closestWaters_);
    return Action::ERR;
  }
  if ( actionArgs.hasKey("oxygen") || actionArgs.hasKey("first") )
    firstAtom_=true;
  useMaskCenter_ = actionArgs.hasKey("center");
  image_.InitImaging( !(actionArgs.hasKey("noimage")) );
  prefix_ = actionArgs.GetStringKey("outprefix");
  parmoutName_ = actionArgs.GetStringKey("parmout");
  // Setup output file and sets if requested.
  // Will keep track of Frame, Mol#, Distance, and first solvent atom
  std::string filename = actionArgs.GetStringKey("closestout");
  if (!filename.empty()) {
    std::string dsetName = actionArgs.GetStringKey("name");
    if (dsetName.empty())
      dsetName = init.DSL().GenerateDefaultName("CLOSEST");
    // Set up datasets
    framedata_ = init.DSL().AddSet(DataSet::INTEGER, MetaData(dsetName, "Frame"));
    moldata_   = init.DSL().AddSet(DataSet::INTEGER, MetaData(dsetName, "Mol"));
    distdata_  = init.DSL().AddSet(DataSet::DOUBLE,  MetaData(dsetName, "Dist"));
    atomdata_  = init.DSL().AddSet(DataSet::INTEGER, MetaData(dsetName, "FirstAtm"));
    if (framedata_==0 || moldata_==0 || distdata_==0 || atomdata_==0) {
      mprinterr("Error: Could not setup data sets for output file %s\n",
                filename.c_str());
      return Action::ERR;
    }
    // Add sets to datafile in list.
    outFile_ = init.DFL().AddDataFile( filename );
    if (outFile_ == 0) {
      mprinterr("Error: Could not set up output file %s\n", filename.c_str());
      return Action::ERR;
    }
    outFile_->AddDataSet(framedata_);
    outFile_->AddDataSet(moldata_);
    outFile_->AddDataSet(distdata_);
    outFile_->AddDataSet(atomdata_);
    outFile_->ProcessArgs("noxcol");
  }

  // Get Masks
  std::string mask1 = actionArgs.GetMaskNext();
  if (mask1.empty()) {
    mprinterr("Error: No mask specified.\n");
    return Action::ERR;
  }
  distanceMask_.SetMaskString(mask1);

  mprintf("    CLOSEST: Finding closest %i solvent molecules to atoms in mask %s\n",
          closestWaters_, distanceMask_.MaskString());
  if (useMaskCenter_)
    mprintf("\tGeometric center of atoms in mask will be used.\n");
  if (!image_.UseImage()) 
    mprintf("\tImaging will be turned off.\n");
  if (firstAtom_)
    mprintf("\tOnly first atom of solvent molecule used for distance calc.\n");
  if (outFile_!=0)
    mprintf("\tClosest molecules will be saved to %s\n",outFile_->DataFilename().base());
  if (!prefix_.empty())
    mprintf("\tStripped topology file will be written with prefix %s\n",
            prefix_.c_str());
  return Action::OK;
}

// Action_Closest::Setup()
/** Like the strip action, closest will modify the current parm keeping info
  * for atoms in mask plus the closestWaters solvent molecules. Set up the
  * vector of MolDist objects, one for every solvent molecule in the original
  * parm file. Atom masks for each solvent molecule will be set up.
  */
Action::RetType Action_Closest::Setup(ActionSetup& setup) {
  // If there are no solvent molecules this action is not valid.
  if (setup.Top().Nsolvent()==0) {
    mprintf("Warning: Parm %s does not contain solvent.\n",setup.Top().c_str());
    return Action::SKIP;
  }
  // If # solvent to keep >= solvent in this parm the action is not valid.
  if (closestWaters_ >= setup.Top().Nsolvent()) {
    mprintf("Warning: # solvent to keep (%i) >= # solvent molecules in '%s' (%i)\n",
            closestWaters_, setup.Top().c_str(), setup.Top().Nsolvent());
    return Action::SKIP;
  }
  image_.SetupImaging( setup.CoordInfo().TrajBox().Type() );
  if (image_.ImagingEnabled())
    mprintf("\tDistances will be imaged.\n");
  else
    mprintf("\tImaging off.\n"); 
  // LOOP OVER MOLECULES
  // 1: Check that all solvent molecules contain same # atoms. Solvent 
  //    molecules must be identical for the command to work properly; 
  //    the prmtop strip occurs only once so the solvent params become fixed.
  // 2: Set up a mask for all solvent molecules.
  SolventMols_.clear();
  // NOTE: May not be necessary to init 'solvent'
  MolDist solvent;
  solvent.D = 0.0;
  solvent.mol = 0;
  SolventMols_.resize(setup.Top().Nsolvent(), solvent);
  std::vector<MolDist>::iterator mdist = SolventMols_.begin();
  // 3: Set up the soluteMask for all non-solvent molecules.
  stripMask_.ResetMask();
  int molnum = 1;
  int newnatom = 0;
  int nclosest = 0;
  int NsolventAtoms = -1;
  keptWaterAtomNum_.resize(closestWaters_);
  for (Topology::mol_iterator Mol = setup.Top().MolStart();
                              Mol != setup.Top().MolEnd(); ++Mol)
  {
    if ( !Mol->IsSolvent() ) { // Not solvent, add to solute mask.
      stripMask_.AddAtomRange( Mol->BeginAtom(), Mol->EndAtom() );
      newnatom += Mol->NumAtoms();
    } else {                   // Solvent, check for same # of atoms.
      if (NsolventAtoms == -1)
        NsolventAtoms = Mol->NumAtoms();
      else if ( NsolventAtoms != Mol->NumAtoms() ) {
        mprinterr("Error: Solvent molecules in '%s' are not of uniform size.\n"
                  "Error:   First solvent mol = %i atoms, solvent mol %i = %i atoms.\n",
                  setup.Top().c_str(), NsolventAtoms, molnum, (*Mol).NumAtoms());
        return Action::ERR;
      }
      // mol here is the output molecule number which is why it starts from 1.
      mdist->mol = molnum;
      // Solvent molecule mask
      mdist->mask.AddAtomRange( Mol->BeginAtom(), Mol->EndAtom() );
      // Atoms in the solvent molecule to actually calculate distances to.
      if (firstAtom_) {
        mdist->solventAtoms.assign(1, Mol->BeginAtom() );
      } else {
        mdist->solventAtoms.clear();
        mdist->solventAtoms.reserve( Mol->NumAtoms() );
        for (int svatom = Mol->BeginAtom(); svatom < Mol->EndAtom(); svatom++)
          mdist->solventAtoms.push_back( svatom );
      }
      // For solvent molecules that will be kept, record what the atom number
      // will be in the new stripped parm.
      if (nclosest < closestWaters_) {
        keptWaterAtomNum_[nclosest] = newnatom;
        stripMask_.AddAtomRange( Mol->BeginAtom(), Mol->EndAtom() );
        newnatom += Mol->NumAtoms();
        ++nclosest;
      }
      //SolventMols[solventMol].mask.PrintMaskAtoms("solvent");
      ++mdist;
    }
    ++molnum;
  }

  // Setup distance atom mask
  // NOTE: Should ensure that no solvent atoms are selected!
  if ( setup.Top().SetupIntegerMask(distanceMask_) ) return Action::ERR;
  if (distanceMask_.None()) {
    mprintf("Warning: Distance mask '%s' contains no atoms.\n",
            distanceMask_.MaskString());
    return Action::SKIP;
  }
  distanceMask_.MaskInfo();

  // Check the total number of solvent atoms to be kept.
  NsolventAtoms *= closestWaters_;
  mprintf("\tKeeping %i solvent atoms.\n",NsolventAtoms);
  if (NsolventAtoms < 1) {
    mprintf("Warning: # of solvent atoms to be kept is < 1.\n");
    return Action::SKIP;
  }
  // Store original # of molecules.
  // NOTE: This is stored so that it can be used in the OpenMP section
  //       of action. I dont think iterators are thread-safe.
  NsolventMolecules_ = setup.Top().Nsolvent();
  // Create stripped Parm
  if (newParm_!=0) delete newParm_;
  newParm_ = setup.Top().modifyStateByMask(stripMask_);
  if (newParm_==0) {
    mprinterr("Error: Could not create new topology.\n");
    return Action::ERR;
  }
  setup.SetTopology( newParm_ );
  newParm_->Brief("Closest topology:");
  // Allocate space for new frame
  newFrame_.SetupFrameV( setup.Top().Atoms(), setup.CoordInfo() );

  // If prefix given then output stripped parm
  if (!prefix_.empty()) {
    ParmFile pfile;
    if ( pfile.WritePrefixTopology(*newParm_, prefix_, ParmFile::AMBERPARM, debug_ ) )
      mprinterr("Error: Could not write out 'closest' topology file.\n");
  }
  if (!parmoutName_.empty()) {
    ParmFile pfile;
    if ( pfile.WriteTopology(*newParm_, parmoutName_, ParmFile::AMBERPARM, debug_) )
      mprinterr("Error: Could not write out topology to file %s\n", parmoutName_.c_str());
  }
# ifdef CUDA
  // Allocate space for simple arrays that will be sent to GPU.
  V_distances_.resize( SolventMols_.size() );
  V_atom_coords_.resize( NsolventAtoms * SolventMols_.size() * 3, 0.0 );
  if (!useMaskCenter_)
    U_atom_coords_.resize( distanceMask_.Nselected() * 3, 0.0 );
# endif

  return Action::MODIFY_TOPOLOGY;
}

// Action_Closest::DoAction()
/** Find the minimum distance between atoms in distanceMask and each 
  * solvent Mask.
  */
Action::RetType Action_Closest::DoAction(int frameNum, ActionFrame& frm) {
  double Dist, maxD;
  Matrix_3x3 ucell, recip;
  AtomMask::const_iterator solute_atom;
  Iarray::const_iterator solvent_atom;

  if (image_.ImagingEnabled()) {
    frm.Frm().BoxCrd().ToRecip(ucell, recip);
    // Calculate max possible imaged distance
    maxD = frm.Frm().BoxCrd().BoxX() + frm.Frm().BoxCrd().BoxY() + 
           frm.Frm().BoxCrd().BoxZ();
    maxD *= maxD;
  } else {
    // If not imaging, set max distance to an arbitrarily large number
    maxD = DBL_MAX;
  }
#ifdef CUDA
// -----------------------------------------------------------------------------
/*
  useMaskCenter_ = false;
  cudaEvent_t start_event, stop_event;
  float elapsed_time_seq;
  mprintf("Solute Center : %s\n", v[k] ? "YES" : "NO");
  useMaskCenter_ = v[k];

  cudaEventCreate(&start_event);
  cudaEventCreate(&stop_event);
  cudaEventRecord(start_event, 0);

  // serial section of the code 
  Action_NoImage(frmIn,maxD);

  cudaThreadSynchronize();
  cudaEventRecord(stop_event, 0);
  cudaEventSynchronize(stop_event);
  cudaEventElapsedTime(&elapsed_time_seq,start_event, stop_event );
  mprintf("Done with kernel SEQ Kernel Time: %.2f\n", elapsed_time_seq);
*/
  float elapsed_time_gpu;
  // Copy solvent atom coords to array
  int NAtoms = SolventMols_[0].solventAtoms.size(); // guaranteed to same size due to setup
  for (int sMol = 0; sMol < NsolventMolecules_; ++sMol) {
    for(int sAtom = 0; sAtom < NAtoms; sAtom++) {
      int index =  (sAtom * 3 ) + (sMol * 3 * NAtoms);
      const double *a = frm.Frm().XYZ(SolventMols_[sMol].solventAtoms[sAtom]);
      V_atom_coords_[index + 0] = a[0];
      V_atom_coords_[index + 1] = a[1];
      V_atom_coords_[index + 2] = a[2];
    }
  }
  // TODO handle imaging
  if (useMaskCenter_) {
    Vec3 maskCenter = frm.Frm().VGeometricCenter( distanceMask_ );
    Action_NoImage_Center( &V_atom_coords_[0], &V_distances_[0], maskCenter.Dptr(),
                           maxD, NsolventMolecules_, NAtoms, elapsed_time_gpu );
  } else {
    int NSAtoms = distanceMask_.Nselected();
    for (int nsAtom = 0; nsAtom < NSAtoms; ++nsAtom) {
      const double* a = frm.Frm().XYZ(distanceMask_[nsAtom]);
      U_atom_coords_[nsAtom*3 + 0] = a[0];
      U_atom_coords_[nsAtom*3 + 1] = a[1];
      U_atom_coords_[nsAtom*3 + 2] = a[2];
    }
    Action_NoImage_no_Center( &V_atom_coords_[0], &V_distances_[0], &U_atom_coords_[0],
                              maxD, NsolventMolecules_, NAtoms, NSAtoms, elapsed_time_gpu );
  }
  // Copy distances back into SolventMols_
  for (int sMol = 0; sMol < NsolventMolecules_; sMol++)
    SolventMols_[sMol].D = V_distances_[sMol];
# ifdef DEBUG_CUDA
  mprintf("CUDA Time: = %0.2f\n", elapsed_time_gpu);
# endif

#else /* Not CUDA */
// -----------------------------------------------------------------------------
  int solventMol; 

  // Loop over all solvent molecules in original frame
  if (useMaskCenter_) {
    Vec3 maskCenter = frm.Frm().VGeometricCenter( distanceMask_ );
#ifdef _OPENMP
#pragma omp parallel private(solventMol,Dist,solvent_atom)
{
#   pragma omp for
#endif
    for (solventMol=0; solventMol < NsolventMolecules_; solventMol++) {
      SolventMols_[solventMol].D = maxD;
      for (solvent_atom = SolventMols_[solventMol].solventAtoms.begin();
           solvent_atom != SolventMols_[solventMol].solventAtoms.end(); ++solvent_atom)
      {
        Dist = DIST2( maskCenter.Dptr(),
                      frm.Frm().XYZ(*solvent_atom), image_.ImageType(),
                      frm.Frm().BoxCrd(), ucell, recip);
        if (Dist < SolventMols_[solventMol].D) 
          SolventMols_[solventMol].D = Dist;
      }
    }
#ifdef _OPENMP
}
#endif
  } else {
#ifdef _OPENMP
#pragma omp parallel private(solventMol,solute_atom,Dist,solvent_atom)
{
    //mprintf("OPENMP: %i threads\n",omp_get_num_threads());
#   pragma omp for
#endif
    for (solventMol=0; solventMol < NsolventMolecules_; solventMol++) {
      //mprintf("[%i] Calculating distance for molecule %i\n",omp_get_thread_num(),solventMol);
      // Set the initial minimum distance for this solvent mol to be the
      // max possible distance.
      SolventMols_[solventMol].D = maxD;
      // Calculate distance between each atom in distanceMask and atoms in solvent Mask
      for (solvent_atom = SolventMols_[solventMol].solventAtoms.begin();
           solvent_atom != SolventMols_[solventMol].solventAtoms.end(); ++solvent_atom)
      {
        for (solute_atom = distanceMask_.begin(); 
             solute_atom != distanceMask_.end(); ++solute_atom)
        {
          Dist = DIST2(frm.Frm().XYZ(*solute_atom),
                       frm.Frm().XYZ(*solvent_atom), image_.ImageType(), 
                       frm.Frm().BoxCrd(), ucell, recip);
          if (Dist < SolventMols_[solventMol].D) 
            SolventMols_[solventMol].D = Dist;
          //mprintf("DBG: SolvMol %i, soluteAtom %i, solventAtom %i, D= %f, minD= %f\n",
          //        solventMol, *solute_atom, *solvent_atom, Dist, sqrt(SolventMols_[solventMol].D));
        }
      }
      // DEBUG - Print distances
      //mprintf("DEBUG:\tMol %8i minD= %lf\n",solventMol, SolventMols[solventMol].D);
    } // END for loop over solventMol
#ifdef _OPENMP
} // END pragma omp parallel
#endif
  }
  // DEBUG
  //mprintf("Closest: End parallel loop for %i, got %i Distances.\n",frameNum,(int)SolventMols.size());
  // DEBUG
#endif /* CUDA */
// -----------------------------------------------------------------------------

  // Sort distances
  std::sort( SolventMols_.begin(), SolventMols_.end(), moldist_cmp() );
  // Add first closestWaters solvent atoms to stripMask
  std::vector<MolDist>::iterator solventend = SolventMols_.begin() + closestWaters_;
  Iarray::const_iterator katom = keptWaterAtomNum_.begin();
  for ( std::vector<MolDist>::const_iterator solvent = SolventMols_.begin();
                                             solvent != solventend;
                                           ++solvent, ++katom ) 
  {
    //mprintf("DEBUG:\tmol %i ",(*solvent).mol);
    //solvent->mask.PrintMaskAtoms("Mask");
    stripMask_.AddMaskAtPosition( solvent->mask, *katom );
    // Record which water molecules are closest if requested
    if (outFile_!=0) {
      int fnum = frm.TrajoutNum() + 1;
      framedata_->Add(Nclosest_, &fnum);
      moldata_->Add(Nclosest_, &(solvent->mol));
      Dist = sqrt( solvent->D );
      distdata_->Add(Nclosest_, &Dist);
      solvent_atom = solvent->mask.begin();
      int solvent_first_atom = *solvent_atom + 1; 
      atomdata_->Add(Nclosest_, &solvent_first_atom);
      ++Nclosest_;
    }
    // DEBUG - print first closestWaters distances
    //mprintf("DEBUG: Mol %i   D2= %lf   Atom0= %i\n",(*it).mol, (*it).D, (*it).mask->Selected[0]);
  }

  // Modify and set frame
  //mprintf("DEBUG:\t");
  //stripMask.PrintMaskAtoms("action_stripMask");
  newFrame_.SetFrame(frm.Frm(), stripMask_);
  frm.SetFrame( &newFrame_ );

  return Action::MODIFY_COORDS;
}

#ifdef MPI
/** Since datasets are actually # frames * closestWaters_ in length, need to
  * sync here.
  */
int Action_Closest::SyncAction() {
  if (outFile_ == 0) return 0;
  // Get total number of closest entries.
  std::vector<int> rank_frames( trajComm_.Size() );
  trajComm_.GatherMaster( &Nclosest_, 1, MPI_INT, &(rank_frames[0]) );
  for (int rank = 1; rank < trajComm_.Size(); rank++)
    Nclosest_ += rank_frames[ rank ];
  framedata_->Sync( Nclosest_, rank_frames, trajComm_ );
  framedata_->SetNeedsSync( false );
  moldata_->Sync( Nclosest_, rank_frames, trajComm_ );
  moldata_->SetNeedsSync( false );
  distdata_->Sync( Nclosest_, rank_frames, trajComm_ );
  distdata_->SetNeedsSync( false );
  atomdata_->Sync( Nclosest_, rank_frames, trajComm_ );
  atomdata_->SetNeedsSync( false );
  return 0;
}
#endif
