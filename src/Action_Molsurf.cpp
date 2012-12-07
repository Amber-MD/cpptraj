// Action_Molsurf
#include <cstring> // strcpy, memset
#include "Action_Molsurf.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_Molsurf::Action_Molsurf() {
  //fprintf(stderr,"Angle Con\n");
  sasa_=NULL;
  atom_ = NULL;
  probe_rad_ = 1.4;
  rad_offset_ = 0;

  upper_neighbors=NULL;
  neighbors=NULL;
  toruslist=NULL;
  probelist=NULL;

  concave_face=NULL;
  saddle_face=NULL;
  convex_face=NULL;
  cone_face=NULL;
  broken_concave_face=NULL;
  concave_cycle=NULL;

  vertexlist=NULL;
  concave_edge_list=NULL;
  convex_edge_list=NULL;
  convex_circle_list=NULL;
  concave_circle_list=NULL;

  cyclelist=NULL;
  low_torus=NULL;
  cusp_edge=NULL;
  cusp_pair=NULL;
} 

// DESTRUCTOR
Action_Molsurf::~Action_Molsurf() {
  ClearMemory();
  if (atom_!=NULL) delete[] atom_;
}

void Action_Molsurf::Help() {
  mprintf("molsurf [<name>] [<mask1>] [out filename] [probe <probe_rad>] [offset <rad_offset>]\n");
  mprintf("\tCalculate Connolly surface area of atoms in <mask1>\n");
}

// MolSurf::ClearMemory()
/// Clear mem used by molsurf data structures
void Action_Molsurf::ClearMemory() {
  if (upper_neighbors!=NULL) delete[] upper_neighbors;
  if (neighbors!=NULL) delete[] neighbors;
  if (probelist!=NULL) delete[] probelist;
  if (toruslist!=NULL) delete[] toruslist;
  if (convex_circle_list!=NULL) delete[] convex_circle_list;
  if (concave_circle_list!=NULL) delete[] concave_circle_list;
  if (concave_face!=NULL) delete[] concave_face;
  if (convex_face!=NULL) delete[] convex_face;
  if (saddle_face!=NULL) delete[] saddle_face;
  if (cone_face!=NULL) delete[] cone_face;
  if (broken_concave_face!=NULL) delete[] broken_concave_face;
  if (concave_cycle!=NULL) delete[] concave_cycle;
  if (cyclelist!=NULL) delete[] cyclelist;
  if (vertexlist!=NULL) delete[] vertexlist;
  if (concave_edge_list!=NULL) delete[] concave_edge_list;
  if (convex_edge_list!=NULL) delete[] convex_edge_list;
  if (low_torus!=NULL) delete[] low_torus;
  if (cusp_edge!=NULL) delete[] cusp_edge;
  if (cusp_pair!=NULL) delete[] cusp_pair;
}

// Action_Molsurf::AllocateMemory()
/// Allocate mem used by molsurf internal data structures
int Action_Molsurf::AllocateMemory() {
  int error_status = 0;
  int natm_sel = Mask1_.Nselected();
  // ---------- Allocate Memory For molsurf routines --------------------
  upper_neighbors = new NEIGHBOR_TORUS[ NUM_NEIGHBOR * natm_sel ];
  neighbors = new NEIGHBOR [ NUM_NEIGHBOR * natm_sel ]; 
  probelist = new PROBE[ NUM_PROBE * natm_sel ];
  toruslist = new TORUS[ NUM_TORUS * natm_sel ];
  convex_circle_list = new CIRCLE[ NUM_CIRCLE * natm_sel ];
  concave_circle_list = new CIRCLE[ NUM_CIRCLE * natm_sel ];
  concave_face = new CONCAVE_FACE[ NUM_FACE * natm_sel ];
  convex_face = new CONVEX_FACE[ NUM_FACE * natm_sel ];
  saddle_face = new SADDLE_FACE[ NUM_FACE * natm_sel ];
  cone_face = new CONE_FACE[ NUM_FACE * natm_sel ];  
  broken_concave_face = new BROKEN_CONCAVE_FACE[ NUM_FACE * natm_sel ];
  concave_cycle = new CONCAVE_CYCLE[ NUM_CYCLE * natm_sel ];
  cyclelist = new CYCLE[ NUM_CYCLE * natm_sel ];
  vertexlist = new VERTEX[ NUM_VERTEX * natm_sel ];
  concave_edge_list = new EDGE[ NUM_EDGE * natm_sel ];
  convex_edge_list = new EDGE[ NUM_EDGE * natm_sel ];
  low_torus = new LOW_TORUS[ NUM_TORUS * natm_sel ];
  // NOTE: cusp_edge must be initialized. Handled in action before 1st call
  cusp_edge = new CUSP_EDGE[ NUM_EDGE * natm_sel ];
  cusp_pair = new CUSP_PAIR[ NUM_CUSP * natm_sel ]; 
  // ------------------------------------------------------------
  return error_status;
}

// Action_Molsurf::init()
Action::RetType Action_Molsurf::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  // Get keywords
  DataFile* outfile = DFL->AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  probe_rad_ = actionArgs.getKeyDouble("probe",1.4);
  rad_offset_ = actionArgs.getKeyDouble("offset",0.0);

  // Get Masks
  Mask1_.SetMaskString( actionArgs.GetMaskNext() );

  // Dataset to store angles
  sasa_ = DSL->AddSet(DataSet::DOUBLE, actionArgs.GetStringNext(),"MSURF");
  if (sasa_==0) return Action::ERR;
  // Add dataset to data file list
  if (outfile != 0) outfile->AddSet(sasa_);

  mprintf("    MOLSURF: [%s] Probe Radius=%.3lf\n",Mask1_.MaskString(),probe_rad_);
  if (rad_offset_>0)
    mprintf("             Radii will be incremented by %.3lf\n",rad_offset_);

  return Action::OK;
}

// Action_Molsurf::setup()
/** Set mask up for this parmtop. Allocate the ATOM structure array used 
  * by the molsurf C routines and set everything but the atom coords.
  */
Action::RetType Action_Molsurf::Setup(Topology* currentParm, Topology** parmAddress) {
  if ( currentParm->SetupIntegerMask(Mask1_) ) return Action::ERR;
  if (Mask1_.None()) {
    mprintf("    Error: Molsurf::setup: Mask contains 0 atoms.\n");
    return Action::ERR;
  }

  mprintf("    MOLSURF: Calculating surface area for %i atoms.\n",Mask1_.Nselected());
  // NOTE: If Mask is * dont include any solvent?

  // The ATOM structure is how molsurf organizes atomic data. Allocate
  // here and fill in parm info. Coords will be filled in during action. 
  if (atom_!=NULL) delete[] atom_;
  atom_ = new ATOM[ Mask1_.Nselected() ];
  if (atom_==NULL) {
    mprinterr("Error: Molsurf::Setup Could not allocate memory for ATOMs.\n");
    return Action::ERR;
  }
  // Set up parm info for atoms in mask
  if ( (*currentParm)[0].Radius() == 0 ) {
    mprinterr("Error: Molsurf::Setup: Molsurf requires radii, but no radii in %s\n",
              currentParm->c_str());
    return Action::ERR;
  }
  ATOM *atm_ptr = atom_;
  for (AtomMask::const_iterator parmatom = Mask1_.begin();
                                parmatom != Mask1_.end();
                                parmatom++)
  {
    atm_ptr->anum = *parmatom + 1; // anum is for debug output only, atoms start from 1
    const Atom patom = (*currentParm)[*parmatom];
    int nres = patom.ResNum();
    atm_ptr->rnum = nres+1; // for debug output only, residues start from 1
    patom.Name().ToBuffer( atm_ptr->anam );
    strcpy(atm_ptr->rnam, currentParm->Res(nres).c_str());
    atm_ptr->pos[0] = 0;
    atm_ptr->pos[1] = 0;
    atm_ptr->pos[2] = 0;
    atm_ptr->q = patom.Charge();
    atm_ptr->rad = patom.Radius() + rad_offset_;
    ++atm_ptr;
  }

  // De-allocate memory first since # atoms may have changed
  ClearMemory();
  if (AllocateMemory()) return Action::ERR;

  //if (debug>0) memory_usage();

  return Action::OK;  
}

// Action_Molsurf::action()
Action::RetType Action_Molsurf::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  // Set up coordinates for atoms in mask
  ATOM *atm_ptr = atom_;
  for (AtomMask::const_iterator maskatom = Mask1_.begin(); maskatom != Mask1_.end(); ++maskatom)
  {
    const double* XYZ = currentFrame->XYZ( *maskatom );
    memcpy(atm_ptr->pos, XYZ, 3*sizeof(double));
    ++atm_ptr;
  }

  // NOTE: cusp_edge is the only data structure that requires initialization 
  memset( cusp_edge, 0, NUM_EDGE * Mask1_.Nselected() * sizeof(CUSP_EDGE));
  double molsurf_sasa = molsurf( probe_rad_, atom_, Mask1_.Nselected(),
                                 upper_neighbors, neighbors,
                                 toruslist, probelist, concave_face,
                                 saddle_face, convex_face,
                                 cone_face, broken_concave_face,
                                 concave_cycle, vertexlist,
                                 concave_edge_list, convex_edge_list,
                                 convex_circle_list, concave_circle_list,
                                 cyclelist, low_torus, cusp_edge,
                                 cusp_pair); 

  sasa_->Add(frameNum, &molsurf_sasa);

  //fprintf(outfile,"%10i %10.4lf\n",frameNum,molsurf_sasa);
  
  return Action::OK;
} 
