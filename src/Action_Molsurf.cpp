// Molsurf
#include <cstring> // strcpy, memset
#include "Action_Molsurf.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Molsurf::Molsurf() {
  //fprintf(stderr,"Angle Con\n");
  sasa=NULL;
  atom = NULL;
  probe_rad = 1.4;
  rad_offset = 0;

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
Molsurf::~Molsurf() {
  ClearMemory();
  if (atom!=NULL) delete[] atom;
}

// MolSurf::ClearMemory()
/// Clear mem used by molsurf data structures
void Molsurf::ClearMemory() {
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

// Molsurf::AllocateMemory()
/// Allocate mem used by molsurf internal data structures
int Molsurf::AllocateMemory() {
  int error_status = 0;
  int natm_sel = Mask1.Nselected;
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

// Molsurf::init()
/** Expected call: molsurf [<name>] [<mask1>] [out filename] [probe <probe_rad>]
                           [offset <rad_offset>]
  * Dataset name will be the last arg checked for. Check order is:
  *    1) Keywords
  *    2) Masks
  *    3) Dataset name
  */
int Molsurf::init() {
  char *mask1;
  char *molsurfFile;

  // Get keywords
  molsurfFile = actionArgs.getKeyString("out",NULL);
  probe_rad = actionArgs.getKeyDouble("probe",1.4);
  rad_offset = actionArgs.getKeyDouble("offset",0.0);

  // Get Masks
  mask1 = actionArgs.getNextMask();
  Mask1.SetMaskString(mask1);

  // Dataset to store angles
  sasa = DSL->Add(DOUBLE, actionArgs.getNextString(),"MSURF");
  if (sasa==NULL) return 1;
  // Add dataset to data file list
  DFL->Add(molsurfFile,sasa);

  mprintf("    MOLSURF: [%s] Probe Radius=%.3lf\n",Mask1.MaskString(),probe_rad);
  if (rad_offset>0)
    mprintf("             Radii will be incremented by %.3lf\n",rad_offset);

  return 0;
}

// Molsurf::setup()
/** Set mask up for this parmtop.
  * currentParm is set in Action::Setup
  */
int Molsurf::setup() {
  if ( currentParm->SetupIntegerMask(Mask1, activeReference) ) return 1;
  if (Mask1.None()) {
    mprintf("    Error: Molsurf::setup: Mask contains 0 atoms.\n");
    return 1;
  }

  mprintf("    MOLSURF: Calculating surface area for %i atoms.\n",Mask1.Nselected);
  // NOTE: If Mask is * dont include any solvent?

  // The ATOM structure is how molsurf organizes atomic data. Allocate
  // here and fill in parm info. Coords will be filled in during action. 
  if (atom!=NULL) delete[] atom;
  atom = new ATOM[ Mask1.Nselected ];
  if (atom==NULL) {
    mprinterr("Error: Molsurf::Setup Could not allocate memory for ATOMs.\n");
    return 1;
  }
  // Set up parm info for atoms in mask
  double *Radii = currentParm->GB_radii_ptr();
  if (Radii==NULL) {
    mprinterr("Error: Molsurf::Setup: Molsurf requires radii, but no radii in %s\n",
              currentParm->parmName);
    return 1;
  }
  for (int maskidx = 0; maskidx < Mask1.Nselected; maskidx++) {
    int parmatom = Mask1.Selected[maskidx];
    int nres = currentParm->atomToResidue(parmatom);
    atom[maskidx].anum = parmatom + 1; // anum is for debug output only, atoms start from 1
    strcpy(atom[maskidx].anam,currentParm->AtomName(parmatom));
    atom[maskidx].rnum = nres + 1; // again for debug output only, residues start from 1
    strcpy(atom[maskidx].rnam,currentParm->ResidueName(nres));
    atom[maskidx].pos[0] = 0;
    atom[maskidx].pos[1] = 0;
    atom[maskidx].pos[2] = 0;
    atom[maskidx].q = currentParm->AtomCharge(parmatom);
    atom[maskidx].rad = Radii[parmatom] + rad_offset;
  }

  // De-allocate memory first since # atoms may have changed
  ClearMemory();
  if (AllocateMemory()) return 1;

  if (debug>0) memory_usage();

  return 0;  
}

// Molsurf::action()
int Molsurf::action() {
  double molsurf_sasa;

  // Set up coordinates for atoms in mask
  for (int maskidx = 0; maskidx < Mask1.Nselected; maskidx++) {
    int i3 = Mask1.Selected[maskidx] * 3;
    atom[maskidx].pos[0] = currentFrame->X[i3  ];
    atom[maskidx].pos[1] = currentFrame->X[i3+1];
    atom[maskidx].pos[2] = currentFrame->X[i3+2];
  }

  // NOTE: cusp_edge is the only data structure that requires initialization 
  memset( cusp_edge, 0, NUM_EDGE * Mask1.Nselected * sizeof(CUSP_EDGE));
  molsurf_sasa = molsurf( probe_rad, atom, Mask1.Nselected,
                          upper_neighbors, neighbors,
                          toruslist, probelist, concave_face,
                          saddle_face, convex_face,
                          cone_face, broken_concave_face,
                          concave_cycle, vertexlist,
                          concave_edge_list, convex_edge_list,
                          convex_circle_list, concave_circle_list,
                          cyclelist, low_torus, cusp_edge,
                          cusp_pair); 

  sasa->Add(frameNum, &molsurf_sasa);

  //fprintf(outfile,"%10i %10.4lf\n",frameNum,molsurf_sasa);
  
  return 0;
} 


