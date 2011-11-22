// Molsurf
#include <cstdlib> // Using malloc since interfacing with C code
#include <cstring> 
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
  if (atom!=NULL) free(atom);
}

// MolSurf::ClearMemory()
/// Clear mem used by molsurf data structures
void Molsurf::ClearMemory() {
  if (upper_neighbors!=NULL) free (upper_neighbors);
  if (neighbors!=NULL) free (neighbors);
  if (probelist!=NULL) free (probelist);
  if (toruslist!=NULL) free (toruslist);
  if (convex_circle_list!=NULL) free (convex_circle_list);
  if (concave_circle_list!=NULL) free (concave_circle_list);
  if (concave_face!=NULL) free (concave_face);
  if (convex_face!=NULL) free (convex_face);
  if (saddle_face!=NULL) free (saddle_face);
  if (cone_face!=NULL) free (cone_face);
  if (broken_concave_face!=NULL) free (broken_concave_face);
  if (concave_cycle!=NULL) free (concave_cycle);
  if (cyclelist!=NULL) free (cyclelist);
  if (vertexlist!=NULL) free (vertexlist);
  if (concave_edge_list!=NULL) free (concave_edge_list);
  if (convex_edge_list!=NULL) free (convex_edge_list);
  if (low_torus!=NULL) free (low_torus);
  if (cusp_edge!=NULL) free (cusp_edge);
  if (cusp_pair!=NULL) free (cusp_pair);
}

// Molsurf::AllocateMemory()
/// Allocate mem used by molsurf internal data structures
int Molsurf::AllocateMemory() {
  int error_status = 0;
  int natm_sel = Mask1.Nselected;
  // ---------- Allocate Memory For molsurf routines --------------------
  if ((upper_neighbors = (NEIGHBOR_TORUS *) malloc (
                NUM_NEIGHBOR * natm_sel * sizeof (NEIGHBOR_TORUS))) == NULL) {
        mprinterr("Unable to allocate space for upper_neighbors\n");
        error_status++;
  }
  if ((neighbors = (NEIGHBOR *) malloc (
                  NUM_NEIGHBOR * natm_sel * sizeof (NEIGHBOR))) == NULL) {
        mprinterr("Unable to allocate space for neighbors\n");
        error_status++;
  }
  if ((probelist = (PROBE *) malloc (
                NUM_PROBE * natm_sel * sizeof (PROBE))) == NULL) {
        mprinterr("Unable to allocate space for probelist\n");
        error_status++;
  }
  if ((toruslist = (TORUS *) malloc (
                NUM_TORUS * natm_sel * sizeof (TORUS))) == NULL) {
        mprinterr("Unable to allocate space for toruslist\n");
        error_status++;
  }
  if ((convex_circle_list = (CIRCLE *) malloc (
                NUM_CIRCLE * natm_sel * sizeof (CIRCLE))) == NULL) {
        mprinterr("Unable to allocate space for convex_circle_list\n");
        error_status++;
  }
  if ((concave_circle_list = (CIRCLE *) malloc (
                NUM_CIRCLE * natm_sel * sizeof (CIRCLE))) == NULL) {
        mprinterr("Unable to allocate space for concave_circle_list\n");
        error_status++;
  }
  if ((concave_face = (CONCAVE_FACE *) malloc (
                NUM_FACE * natm_sel * sizeof (CONCAVE_FACE))) == NULL) {
        mprinterr("Unable to allocate space for concave_face\n");
        error_status++;
  }
  if ((convex_face = (CONVEX_FACE *) malloc (
           NUM_FACE * natm_sel * sizeof (CONVEX_FACE))) == NULL) {
        mprinterr("Unable to allocate space for convex_face\n");
        error_status++;
  }
  if ((saddle_face = (SADDLE_FACE *) malloc (
           NUM_FACE * natm_sel * sizeof (SADDLE_FACE))) == NULL) {
        mprinterr("Unable to allocate space for saddle_face\n");
        error_status++;
  }
  if ((cone_face = (CONE_FACE *) malloc (
                 NUM_FACE * natm_sel * sizeof (CONE_FACE))) == NULL) {
        mprinterr("Unable to allocate space for cone_face\n");
        error_status++;
  }
  if ((broken_concave_face = (BROKEN_CONCAVE_FACE *) malloc (
                   NUM_FACE * natm_sel * sizeof (BROKEN_CONCAVE_FACE))) == NULL) {
        mprinterr("Unable to allocate space for broken_concave_face\n");
        error_status++;
  }
  if ((concave_cycle = (CONCAVE_CYCLE *) malloc (
                         NUM_CYCLE * natm_sel * sizeof (CONCAVE_CYCLE))) == NULL) {
        mprinterr("Unable to allocate space for concave_cycle\n");
        error_status++;
  }
  if ((cyclelist = (CYCLE *) malloc (
                         NUM_CYCLE * natm_sel * sizeof (CYCLE))) == NULL) {
        mprinterr("Unable to allocate space for cyclelist\n");
        error_status++;
  }
  if ((vertexlist = (VERTEX *) malloc (
                NUM_VERTEX * natm_sel * sizeof (VERTEX))) == NULL) {
        mprinterr("Unable to allocate space for vertexlist\n");
        error_status++;
  }
  if ((concave_edge_list = (EDGE *) malloc (
          NUM_EDGE * natm_sel * sizeof (EDGE))) == NULL) {
        mprinterr("Unable to allocate space for concave_edge_list\n");
        error_status++;
  }
  if ((convex_edge_list = (EDGE *) malloc (
          NUM_EDGE * natm_sel * sizeof (EDGE))) == NULL) {
        mprinterr("Unable to allocate space for convex_edge_list\n");
        error_status++;
  }
  if ((low_torus = (LOW_TORUS *) malloc (
                 NUM_TORUS * natm_sel * sizeof (LOW_TORUS))) == NULL) {
        mprinterr("Unable to allocate space for low_torus\n");
        error_status++;
  }
  // NOTE: cusp_edge must be initialized. Handled in action before 1st call 
  if ((cusp_edge = (CUSP_EDGE *) malloc (
                 NUM_EDGE * natm_sel * sizeof (CUSP_EDGE))) == NULL) {
        mprinterr("Unable to allocate space for cusp_edge\n");
        error_status++;
  }
  if ((cusp_pair = (CUSP_PAIR *) malloc (
                 NUM_CUSP * natm_sel * sizeof (CUSP_PAIR))) == NULL) {
        mprinterr("Unable to allocate space for cusp_pair\n");
        error_status++;
  }
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

  mprintf("    MOLSURF: [%s] Probe Radius=%.3lf\n",Mask1.maskString,probe_rad);
  if (rad_offset>0)
    mprintf("             Radii will be incremented by %.3lf\n",rad_offset);

  return 0;
}

// Molsurf::setup()
/** Set mask up for this parmtop.
  * currentParm is set in Action::Setup
  */
int Molsurf::setup() {
  if ( Mask1.SetupMask(currentParm,activeReference,debug) ) return 1;
  if (Mask1.None()) {
    mprintf("    Error: Molsurf::setup: Mask contains 0 atoms.\n");
    return 1;
  }

  mprintf("    MOLSURF: Calculating surface area for %i atoms.\n",Mask1.Nselected);
  // NOTE: If Mask is * dont include any solvent?

  // The ATOM structure is how molsurf organizes atomic data. Allocate
  // here and fill in parm info. Coords will be filled in during action. 
  if (atom!=NULL) free(atom);
  atom = (ATOM*) malloc(Mask1.Nselected * sizeof(ATOM));
  if (atom==NULL) {
    mprinterr("Error: Molsurf::Setup Could not allocate memory for ATOMs.\n");
    return 1;
  }
  // Set up parm info for atoms in mask
  double *Radii = currentParm->GB_radii();
  for (int maskidx = 0; maskidx < Mask1.Nselected; maskidx++) {
    int parmatom = Mask1.Selected[maskidx];
    int nres = currentParm->atomToResidue(parmatom);
    atom[maskidx].anum = parmatom + 1; // anum is for debug output only, atoms start from 1
    strcpy(atom[maskidx].anam,currentParm->names[parmatom]);
    atom[maskidx].rnum = nres + 1; // again for debug output only, residues start from 1
    strcpy(atom[maskidx].rnam,currentParm->resnames[nres]);
    atom[maskidx].pos[0] = 0;
    atom[maskidx].pos[1] = 0;
    atom[maskidx].pos[2] = 0;
    atom[maskidx].q = currentParm->charge[parmatom];
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


