#ifndef INC_MOLSURF_H
#define INC_MOLSURF_H
#ifdef __cplusplus
extern "C" {
#endif

// ---------- CPPTRAJ includes
#include "Name.h" // NAME

// ---------- DEFINES ----------------------------------------------------------
#define REAL_T double
#define MAXAT_EDGE 30
#define NAME_SIZE 8

// ---------- DATA structures --------------------------------------------------
typedef REAL_T POINT[3];
typedef struct atom {
        POINT pos;
        REAL_T q, rad;
        char anam[NAME_SIZE];
        char rnam[NAME_SIZE];
        int anum;
        int rnum;
        int buried;                   // 1 = buried
        int neighbor_start, n_neighbors;
        int upper_start, n_upper;     // points to neighbor_torus struct
        int torus_start;              // points to torus struct
        // note that the toruslist[] only contains pairs of atoms 
        // in ascending order if you want all tori that contain an 
        // atom you have to look through the whole toruslist 
        int ntorus; 
        int n_convex_edges;           // convex edges associated with atom
        int convex_edges[MAXAT_EDGE]; // convex edges associated with atom 
        int n_cycles;                 // cycles of edges associated with atom
        int cycle_start;              // points to start of cyclelist for the atom
        REAL_T area;                  // accessible surf area associated with the atom
} ATOM;

// ---------- FUNCTIONS --------------------------------------------------------
REAL_T molsurf(REAL_T, ATOM*, int);
 
#ifdef __cplusplus
}
#endif
#endif
