/*
Copyright 2009, D. E. Shaw Research, LLC
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

* Redistributions of source code must retain the above copyright
  notice, this list of conditions, and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions, and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

* Neither the name of D. E. Shaw Research, LLC nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "molfile_plugin.h"
#include <sqlite3.h>

#include <ctype.h> /* for isspace() */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <vector>
#include <map>


#ifdef _MSC_VER
#include <windows.h>
static int fs_remove_file( const char * path ) {
  DeleteFileA( path );
  return 1; /* FIXME: check for failure */
}
#else
#include <unistd.h>
#include <errno.h>
static int fs_remove_file( const char * path ) {
    if (unlink(path)!=0 && errno!=ENOENT) {
        fprintf(stderr, "Removing %s failed: %s\n", path, strerror(errno));
        return 0;
    }
    return 1;
}
#endif

#ifndef M_PI
#define M_PI (3.1415926535897932385)
#endif

#ifndef M_PI_2
#define M_PI_2 (1.5707963267948966192)
#endif

struct element {
  double daltons;
  const char* abbreviation;
  const char* name;
};

/* url = "http://physics.nist.gov/cgi-bin/Elements/elInfo.pl?element=%d&context=noframes"%element */
static struct element amu[] = {
  {1.00794,"H","Hydrogen"},
  {4.002602,"He","Helium"},
  {6.941,"Li","Lithium"},
  {9.012182,"Be","Beryllium"},
  {10.811,"B","Boron"},
  {12.0107,"C","Carbon"},
  {14.0067,"N","Nitrogen"},
  {15.9994,"O","Oxygen"},
  {18.9984032,"F","Fluorine"},
  {20.1797,"Ne","Neon"},
  {22.989770,"Na","Sodium"},
  {24.3050,"Mg","Magnesium"},
  {26.981538,"Al","Aluminum"},
  {28.0855,"Si","Silicon"},
  {30.973761,"P","Phosphorus"},
  {32.065,"S","Sulfur"},
  {35.453,"Cl","Chlorine"},
  {39.0983,"K","Potassium"},
  {39.948,"Ar","Argon"},
  {40.078,"Ca","Calcium"},
  {44.955910,"Sc","Scandium"},
  {47.867,"Ti","Titanium"},
  {50.9415,"V","Vanadium"},
  {51.9961,"Cr","Chromium"},
  {54.938049,"Mn","Manganese"},
  {55.845,"Fe","Iron"},
  {58.6934,"Ni","Nickel"},
  {58.933200,"Co","Cobalt"},
  {63.546,"Cu","Copper"},
  {65.409,"Zn","Zinc"},
  {69.723,"Ga","Gallium"},
  {72.64,"Ge","Germanium"},
  {74.92160,"As","Arsenic"},
  {78.96,"Se","Selenium"},
  {79.904,"Br","Bromine"},
  {83.798,"Kr","Krypton"},
  {85.4678,"Rb","Rubidium"},
  {87.62,"Sr","Strontium"},
  {88.90585,"Y","Yttrium"},
  {91.224,"Zr","Zirconium"},
  {92.90638,"Nb","Niobium"},
  {95.94,"Mo","Molybdenum"},
  {101.07,"Ru","Ruthenium"},
  {102.90550,"Rh","Rhodium"},
  {106.42,"Pd","Palladium"},
  {107.8682,"Ag","Silver"},
  {112.411,"Cd","Cadmium"},
  {114.818,"In","Indium"},
  {118.710,"Sn","Tin"},
  {121.760,"Sb","Antimony"},
  {126.90447,"I","Iodine"},
  {127.60,"Te","Tellurium"},
  {131.293,"Xe","Xenon"},
  {132.90545,"Cs","Cesium"},
  {137.327,"Ba","Barium"},
  {138.9055,"La","Lanthanum"},
  {140.116,"Ce","Cerium"},
  {140.90765,"Pr","Praseodymium"},
  {144.24,"Nd","Neodymium"},
  {150.36,"Sm","Samarium"},
  {151.964,"Eu","Europium"},
  {157.25,"Gd","Gadolinium"},
  {158.92534,"Tb","Terbium"},
  {162.500,"Dy","Dysprosium"},
  {164.93032,"Ho","Holmium"},
  {167.259,"Er","Erbium"},
  {168.93421,"Tm","Thulium"},
  {173.04,"Yb","Ytterbium"},
  {174.967,"Lu","Lutetium"},
  {178.49,"Hf","Hafnium"},
  {180.9479,"Ta","Tantalum"},
  {183.84,"W","Tungsten"},
  {186.207,"Re","Rhenium"},
  {190.23,"Os","Osmium"},
  {192.217,"Ir","Iridium"},
  {195.078,"Pt","Platinum"},
  {196.96655,"Au","Gold"},
  {200.59,"Hg","Mercury"},
  {204.3833,"Tl","Thallium"},
  {207.2,"Pb","Lead"},
  {208.98038,"Bi","Bismuth"},
  {231.03588,"Pa","Protactinium"},
  {232.0381,"Th","Thorium"},
  {238.02891,"U","Uranium"}
};

static const int nelements = sizeof(amu)/sizeof(amu[0]);

static const char * 
find_element_by_atomic_number(int target) {
  if (target < 1) target=1;
  if (target >= nelements) target = nelements;
  return amu[target-1].abbreviation;
}

static int 
find_element_by_amu(double target) {
  int left = 0;
  int right = nelements-1;
  int swap;

  /* -----------------------------------------------
  // Knuth's binary search
  // ----------------------------------------------- */
  while(left <= right) {
    int mid = (left+right)/2;
    if (target> amu[mid].daltons) {
      left = mid + 1;
    } else if (target< amu[mid].daltons) {
      right = mid - 1;
    } else {
      /* Exact match (unlikely) */
      left = right = mid;
      return left+1;
    }
  }

  /* -----------------------------------------------
  // CAUTION: at this point, the meanings of 
  // left and right are switched (i.e. left >= right,
  // see the while() loop above if you don't believe me!
  // ----------------------------------------------- */
  swap = left;
  left = right;
  right = swap;

  if (left < 0) left = right;
  if (right > nelements-1) right = left;

  if (target - amu[left].daltons < amu[right].daltons - target) {
    return left+1;
  }

  return right+1;
}

static int table_size( sqlite3 * db, const char * tname, 
                       sqlite_int64 * count ) {
    sqlite3_stmt *stmt;
    char * buf;

    if (!tname) return 0;
    buf = (char *)malloc(strlen(tname) + 100);
    sprintf(buf, "select count() from '%s'", tname);
    if (sqlite3_prepare_v2(db, buf, -1, &stmt, NULL)) {
        free(buf);
        return 0;
    }
    if (SQLITE_ROW != sqlite3_step(stmt)) {
        sqlite3_finalize(stmt);
        free(buf);
        return 0;
    }
    if (count) *count = sqlite3_column_int64(stmt, 0);
    sqlite3_finalize(stmt);
    return 1; /* success */
}

static molfile_plugin_t plugin;
//static molfile_plugin_t append_plugin;

namespace {
    struct Handle {
        sqlite3 * db;
        int natoms;
        int frames_read;
        std::vector<int> from, to, gids;
        std::vector<int> glue_from, glue_to;
        std::vector<float> order;
        Handle(sqlite3 * _db, int n) 
        : db(_db), natoms(n), frames_read(0) {}
            
        ~Handle() {
            if (db && sqlite3_close(db)) 
                fprintf(stderr, "Error closing db: %s\n", sqlite3_errmsg(db));
        }
    };
}


  static void *open_file_read( const char *filename, const char *filetype,
                        int *natoms ) {

      sqlite3 * db;
      sqlite_int64 count;
      Handle * h;
      if (sqlite3_open_v2( filename, &db, SQLITE_OPEN_READONLY, NULL)) {
          fprintf(stderr, "Error opening dms at %s: %s", 
                  filename, sqlite3_errmsg(db));
          return NULL;
      }
      if (!table_size( db, "particle", &count )) return NULL;
      h = new Handle(db, count);
      *natoms = count;
      return h;
  }

static void close_file_read(void *v) { delete (Handle *)v; }

static int has_nonwhitespace( const char * p ) {
    for (; *p; ++p) {
        if (!isspace(*p)) return 1;
    }
    return 0;
}

#define GET_STRING( col, field ) do { \
    if ( col >= 0 ) { \
        const char * val = (const char *)sqlite3_column_text(stmt, col); \
        if (val) { \
            strncpy( atom->field, val, sizeof(atom->field) ); \
            atom->field[sizeof(atom->field)-1] = '\0'; \
        } \
    } \
} while(0)

#define GET_INT( col, field ) do { \
    if ( col >= 0 ) { \
        atom->field = sqlite3_column_int(stmt, col); \
    } \
} while (0)

#define GET_DOUBLE( col, field ) do { \
    if ( col >= 0 ) { \
        atom->field = sqlite3_column_double(stmt, col); \
    } \
} while (0)

  static
  int read_structure( void *v, int *optflags, molfile_atom_t *atoms ) {
      Handle *h = (Handle *)v;
      sqlite3 * db = h->db;
      sqlite3_stmt * stmt;
      static const char infosql[] = "pragma table_info(particle)";
      static const char loadsql[] = "select * from particle";
      molfile_atom_t * atom = NULL;

      int gid_col = -1;
      int name_col = -1;
      int anum_col = -1;
      int resid_col = -1;
      int resname_col = -1;
      int chain_col = -1;
      int segid_col = -1;
      int mass_col = -1;
      int charge_col = -1;
      int x_col = -1;
      int y_col = -1;
      int z_col = -1;
      int vx_col = -1;
      int vy_col = -1;
      int vz_col = -1;
      int occ_col = -1;
  
      int ncols=0;
      *optflags = 0;
  
      if (sqlite3_prepare_v2(db, infosql, -1, &stmt, NULL)) {
          fprintf(stderr, "Error getting particle info: %s", 
                  sqlite3_errmsg(db));
          return MOLFILE_ERROR;
      }
      atom = atoms;
      while (SQLITE_ROW==sqlite3_step(stmt)) {
          const char * name = (const char *)sqlite3_column_text(stmt,1);

          if      (!strcmp(name, "anum")) anum_col = ncols;
          else if (!strcmp(name, "id")) gid_col = ncols;
          else if (!strcmp(name, "name")) name_col = ncols;
          else if (!strcmp(name, "resid")) resid_col = ncols;
          else if (!strcmp(name, "resname")) resname_col = ncols;
          else if (!strcmp(name, "chain")) chain_col = ncols;
          else if (!strcmp(name, "segid")) segid_col = ncols;
          else if (!strcmp(name, "mass")) mass_col = ncols;
          else if (!strcmp(name, "charge")) charge_col = ncols;
          else if (!strcmp(name, "x")) x_col = ncols;
          else if (!strcmp(name, "y")) y_col = ncols;
          else if (!strcmp(name, "z")) z_col = ncols;
          else if (!strcmp(name, "vx")) vx_col = ncols;
          else if (!strcmp(name, "vy")) vy_col = ncols;
          else if (!strcmp(name, "vz")) vz_col = ncols;
          else if (!strcmp(name, "occupancy")) occ_col = ncols;
          ++ncols;
      }  
      sqlite3_finalize(stmt);

      if (anum_col >= 0)   *optflags |= MOLFILE_ATOMICNUMBER;
      if (charge_col >= 0) *optflags |= MOLFILE_CHARGE;
      if (mass_col >= 0)   *optflags |= MOLFILE_MASS;
      if (occ_col >= 0)    *optflags |= MOLFILE_OCCUPANCY;

      if (gid_col<0) {
          fprintf(stderr, "Missing id column in particle table\n");
          return MOLFILE_ERROR;
      }

      if (sqlite3_prepare_v2(db, loadsql, -1, &stmt, NULL)) {
          fprintf(stderr, "Error loading particle structure: %s",
                  sqlite3_errmsg(db));
          return MOLFILE_ERROR;
      }
      while (SQLITE_ROW==sqlite3_step(stmt)) {
          memset(atom, 0, sizeof(*atom));
          GET_STRING(    name_col, name );
          GET_INT(       anum_col, atomicnumber );
          GET_INT(      resid_col, resid );
          GET_STRING( resname_col, resname );
          GET_DOUBLE(    mass_col, mass );
          GET_DOUBLE(  charge_col, charge );
          GET_DOUBLE(     occ_col, occupancy );
          GET_STRING(   chain_col, chain );
          GET_STRING(   segid_col, segid );

          h->gids.push_back(sqlite3_column_int(stmt, gid_col));

          /* if the name is all space, and we have an atomic number, use the
          element name instead. */
          if (atom->atomicnumber) {
            if (!has_nonwhitespace(atom->name)) {
              const char *nm = 
                find_element_by_atomic_number(atom->atomicnumber);
              strncpy( atom->name, nm, sizeof(atom->name));
            }
          }
          ++atom;
      }
      sqlite3_finalize(stmt);
      return MOLFILE_SUCCESS;
  }

static
#if vmdplugin_ABIVERSION > 14
int read_bonds(void *v, int *nbonds, int **from, int **to, float** order,
               int **bondtype, int *nbondtypes, char ***bondtypename) {
#else
int read_bonds(void *v, int *nbonds, int **from, int **to, float** order) {
#endif
    Handle *h = (Handle *)v;
    sqlite3 * db = h->db;
    sqlite_int64 count = 0;
    if (table_size(h->db, "bond", &count) && count>0) {
        sqlite3_stmt * stmt;
        int with_order = 1;
        /* Try including order first */
        static const char * sql = "select p0,p1,\"order\" from bond";
        if (sqlite3_prepare_v2(db, sql, -1, &stmt, NULL)) {
            /* Fall back to just p0, p1 */
            static const char * sql = "select p0,p1 from bond";
            with_order = 0;
            if (sqlite3_prepare_v2(db, sql, -1, &stmt, NULL)) {
                fprintf(stderr, "Error reading bonds: %s", sqlite3_errmsg(db));
                return MOLFILE_ERROR;
            }
        }
        h->from.resize(count);
        h->to.resize(count);
        h->order.resize(count);

        *from = &h->from[0];
        *to   = &h->to[0];
        *order = &h->order[0];
#if vmdplugin_ABIVERSION > 14
        *bondtype = NULL;
        *nbondtypes = 0;
        *bondtypename = NULL;
  #endif

        std::map<int,int> gidmap;
        for (unsigned i=0; i<h->gids.size(); i++) gidmap[h->gids[i]]=i+1;

        for (int i=0; sqlite3_step(stmt)==SQLITE_ROW; i++) {
            std::map<int,int>::const_iterator iter;

            int from = sqlite3_column_int(stmt,0);
            int to   = sqlite3_column_int(stmt,1);
            h->order[i] = with_order ? sqlite3_column_int(stmt,2) : 1;
            if ((iter=gidmap.find(from))==gidmap.end()) {
                fprintf(stderr, "Illegal bond in row %d: %d-%d\n", 
                        i+1, from, to);
                return MOLFILE_ERROR;
            }
            h->from[i] = iter->second;
            if ((iter=gidmap.find(to))==gidmap.end()) {
                fprintf(stderr, "Illegal bond in row %d: %d-%d\n", 
                        i+1, from, to);
                return MOLFILE_ERROR;
            }
            h->to[i] = iter->second;
        }
        sqlite3_finalize(stmt);
    }
    *nbonds = count;
    return MOLFILE_SUCCESS;
}

static
int read_fictitious_bonds(void *v, int *nbonds, int **from, int **to) {
    Handle *h = (Handle *)v;
    sqlite3 * db = h->db;
    sqlite_int64 count = 0;
    if (table_size(h->db, "glue", &count) && count>0) {
        sqlite3_stmt * stmt;
        static const char * sql = "select p0,p1 from glue";
        if (sqlite3_prepare_v2(db, sql, -1, &stmt, NULL)) {
            fprintf(stderr, "Error reading bonds: %s", sqlite3_errmsg(db));
            return MOLFILE_ERROR;
        }
        h->glue_from.resize(count);
        h->glue_to.resize(count);

        *from = &h->glue_from[0];
        *to   = &h->glue_to[0];

        std::map<int,int> gidmap;
        for (unsigned i=0; i<h->gids.size(); i++) gidmap[h->gids[i]]=i+1;

        for (int i=0; sqlite3_step(stmt)==SQLITE_ROW; i++) {
            std::map<int,int>::const_iterator iter;

            int from = sqlite3_column_int(stmt,0);
            int to   = sqlite3_column_int(stmt,1);
            if ((iter=gidmap.find(from))==gidmap.end()) {
                fprintf(stderr, "Illegal glue in row %d: %d-%d\n", 
                        i+1, from, to);
                return MOLFILE_ERROR;
            }
            h->glue_from[i] = iter->second;
            if ((iter=gidmap.find(to))==gidmap.end()) {
                fprintf(stderr, "Illegal glue in row %d: %d-%d\n", 
                        i+1, from, to);
                return MOLFILE_ERROR;
            }
            h->glue_to[i] = iter->second;
        }
        sqlite3_finalize(stmt);
    }
    *nbonds = count;
    return MOLFILE_SUCCESS;
}

  static
  double dotprod(const double *x, const double *y) {
    return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
  }

  static
  int read_global_cell(Handle * h, molfile_timestep_t * ts) {
      sqlite3_stmt * stmt;
      sqlite_int64 count;
      int i,j;
      double cell[3][3];
      double cosAB, cosAC, cosBC;
      const double * A, * B, * C;
      if (!table_size(h->db, "global_cell", &count)) return MOLFILE_SUCCESS;
      if (count != 3) {
          fprintf(stderr, "Error: expected global_cell size of 3, got %d",
                  (int)count);
          return MOLFILE_ERROR;
      }
      if (sqlite3_prepare_v2(h->db,"select x,y,z from global_cell", -1, &stmt, NULL))  {
          fprintf(stderr, "Error compiling global_cell reader: %s",
                  sqlite3_errmsg(h->db));
          return MOLFILE_ERROR;
      }
      for (i=0; i<3; i++) {
          sqlite3_step(stmt);
          for (j=0; j<3; j++) {
              cell[i][j] = sqlite3_column_double(stmt,j);
          }
      }
      sqlite3_finalize(stmt);

      A = cell[0];
      B = cell[1];
      C = cell[2];

      /* store lengths */
      ts->A = sqrt(dotprod(A,A));
      ts->B = sqrt(dotprod(B,B));
      ts->C = sqrt(dotprod(C,C));

      if (ts->A==0 || ts->B==0 || ts->C==0) {
          ts->alpha = ts->beta = ts->gamma = 90.0;

      } else {
        /* compute angles */
        cosAB = dotprod(A,B)/(ts->A * ts->B);
        cosAC = dotprod(A,C)/(ts->A * ts->C);
        cosBC = dotprod(B,C)/(ts->B * ts->C);

        /* clamp */
        if (cosAB > 1.0) cosAB = 1.0; else if (cosAB < -1.0) cosAB = -1.0;
        if (cosAC > 1.0) cosAC = 1.0; else if (cosAC < -1.0) cosAC = -1.0;
        if (cosBC > 1.0) cosBC = 1.0; else if (cosBC < -1.0) cosBC = -1.0;

        /* convert to angles using asin to avoid nasty rounding when we are
        close to 90 degree angles. */
        ts->alpha = 90.0 - asin(cosBC) * 90.0 / M_PI_2; /* cosBC */
        ts->beta  = 90.0 - asin(cosAC) * 90.0 / M_PI_2; /* cosAC */
        ts->gamma = 90.0 - asin(cosAB) * 90.0 / M_PI_2; /* cosAB */
      }
      return MOLFILE_SUCCESS;
  }

  static
  int read_timestep(Handle * h, molfile_timestep_t *ts) {

      sqlite3_stmt * stmt;
      float * pos = ts->coords;
      float * vel = ts->velocities;
      const char * sql = vel ? "select x,y,z,vx,vy,vz from particle"
                             : "select x,y,z          from particle";
      if (sqlite3_prepare_v2( h->db, sql, -1, &stmt, NULL)) {
          fprintf(stderr, "Error reading timestep: %s\n", 
                  sqlite3_errmsg(h->db));
          return MOLFILE_ERROR;
      }
      while (SQLITE_ROW==sqlite3_step(stmt)) {
          int i;
          for (i=0; i<3; i++) *pos++ = sqlite3_column_double(stmt,i);
          if (vel) {
            for (i=0; i<3; i++) *vel++ = sqlite3_column_double(stmt,i+3);
          }
      }
      sqlite3_finalize(stmt);
      return read_global_cell(h,ts);
  }


  static
  int read_next_timestep( void *v, int atoms, molfile_timestep_t *ts) {
    Handle *h = (Handle *)v;
    if (h->frames_read++) return MOLFILE_EOF;
    return read_timestep(h, ts);
  }

  static
  int read_timestep2(void *v, ssize_t n, molfile_timestep_t *ts) {
    Handle *h = (Handle *)v;
    if (n!=0) return MOLFILE_EOF;
    return read_timestep(h, ts);
  }

  static
  int read_timestep_metadata( void *v, molfile_timestep_metadata_t *m) {
    /* FIXME: assume velocities */
    m->has_velocities = 1;
    m->count = 1; /* number of timesteps */
    m->avg_bytes_per_timestep = 10000; /* FIXME */
    return MOLFILE_SUCCESS;
  }

  static
  void *open_file_write(const char *path, const char *type, int natoms) {

      sqlite3 * db;
      Handle *h=NULL;
      /* delete existing file so we don't just append */
      if (!fs_remove_file(path)) return NULL;
      if (sqlite3_open(path, &db)) {
          fprintf(stderr, "Failed opening %s: %s\n", 
                  path, sqlite3_errmsg(db));
          return NULL;
      }
      h = new Handle(db, natoms);
      return h;
  }

  static
#if vmdplugin_ABIVERSION > 14
  int write_bonds(void *v, int nbonds, int *from, int *to, float *order,
                  int *bondtype, int nbondtypes, char **bondtypename) {
#else
  int write_bonds(void *v, int nbonds, int *from, int *to, float *order) {
#endif
      Handle *h = (Handle *)v;
      h->from.resize(nbonds);
      h->to.resize(nbonds);
      
      std::copy(from, from+nbonds, h->from.begin());
      std::copy(to,   to+nbonds,   h->to.begin());
      if (order) {
          h->order.resize(nbonds);
          std::copy(order, order+nbonds, h->order.begin());
      } else {
          h->order.clear();
          h->order.insert(h->order.begin(), nbonds, 1);
      }
      return MOLFILE_SUCCESS;
  }

  static
  int write_structure(void *v, int optflags, const molfile_atom_t *atoms) {
    Handle *h = (Handle *)v;
    sqlite3 * db = h->db;
    sqlite3_stmt * stmt;
    int i;

    sqlite3_exec(db, "begin", NULL, NULL, NULL);
    if (sqlite3_exec(db, 
        " create table if not exists particle (\n"
        "id integer primary key, \n"
        "anum integer,\n"
        "x float, y float, z float,\n"
        "vx float, vy float, vz float,\n"
        "mass float, charge float,\n"
        "name text, resname text, resid integer, chain text, segid text)", 
        NULL, NULL, NULL)) {
        fprintf(stderr, "Error creating particle table %s\n",
                sqlite3_errmsg(db));
        return MOLFILE_ERROR;
    }
    if (sqlite3_prepare_v2(db, "insert into particle values ("
                "?," /* id */
                "?," /* anum */
                "?,?,?,?,?,?," /* x,y,z, vx,vy,vz */
                "?,?,"  /* mass, charge */
                "?,?,?,?,?)", /* name, resname, resid, chain, segid */
                -1, &stmt, NULL)) {
        fprintf(stderr, "Error compiling insert particle statement %s\n",
                sqlite3_errmsg(db));
        return MOLFILE_ERROR;
    }
    for (i=0; i<h->natoms; i++) {
        const molfile_atom_t * a = atoms+i;
        int anum=0;
        if (optflags & MOLFILE_ATOMICNUMBER) {
            anum=a->atomicnumber;
        } else if (optflags & MOLFILE_MASS) {
            anum=find_element_by_amu(a->mass);
        }

        sqlite3_bind_int(stmt, 1, i );
        sqlite3_bind_int(stmt, 2, anum );
        sqlite3_bind_double( stmt, 9, 
                optflags & MOLFILE_MASS ?  a->mass : 0.0 );
        sqlite3_bind_double( stmt, 10, 
                optflags & MOLFILE_CHARGE ? a->charge : 0.0 );
        sqlite3_bind_text( stmt, 11, a->name, -1, SQLITE_STATIC );
        sqlite3_bind_text( stmt, 12, a->resname, -1, SQLITE_STATIC );
        sqlite3_bind_int(  stmt, 13, a->resid );
        sqlite3_bind_text( stmt, 14, a->chain, -1, SQLITE_STATIC );
        sqlite3_bind_text( stmt, 15, a->segid, -1, SQLITE_STATIC );
        
        if (sqlite3_step(stmt)!=SQLITE_DONE) {
            fprintf(stderr, "Error adding particle: %s\n",
                    sqlite3_errmsg(db));
            sqlite3_finalize(stmt);
            return MOLFILE_ERROR;
        }
        sqlite3_reset(stmt);
    }
    sqlite3_finalize(stmt);
    
    /* bonds */
    if (sqlite3_exec(db,
        "  create table if not exists bond (\n"
        "    p0 integer, p1 integer, 'order' integer)", NULL, NULL, NULL )) {
        fprintf(stderr, "Error creating bond table: %s\n",
                sqlite3_errmsg(db));
        return MOLFILE_ERROR;
    }
    if (sqlite3_prepare_v2(db, "insert into bond values (?,?,?)", -1, &stmt, NULL)) {
        fprintf(stderr, "Error compiling bond insert statement: %s\n",
                sqlite3_errmsg(db));
        return MOLFILE_ERROR;
    }
    for (unsigned i=0; i<h->from.size(); i++) {
        sqlite3_bind_int( stmt, 1, h->from[i]-1);
        sqlite3_bind_int( stmt, 2, h->to[i]-1);
        sqlite3_bind_int( stmt, 3, h->order[i]);
        if (sqlite3_step(stmt)!=SQLITE_DONE) {
            fprintf(stderr, "Error adding bond: %s\n", 
                    sqlite3_errmsg(db));
            sqlite3_finalize(stmt);
            return MOLFILE_ERROR;
        }
        sqlite3_reset(stmt);
    }
    sqlite3_finalize(stmt);
    sqlite3_exec(db, "commit", NULL, NULL, NULL);
    return MOLFILE_SUCCESS;
  }

  static
  void convert_to_homebox( const molfile_timestep_t * ts, double box[9]) {

    double A[3], B[3], C[3];

    /* Convert VMD's unit cell information */
    double cosBC = sin( ((90 - ts->alpha ) / 180) * M_PI );
    double cosAC = sin( ((90 - ts->beta  ) / 180) * M_PI );
    double cosAB = sin( ((90 - ts->gamma ) / 180) * M_PI );
    double sinAB = cos( ((90 - ts->gamma ) / 180) * M_PI );

    double Ax = ts->A;
    double Ay = 0;
    double Az = 0;
    double Bx = ts->B * cosAB;
    double By = ts->B * sinAB;
    double Bz = 0;
    double Cx,Cy,Cz;
    if (sinAB != 0) {
      Cx = cosAC;
      Cy = (cosBC - cosAC*cosAB) / sinAB;
      Cz = sqrt(1-Cx*Cx-Cy*Cy);
      Cx *= ts->C;
      Cy *= ts->C;
      Cz *= ts->C;
    } else {
      Cx=Cy=Cz=0;
    }
    A[0] = Ax; A[1] = Ay; A[2] = Az;
    B[0] = Bx; B[1] = By; B[2] = Bz;
    C[0] = Cx; C[1] = Cy; C[2] = Cz;

    /* put vectors in rows of box */
    box[0] = A[0]; box[1] = A[1]; box[2] = A[2];
    box[3] = B[0]; box[4] = B[1]; box[5] = B[2];
    box[6] = C[0]; box[7] = C[1]; box[8] = C[2];
  }

  static
  int write_timestep(void *v, const molfile_timestep_t *ts) {
    Handle *h = (Handle *)v;
    sqlite3 * db = h->db;
    sqlite3_stmt * stmt;
    double box[9];
    int i;
    if (h->frames_read) return MOLFILE_EOF;
    static const char sql[]=
        "update particle set x=?,y=?,z=?,vx=?,vy=?,vz=? where id=?";
    if (sqlite3_prepare_v2(db, sql, -1, &stmt, NULL)) {
        fprintf(stderr, "Error compiling update positions statement: %s\n",
                sqlite3_errmsg(db));
        return MOLFILE_ERROR;
    }
    sqlite3_exec(db, "begin", NULL, NULL, NULL);
    for (i=0; i<h->natoms; i++) {
        sqlite3_bind_double( stmt, 1, ts->coords[3*i  ] );
        sqlite3_bind_double( stmt, 2, ts->coords[3*i+1] );
        sqlite3_bind_double( stmt, 3, ts->coords[3*i+2] );
        double vx=0, vy=0, vz=0;
        if (ts->velocities) {
            vx=ts->velocities[3*i];
            vy=ts->velocities[3*i+1];
            vz=ts->velocities[3*i+2];
        }
        sqlite3_bind_double( stmt, 4, vx );
        sqlite3_bind_double( stmt, 5, vy );
        sqlite3_bind_double( stmt, 6, vz );
        sqlite3_bind_int( stmt, 7, i );
        if (sqlite3_step(stmt)!=SQLITE_DONE) {
            fprintf(stderr, "Error updating position: %s\n",
                    sqlite3_errmsg(db));
            sqlite3_finalize(stmt);
            return MOLFILE_ERROR;
        }
        sqlite3_reset(stmt);
    }
    sqlite3_finalize(stmt);

    if (sqlite3_exec(db,
        "  create table if not exists global_cell (\n"
        "  id integer primary key, x float, y float, z float)",
        NULL, NULL, NULL)) {
        fprintf(stderr, "Error creating global cell table %s\n",
                sqlite3_errmsg(db));
        return MOLFILE_ERROR;
    }
    sqlite3_exec(db, "delete from global_cell", NULL, NULL, NULL);
    if (sqlite3_prepare_v2(db, "insert into global_cell (x,y,z) values (?,?,?)",
                -1, &stmt, NULL)) {
        fprintf(stderr, "Error compiling global_cell statement: %s\n",
                sqlite3_errmsg(db));
        return MOLFILE_ERROR;
    }
    convert_to_homebox( ts, box );
    for (i=0; i<3; i++) {
        sqlite3_bind_double(stmt, 1, box[3*i  ]);
        sqlite3_bind_double(stmt, 2, box[3*i+1]);
        sqlite3_bind_double(stmt, 3, box[3*i+2]);
        if (sqlite3_step(stmt)!=SQLITE_DONE) {
            fprintf(stderr, "Error updating global cell: %s\n", 
                    sqlite3_errmsg(db));
            sqlite3_finalize(stmt);
            return MOLFILE_ERROR;
        }
        sqlite3_reset(stmt);
    }
    sqlite3_finalize(stmt);
    sqlite3_exec(db, "commit", NULL, NULL, NULL);
    return MOLFILE_SUCCESS;
  }

  static void close_file_write(void *v) { delete (Handle *)v; }


VMDPLUGIN_API int VMDPLUGIN_init(void) {
  memset(&plugin, 0, sizeof(plugin));
  plugin.abiversion = vmdplugin_ABIVERSION;
  plugin.type = MOLFILE_PLUGIN_TYPE;
  plugin.name = "dms";
  plugin.prettyname = "DESRES Molecular Structure";
  plugin.author = "D.E. Shaw Research, LLC: Justin Gullingsrud";
  plugin.majorv = 1;
  plugin.minorv = 0;
  plugin.is_reentrant = VMDPLUGIN_THREADSAFE;

  plugin.filename_extension = "dms";
  plugin.open_file_read = open_file_read;
  plugin.read_structure = read_structure;
  plugin.read_bonds = read_bonds;
  plugin.read_next_timestep = read_next_timestep;
#if defined(DESRES_READ_TIMESTEP2)
  plugin.read_timestep2 = read_timestep2;
  plugin.read_fictitious_bonds = read_fictitious_bonds;
#endif
  plugin.close_file_read = close_file_read;

  plugin.read_timestep_metadata = read_timestep_metadata;

  plugin.open_file_write = open_file_write;
  plugin.write_structure = write_structure;
  plugin.write_bonds = write_bonds;
  plugin.write_timestep = write_timestep;
  plugin.close_file_write = close_file_write;

  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_register( void *v, vmdplugin_register_cb cb) {
  cb( v, (vmdplugin_t *)&plugin);
  return VMDPLUGIN_SUCCESS;
}

VMDPLUGIN_API int VMDPLUGIN_fini(void) {
  return VMDPLUGIN_SUCCESS;
}

