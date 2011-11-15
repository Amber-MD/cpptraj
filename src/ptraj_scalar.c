#include <string.h> // strcmp
#include "ptraj_scalar.h"
#include "ptraj_common.h"
// -----------------------------------------------------------------------------
// ptrajCopyScalar()
scalarInfo *ptrajCopyScalar(scalarInfo **scalarinp) {
  scalarInfo *scalar, *scalarin;
  /*
   *  Make a copy of the input scalar pointer.
   */
  scalarin = *scalarinp;
  scalar = (scalarInfo*)safe_malloc(sizeof(scalarInfo));
  INITIALIZE_scalarInfo(scalar);
  scalar->mode = scalarin->mode;
  scalar->totalFrames = scalarin->totalFrames;
  scalar->frame = scalarin->frame;
  scalar->mean = scalarin->mean;
  scalar->stddev = scalarin->stddev;
  scalar->max = scalarin->max;
  scalar->min = scalarin->min;
  scalar->atom1 = scalarin->atom1;
  scalar->atom2 = scalarin->atom2;
  scalar->atom3 = scalarin->atom3;
  scalar->atom4 = scalarin->atom4;
  scalar->atom5 = scalarin->atom5;
  scalar->mask1 = scalarin->mask1;
  scalar->mask2 = scalarin->mask2;
  scalar->mask3 = scalarin->mask3;
  scalar->mask4 = scalarin->mask4;
  scalar->mask5 = scalarin->mask5;
  scalar->cos = scalarin->cos;
  scalar->sin = scalarin->sin;
  scalar->value = scalarin->value;
  scalar->results = scalarin->results;
  //scalar->action = scalarin->action;
  scalar->name = scalarin->name;
  scalar->filename = scalarin->filename;
  return scalar;
}

// scalarStackGetName()
scalarInfo *scalarStackGetName(stackType **scalarStackP, char *name)
{
  stackType *s;
  scalarInfo *info, *match;

  match = NULL;
  for (s = *scalarStackP; s != NULL; s = s->next) {
    info = (scalarInfo *) s->entry;
    if ( strcmp(info->name, name) == 0 )
      match = info;
  }
  return (match);
}

// matrixInfoStackGetName()
transformMatrixInfo *matrixInfoStackGetName(stackType **matrixStackP, char *name)
{
  stackType *s;
  transformMatrixInfo *info, *match;

  match = NULL;
  for (s = *matrixStackP; s != NULL; s = s->next) {
    info = (transformMatrixInfo *) s->entry;
    if ( strcmp(info->name, name) == 0 )
      match = info;
  }
  return (match);
}

// vectorInfoStackGetName()
transformVectorInfo *vectorInfoStackGetName(stackType **vectorStackP, char *name) {
  stackType *vStack;
  transformVectorInfo *vinfo, *match;

  match = NULL;
  for (vStack = *vectorStackP; vStack != NULL; vStack = vStack->next) {
    vinfo = (transformVectorInfo *) vStack->entry;
    if ( strcmp(vinfo->name, name) == 0 ) {
      match = vinfo;
      break;
    }
  }
  return match;
}

// modesInfoStackGetName()
modesInfo *modesInfoStackGetName(stackType **modesStackP, char *name)
{
  stackType *s;
  modesInfo *info, *match;

  match = NULL;
  for (s = *modesStackP; s != NULL; s = s->next) {
    info = (modesInfo *) s->entry;
    if ( strcmp(info->name, name) == 0 )
      match = info;
  }
  return (match);
}

// distindex()
int distindex(int mask1tot, int i, int j) {

  /* Assure in call that i < j.
   * Returns index 0 for A-B, 1 for A-C, 2 for A-D, ...,
   *   NOT including main diagonal (i.e. A-A, ...)
   */

  return (i * mask1tot - (i * (i+1) / 2) + (j - i - 1));
}


