// NAstruct 
#include <cstdlib>
#include <cstring>
#include "Action_NAstruct.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
NAstruct::NAstruct() {
  //fprintf(stderr,"NAstruct Con\n");
  resRange=NULL;
  outFilename=NULL;
} 

// DESTRUCTOR
NAstruct::~NAstruct() { 
  if (resRange!=NULL) delete resRange;
  for (std::list<AxisType*>::iterator it=RefCoords.begin(); it!=RefCoords.end(); it++) {
    if ( (*it)->F!=NULL ) delete (*it)->F;
    if ( (*it)->Name!=NULL) free( (*it)->Name );
    free(*it);
  }
}

// ------------------------- PRIVATE FUNCTIONS --------------------------------
/*
 * NAstruct::ID_base()
 * Return a number indicating if this is a NA base. If not recognized, 
 * return -1. 
 * NOTE: Currently based only on amber residue names. Will not recognize
 * non-standard bases.
 */
NAstruct::NAbaseType NAstruct::ID_base(char *resname) {
  //mprintf("        [%s]\n",resname);
  if (resname==NULL) return UNKNOWN_BASE;
  if (resname[0]=='D') {
    switch (resname[1]) {
      case 'A': return DA;
      case 'C': return DC;
      case 'G': return DG;
      case 'T': return DT;
    }
  } else if (resname[0]=='R') {
    switch (resname[1]) {
      case 'A': return RA;
      case 'C': return RC;
      case 'G': return RG;
      case 'U': return RU;
    }
  }
  return UNKNOWN_BASE;
}

/*
 * NAstruct::getRefCoords()
 * Allocate and set the coordinates to a standard ref. for the given base type.
 * Also set target atom names.
 * Coords taken from Olson et al. JMB (2001) 313, 229-237.
 */
NAstruct::AxisType *NAstruct::getRefCoords( NAbaseType btype) {
  AxisType *axis;
  Frame *F;

  F=NULL;
  axis = (AxisType *) malloc(sizeof(AxisType));
  if (axis==NULL) {
    mprintf("Error: NAstruct::getRefCoords: Could not allocate memory for coords.\n");
    return NULL;
  }
  axis->Name=NULL;
  switch (btype) {
    case DA :
      F = new Frame(11,NULL);
      axis->Name = (NAME*) malloc(11*sizeof(NAME));
      F->X[0 ]=-2.479000; F->X[1 ]= 5.346000; F->X[2 ]= 0.000000; strcpy(axis->Name[0 ],"C1'");
      F->X[3 ]=-1.291000; F->X[4 ]= 4.498000; F->X[5 ]= 0.000000; strcpy(axis->Name[1 ],"N9");
      F->X[6 ]= 0.024000; F->X[7 ]= 4.897000; F->X[8 ]= 0.000000; strcpy(axis->Name[2 ],"C8");
      F->X[9 ]= 0.877000; F->X[10]= 3.902000; F->X[11]= 0.000000; strcpy(axis->Name[3 ],"N7");
      F->X[12]= 0.071000; F->X[13]= 2.771000; F->X[14]= 0.000000; strcpy(axis->Name[4 ],"C5");
      F->X[15]= 0.369000; F->X[16]= 1.398000; F->X[17]= 0.000000; strcpy(axis->Name[5 ],"C6");
      F->X[18]= 1.611000; F->X[19]= 0.909000; F->X[20]= 0.000000; strcpy(axis->Name[6 ],"N6");
      F->X[21]=-0.668000; F->X[22]= 0.532000; F->X[23]= 0.000000; strcpy(axis->Name[7 ],"N1");
      F->X[24]=-1.912000; F->X[25]= 1.023000; F->X[26]= 0.000000; strcpy(axis->Name[8 ],"C2");
      F->X[27]=-2.320000; F->X[28]= 2.290000; F->X[29]= 0.000000; strcpy(axis->Name[9 ],"N3");
      F->X[30]=-1.267000; F->X[31]= 3.124000; F->X[32]= 0.000000; strcpy(axis->Name[10],"C4");
      break;
    case DC :
      F = new Frame(9,NULL);
      axis->Name = (NAME*) malloc(9*sizeof(NAME));
      F->X[0 ]=-2.477000; F->X[1 ]= 5.402000; F->X[2 ]= 0.000000; strcpy(axis->Name[0],"C1'");
      F->X[3 ]=-1.285000; F->X[4 ]= 4.542000; F->X[5 ]= 0.000000; strcpy(axis->Name[1],"N1");
      F->X[6 ]=-1.472000; F->X[7 ]= 3.158000; F->X[8 ]= 0.000000; strcpy(axis->Name[2],"C2");
      F->X[9 ]=-2.628000; F->X[10]= 2.709000; F->X[11]= 0.000000; strcpy(axis->Name[3],"O2");
      F->X[12]=-0.391000; F->X[13]= 2.344000; F->X[14]= 0.000000; strcpy(axis->Name[4],"N3");
      F->X[15]= 0.837000; F->X[16]= 2.868000; F->X[17]= 0.000000; strcpy(axis->Name[5],"C4");
      F->X[18]= 1.875000; F->X[19]= 2.027000; F->X[20]= 0.000000; strcpy(axis->Name[6],"N4");
      F->X[21]= 1.056000; F->X[22]= 4.275000; F->X[23]= 0.000000; strcpy(axis->Name[7],"C5");
      F->X[24]=-0.023000; F->X[25]= 5.068000; F->X[26]= 0.000000; strcpy(axis->Name[8],"C6");
      break;
    case DG :
      F = new Frame(12,NULL);
      axis->Name = (NAME*) malloc(12*sizeof(NAME));
      F->X[0 ]=-2.477000; F->X[1 ]= 5.399000; F->X[2 ]= 0.000000; strcpy(axis->Name[0],"C1'");
      F->X[3 ]=-1.289000; F->X[4 ]= 4.551000; F->X[5 ]= 0.000000; strcpy(axis->Name[1],"N9");
      F->X[6 ]= 0.023000; F->X[7 ]= 4.962000; F->X[8 ]= 0.000000; strcpy(axis->Name[2],"C8");
      F->X[9 ]= 0.870000; F->X[10]= 3.969000; F->X[11]= 0.000000; strcpy(axis->Name[3],"N7");
      F->X[12]= 0.071000; F->X[13]= 2.833000; F->X[14]= 0.000000; strcpy(axis->Name[4],"C5");
      F->X[15]= 0.424000; F->X[16]= 1.460000; F->X[17]= 0.000000; strcpy(axis->Name[5],"C6");
      F->X[18]= 1.554000; F->X[19]= 0.955000; F->X[20]= 0.000000; strcpy(axis->Name[6],"O6");
      F->X[21]=-0.700000; F->X[22]= 0.641000; F->X[23]= 0.000000; strcpy(axis->Name[7],"N1");
      F->X[24]=-1.999000; F->X[25]= 1.087000; F->X[26]= 0.000000; strcpy(axis->Name[8],"C2");
      F->X[27]=-2.949000; F->X[28]= 0.139000; F->X[29]=-0.001000; strcpy(axis->Name[9],"N2");
      F->X[30]=-2.342000; F->X[31]= 2.364000; F->X[32]= 0.001000; strcpy(axis->Name[10],"N3");
      F->X[33]=-1.265000; F->X[34]= 3.177000; F->X[35]= 0.000000; strcpy(axis->Name[11],"C4");
      break;
    case DT :
      F = new Frame(10,NULL);
      axis->Name = (NAME*) malloc(10*sizeof(NAME));
      F->X[0 ]=-2.481000; F->X[1 ]= 5.354000; F->X[2 ]=0.000000; strcpy(axis->Name[0],"C1'");
      F->X[3 ]=-1.284000; F->X[4 ]= 4.500000; F->X[5 ]=0.000000; strcpy(axis->Name[1],"N1");
      F->X[6 ]=-1.462000; F->X[7 ]= 3.135000; F->X[8 ]=0.000000; strcpy(axis->Name[2],"C2");
      F->X[9 ]=-2.562000; F->X[10]= 2.608000; F->X[11]=0.000000; strcpy(axis->Name[3],"O2");
      F->X[12]=-0.298000; F->X[13]= 2.407000; F->X[14]=0.000000; strcpy(axis->Name[4],"N3");
      F->X[15]= 0.994000; F->X[16]= 2.897000; F->X[17]=0.000000; strcpy(axis->Name[5],"C4");
      F->X[18]= 1.944000; F->X[19]= 2.119000; F->X[20]=0.000000; strcpy(axis->Name[6],"O4");
      F->X[21]= 1.106000; F->X[22]= 4.338000; F->X[23]=0.000000; strcpy(axis->Name[7],"C5");
      F->X[24]= 2.466000; F->X[25]= 4.961000; F->X[26]=0.001000; strcpy(axis->Name[8],"C5M");
      F->X[27]=-0.024000; F->X[28]= 5.057000; F->X[29]=0.000000; strcpy(axis->Name[9],"C6");
      break;
    case RA:
    case RC:
    case RG:
    case RU:
    case UNKNOWN_BASE:
      mprintf("Warning: NAstruct missing parameters for residue.\n");
  }
  
  axis->F = F;

  return axis;
}
// ----------------------------------------------------------------------------

/*
 * NAstruct::init()
 * Expected call: nastruct [resrange <range>] [out <filename>]
 * Dataset name will be the last arg checked for. Check order is:
 *    1) Keywords
 *    2) Masks
 *    3) Dataset name
 */
int NAstruct::init() {
  char *rangeArg;
  // Get keywords
  outFilename = A->getKeyString("out",NULL);
  rangeArg = A->getKeyString("resrange",NULL); 
  resRange = A->NextArgToRange(rangeArg); 

  // Get Masks
  //mask1 = A->getNextMask();
  //Mask1.SetMaskString(mask1);

  // Dataset
  // Add dataset to data file list

  //mprintf("    NAstruct: %s\n",Mask1.maskString);
  mprintf("    NAstruct: ");
  if (resRange==NULL)
    mprintf("Scanning all NA residues");
  else
    mprintf("Scanning residues %s",rangeArg);
  if (outFilename!=NULL)
    mprintf(", output to file %s",outFilename);
  mprintf("\n");

  return 0;
}

/*
 * NAstruct::setup()
 * Set up for this parmtop. Get masks etc.
 * P is set in Action::Setup
 */
int NAstruct::setup() {
  int res, refAtom, atom;
  std::list<int>::iterator it;
  AxisType *axis; 
  //AtomMask *Mask;

  // If range arg is NULL look for all NA residues.
  if (resRange==NULL) {
    resRange = new std::list<int>();
    for (res=0; res < P->nres; res++) {
      if (ID_base(P->resnames[res])!=UNKNOWN_BASE)
        resRange->push_back(res);
    }

  // For each residue in resRange determine if it is a NA
  } else {
    it=resRange->begin();
    while (it!=resRange->end()) {
      // User residues numbers start from 1
      (*it) = (*it) - 1;
      if (ID_base(P->ResidueName(*it))==UNKNOWN_BASE) 
        it = resRange->erase(it);
      else 
        it++;
    }
  }

  if (resRange->empty()) {
    mprintf("Error: NAstruct::setup: No NA residues found for %s\n",P->parmName);
    return 1;
  }

  // DEBUG
  mprintf("    NAstruct: NA res:");
  for (it=resRange->begin(); it!=resRange->end(); it++)
    mprintf(" %i",(*it)+1);
  mprintf("\n");

  // Set up reference coords for each NA residue
  // NOTE: Should just be combined above later
  for (it=resRange->begin(); it!=resRange->end(); it++) {
    axis = getRefCoords( ID_base(P->ResidueName(*it)) );
    if (axis==NULL) {
      mprintf("Error: NAstruct::setup: Could not get ref coords for %s\n",P->ResidueName(*it));
      return 1;
    }
    RefCoords.push_back( axis );

    // Set up a mask for this NA residue in this parm. The mask will contain
    // only those atoms which are defined in the reference coords.
    //Mask = new AtomMask();
    for (refAtom=0; refAtom < axis->F->natom; refAtom++) {
      res = -1; // Target atom
      for (atom=P->resnums[*it]; atom < P->resnums[(*it)+1]; atom++) {
        if ( strcmp(axis->Name[refAtom], P->names[atom])==0 ) {
          res = atom;  
          break;
        }
      }
      if (res==-1) {
        mprintf("Error:: NAstruct::setup: Ref atom %s not found in residue %i:%s\n",
                 axis->Name[refAtom], *it, P->resnames[(*it)]);
        return 1;
      }
      //Mask->AddAtom(res);
    } // End Loop over reference atoms
  } // End Loop over NA residues

  return 0;  
}

/*
 * NAstruct::action()
 */
int NAstruct::action() {
  //double Ang;

  //Ang=F->ANGLE(&Mask1,&Mask2,&Mask3);

  //ang->Add(currentFrame, &Ang);

  //fprintf(outfile,"%10i %10.4lf\n",currentFrame,Ang);
  
  return 0;
} 


