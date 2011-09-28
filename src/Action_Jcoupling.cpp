// Jcoupling
#include "Action_Jcoupling.h"
#include <cstdlib> //getenv
#include <cstdio> //sscanf,sprintf
#include <cstring> //strcpy, strlen
#include <cmath> //cos
#include "CpptrajStdio.h"
#include "Constants.h" // DEGRAD, RADDEG

// CONSTRUCTOR
Jcoupling::Jcoupling() {
  //fprintf(stderr,"Jcoupling Con\n");
  Nconstants=0;
} 

// DESTRUCTOR
Jcoupling::~Jcoupling() {
  //fprintf(stderr,"Jcoupling Destructor.\n");
  for (std::map<std::string,karplusConstantList>::iterator reslist = KarplusConstants.begin();
                                                    reslist != KarplusConstants.end();
                                                    reslist++)
  {
    std::vector<karplusConstant> *currentList = (*reslist).second;
    delete currentList;
  }
  // DEBUG - Close output
  outputfile.CloseFile();
}

/* Jcoupling::loadKarplus()
 * Load Karplus parameters from input file.
 * Expected format:
 *   {type}<+|-| ><a[4]><+|-| ><b[4]><+|-| ><c[4]><+|-| ><d[4]><A[6]><B[6]><C[6]>{<D[6]>}
 *   <reslabel[4]>* 
 */
int Jcoupling::loadKarplus(char* filename) {
  char buffer[512],residue[5];
  char *end, *ptr;
  int i;
  PtrajFile KarplusFile;
  karplusConstant KC;
  std::vector<karplusConstant> *currentResList=NULL;
  std::string CurrentRes;
  std::map<std::string,karplusConstantList>::iterator reslist;

  if (KarplusFile.SetupFile(filename,READ,debug)) {
    mprinterr("Error: jcoupling: Could not read Karplus file %s\n",filename);
    return 1;
  }
  // residue is only for reading in 4 chars for residue names
  residue[4]='\0'; 
  KarplusFile.OpenFile();
  // Read through all lines of the file
  while (KarplusFile.IO->Gets(buffer,512)==0) {
    // Skip newlines and comments
    if (buffer[0]=='\n' || buffer[0]=='#') continue;
    ptr=buffer;
    // First char is optional type. If optional type is C, then the Karplus 
    // function specified in Perez et al. JACS (2001) 123 will be used, and 
    // A, B, and C will be taken as C0, C1, and C2.
    if(ptr[0]=='C') {
      KC.type=1;
      ptr++;
    } else {
      KC.type=0;
    }
    // Read atom names with optional preceding character (+, -)
    for (i=0; i<4; i++) {
      if      (*ptr=='+') KC.offset[i]=1;
      else if (*ptr=='-') KC.offset[i]=-1;
      else                     KC.offset[i]=0;
      ptr++;
      KC.atomName[i][0]=*ptr; ptr++;
      KC.atomName[i][1]=*ptr; ptr++;
      KC.atomName[i][2]=*ptr; ptr++;
      KC.atomName[i][3]=*ptr; ptr++;
      KC.atomName[i][4]='\0';
      //mprintf("DEBUG:\tAtomName %i [%s]\n",i,KC.atomName[i]);
    }
    // Read parameters
    // NOTE: Using sscanf here instead of atof since the 4th parameter is
    //       optional, behavior is undefined for accessing uninitialized
    //       portion of buffer.
    i = sscanf(ptr, "%6lf%6lf%6lf%6lf",KC.C,KC.C+1,KC.C+2,KC.C+3);
    if (i<3) {
      mprintf("Error: jcoupling: Expected at least 3 Karplus parameters, got %i\n",i);
      mprintf("       Line: [%s]\n",buffer);
      return 1;
    } else if (i==3) KC.C[3]=0.0;
    KC.C[3]*=DEGRAD;
    // Place the read-in karplus constants in a map indexed by residue name 
    // so that all karplus constants for a given residue are in one place. 
    KarplusFile.IO->Gets(buffer,512);
    // end will hold the end of the read-in buffer string
    end = buffer + strlen(buffer);
    for (ptr = buffer; ptr < end; ptr+=4) {
      if (*ptr=='\n') continue;
      residue[0] = ptr[0];
      residue[1] = ptr[1];
      residue[2] = ptr[2];
      residue[3] = ptr[3];
      CurrentRes.assign(residue);
      //mprintf("DEBUG:\t[%s]\n",CurrentRes.c_str());
      reslist = KarplusConstants.find(CurrentRes);
      // If list does not exist for residue yet, create it.
      // Otherwise, retrieve it
      if (reslist == KarplusConstants.end() ) {
        currentResList = new std::vector<karplusConstant>;
        KarplusConstants.insert(reslist, 
                                std::pair<std::string,karplusConstantList>(
                                  CurrentRes,currentResList));
      } else
        currentResList = (*reslist).second;

      currentResList->push_back(KC);
      Nconstants++;
    } // END loop over residues in residue line 
  } // END Gets over input file
  KarplusFile.CloseFile();
  // DEBUG - Print out all parameters
  if (debug>0) {
      mprintf("    KARPLUS PARAMETERS:\n");
      for (reslist=KarplusConstants.begin(); reslist!=KarplusConstants.end(); reslist++) {
        mprintf("\t[%4s]\n",(*reslist).first.c_str());
        for (std::vector<karplusConstant>::iterator kc=currentResList->begin();
                                                    kc!=currentResList->end();
                                                    kc++) {
          mprintf("\t\t%1i",(*kc).type);
          mprintf(" %4s",(*kc).atomName[0]);
          mprintf(" %4s",(*kc).atomName[1]);
          mprintf(" %4s",(*kc).atomName[2]);
          mprintf(" %4s",(*kc).atomName[3]);
          mprintf(" %i %i %i %i",(*kc).offset[0],(*kc).offset[1],(*kc).offset[2],(*kc).offset[3]);
          mprintf(" %6.2lf %6.2lf %6.2lf %6.2lf\n",(*kc).C[0],(*kc).C[1],(*kc).C[2],(*kc).C[3]);
        }
      }
  }
  return 0;
}

// -----------------------------------------------------------------------------
/* Jcoupling::init()
 * Expected call: jcoupling <mask1> [outfile <filename>]
 * Dataset name will be the last arg checked for. Check order is:
 *    1) Keywords
 *    2) Masks
 *    3) Dataset name
 */
int Jcoupling::init( ) {
  char *mask1;
  char *outfilename;
  char *env;
  char *karpluspath;

  // Get Keywords
  outfilename = A->getKeyString("outfile",NULL);

  // Get Masks
  mask1 = A->getNextMask();
  //fprintf(stdout,"    Mask 1: %s\n",mask1);
  Mask1.SetMaskString(mask1);

  // Dataset setup 
  // Add dataset to data file list

  // Get Karplus parameters from file.
  karpluspath=NULL;
  // First look in $AMBERHOME/dat/Karplus.txt
  env = getenv("AMBERHOME");
  if (env==NULL) {
    mprintf("    Warning: jcoupling: AMBERHOME not set.\n");
  } else {
    karpluspath = (char*) malloc( (strlen(env) + 17) * sizeof(char));
    sprintf(karpluspath,"%s/dat/Karplus.txt",env);
    if (!fileExists(karpluspath)) {
      mprintf("    Warning: jcoupling: %s not found\n",karpluspath);
      free(karpluspath);
      karpluspath=NULL;
    }
  }
  // Second look for $KARPLUS
  if (karpluspath==NULL) {
    env = getenv("KARPLUS");
    if (env==NULL) {
      mprinterr("Error: jcoupling: Karplus parameters $AMBERHOME/dat/Karplus.txt not found\n");
      mprinterr("       and KARPLUS environment variable not set.\n");
      return 1;
    }
    karpluspath = (char*) malloc( (strlen(env) + 1) * sizeof(char));
    strcpy(karpluspath, env);
  }
  // Load Karplus parameters
  if (loadKarplus(karpluspath)) {
    free(karpluspath);
    return 1;
  }

  mprintf("    J-COUPLING: Searching for dihedrals in mask [%s].\n",Mask1.maskString);
  mprintf("                Using Karplus parameters in \"%s\"\n",karpluspath);
  mprintf("                %i parameters found for %i residues.\n",Nconstants,
          KarplusConstants.size());
  mprintf("                Writing output to %s\n",outfilename);
  free(karpluspath);

  // DEBUG - Open output
  outputfile.SetupFile(outfilename,WRITE,debug);
  outputfile.OpenFile();

  return 0;
}

/* Jcoupling::setup()
 * Set up a j-coupling calculation for dihedrals defined by atoms within
 * the mask.
 */
int Jcoupling::setup() {
  std::string resName;
  std::vector<karplusConstant> *currentResList=NULL;
  int startatom,endatom,MaxResidues;
  jcouplingInfo JC;

  if ( Mask1.SetupCharMask(P,debug) ) return 1;
  if (Mask1.None()) {
    mprintf("    Error: Jcoupling::setup: Mask specifies no atoms.\n");
    return 1;
  }

  // For each residue, set up 1 jcoupling calc for each parameter defined in
  // KarplusConstants for this resdiue. Only set up the Jcoupling calc if all
  // atoms involved are present in the mask.
  MaxResidues = P->nres;
  if (P->finalSoluteRes > 0) MaxResidues = P->finalSoluteRes;
  for (int residue=0; residue < MaxResidues; residue++) {
    resName.assign(P->resnames[residue]);
    std::map<std::string,karplusConstantList>::iterator reslist = KarplusConstants.find(resName);
    // If list does not exist for residue, skip it.
    if (reslist == KarplusConstants.end() ) {
      mprintf("    Warning: Jcoupling::setup: Karplus parameters not found for residue [%i:%s]\n",
              residue+1, resName.c_str());
      continue;
    }
    startatom = P->resnums[residue];
    endatom = P->resnums[residue+1];
    currentResList = (*reslist).second;
    // For each parameter set in the list find the corresponding atoms.
    for (std::vector<karplusConstant>::iterator kc = currentResList->begin();
                                                kc != currentResList->end();
                                                kc++) 
    {
      // Init jcoupling info
      // NOTE: Should C[] just point to inside KarplusConstants?
      JC.residue = residue;
      JC.atom[0] = -1;
      JC.atom[1] = -1;
      JC.atom[2] = -1;
      JC.atom[3] = -1;
      JC.C[0]=0;
      JC.C[1]=0;
      JC.C[2]=0;
      JC.C[3]=0;
      JC.type=(*kc).type;
      // For each atom in the dihedral specified in this Karplus constant, find
      // corresponding atoms in parm.
      for (int idx=0; idx < 4; idx++) {
        for (int resatom = startatom; resatom < endatom; resatom++) {
          if (strcmp(P->names[resatom], (*kc).atomName[idx])==0) {
            JC.atom[idx] = resatom;
          }
        }
        // At the same time, set the Karplus constant
        JC.C[idx] = (*kc).C[idx];
      }
      // Check that all atoms were found
      for (int idx=0; idx < 4; idx++) {
        if (JC.atom[idx]==-1) {
          mprinterr("Error: jcoupling::setup: Atom %4s:%i not found for residue %i\n",
                    (*kc).atomName[idx], idx, residue);
          return 1;
        }
      }
      // Check that all the atoms involved in this Jcouple dihedral are
      // in the atom mask.
      if (!Mask1.AtomInCharMask(JC.atom[0])) continue;
      if (!Mask1.AtomInCharMask(JC.atom[1])) continue;
      if (!Mask1.AtomInCharMask(JC.atom[2])) continue;
      if (!Mask1.AtomInCharMask(JC.atom[3])) continue;
      // Add this jcoupling info to the list
      JcouplingInfo.push_back(JC);
    } // END loop over karplus parameters for this residue
  } // END loop over all residues

  // Print info for this parm
  mprintf("    J-COUPLING: [%s] Will calculate J-coupling for %u dihedrals.\n",Mask1.maskString,
          JcouplingInfo.size());
  if (JcouplingInfo.size()==0) {
    mprintf("    Warning: No dihedrals found for J-coupling calculation!\n");
    mprintf("             Check that all atoms of dihedrals are included in mask [%s]\n",
            Mask1.maskString);
    mprintf("             and/or that dihedrals are defined in Karplus parameter file.\n");
    return 1;
  }
  // DEBUG
  if (debug>0) {
    startatom=0;
    for (std::vector<jcouplingInfo>::iterator jc = JcouplingInfo.begin();
                                              jc !=JcouplingInfo.end();
                                              jc++) 
    {
      mprintf("%8i [%i:%4s]",startatom+1,(*jc).residue,P->resnames[(*jc).residue]);
      mprintf(" %6i:%-4s",(*jc).atom[0],P->names[(*jc).atom[0]]);
      mprintf(" %6i:%-4s",(*jc).atom[1],P->names[(*jc).atom[1]]);
      mprintf(" %6i:%-4s",(*jc).atom[2],P->names[(*jc).atom[2]]);
      mprintf(" %6i:%-4s",(*jc).atom[3],P->names[(*jc).atom[3]]);
      mprintf(" %6.2lf%6.2lf%6.2lf%6.2lf %i\n",(*jc).C[0],(*jc).C[1],(*jc).C[2],(*jc).C[3],
              (*jc).type);
      startatom++;
    }
  }    
  return 0;  
}

/* Jcoupling::action()
 * For each dihedral defined in JcouplingInfo, perform the dihedral and
 * Jcoupling calculation.
 */
int Jcoupling::action() {
  double phi,J;
  double phitemp,C0,C1,C2,C3;
  int residue;
  char buffer[53];

  for (std::vector<jcouplingInfo>::iterator jc = JcouplingInfo.begin();
                                            jc !=JcouplingInfo.end();
                                            jc++)
  {
    phi = F->DIHEDRAL( (*jc).atom[0], (*jc).atom[1], (*jc).atom[2], (*jc).atom[3] );
    C0 = (*jc).C[0];
    C1 = (*jc).C[1];
    C2 = (*jc).C[2];
    C3 = (*jc).C[3];
    if ((*jc).type==1) {
      //J = JcouplingC((*jc).C, phi);
      phitemp = phi + C3;
      J = C0 + (C1 * cos(phitemp)) + (C2 * cos(phitemp * 2.0)); 
    } else {
      //J = JcouplingABC((*jc).C, phi);
      phitemp = cos( phi + C3 );
      J = (C0 * phitemp * phitemp) + (C1 * phitemp) + C2;
    }

    residue = (*jc).residue;
    // DEBUG - output
    sprintf(buffer,"%5i %4s%4s%4s%4s%4s%12lf%12lf\n",residue+1,P->resnames[residue],
            P->names[(*jc).atom[0]],P->names[(*jc).atom[1]],P->names[(*jc).atom[2]],
            P->names[(*jc).atom[3]],phi*RADDEG,J);
    outputfile.IO->Write(buffer,1,51);
    //mprintf("%5i %4s",residue+1,P->resnames[residue]);
    //mprintf("%4s",P->names[(*jc).atom[0]]);
    //mprintf("%4s",P->names[(*jc).atom[1]]);
    //mprintf("%4s",P->names[(*jc).atom[2]]);
    //mprintf("%4s",P->names[(*jc).atom[3]]);
    //mprintf("%12lf%12lf\n",phi,J);
  }

  //fprintf(outfile,"%10i %10.4lf\n",currentFrame,D);
  
  return 0;
} 
