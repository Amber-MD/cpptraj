// Jcoupling
#include "Action_Jcoupling.h"
#include <cstdlib> //getenv
#include <cstdio> //sscanf
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

// Jcoupling::loadKarplus()
/** Load Karplus parameters from input file.
  * Expected format:
  * - {type}<+|-| ><a[4]><+|-| ><b[4]><+|-| ><c[4]><+|-| ><d[4]><A[6]><B[6]><C[6]>{<D[6]>}
  *   <reslabel[4]>* 
  * \return 0 on success, 1 on error
  */
int Jcoupling::loadKarplus(std::string filename) {
  char buffer[512],residue[5];
  char *end, *ptr;
  int i;
  CpptrajFile KarplusFile;
  karplusConstant KC;
  std::vector<karplusConstant> *currentResList=NULL;
  std::string CurrentRes;
  std::map<std::string,karplusConstantList>::iterator reslist;

  if (filename.empty()) {
    mprinterr("Error: jcoupling: Could not find Karplus parameter file.\n");
    return 1;
  }
  if (KarplusFile.SetupFile((char*)filename.c_str(),READ,debug)) {
    mprinterr("Error: jcoupling: Could not read Karplus parameter file %s\n",
              filename.c_str());
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
// Jcoupling::init()
/** Expected call: jcoupling <mask1> [outfile <filename>]
  */
int Jcoupling::init( ) {
  char *mask1;
  char *outfilename;
  char *env;
  std::string karpluspath;

  // Get Keywords
  outfilename = actionArgs.getKeyString("outfile",NULL);

  // Get Masks
  mask1 = actionArgs.getNextMask();
  //fprintf(stdout,"    Mask 1: %s\n",mask1);
  Mask1.SetMaskString(mask1);

  // Dataset setup 
  // Add dataset to data file list

  // Get Karplus parameters from file.
  karpluspath.clear();
  // First look in $AMBERHOME/dat/Karplus.txt
  env = getenv("AMBERHOME");
  if (env==NULL) {
    mprintf("    Warning: jcoupling: AMBERHOME not set.\n");
  } else {
    karpluspath.assign(env);
    karpluspath += "/dat/Karplus.txt";
    if (!fileExists((char*)karpluspath.c_str())) {
      mprintf("    Warning: jcoupling: %s not found\n",karpluspath.c_str());
      karpluspath.clear(); 
    }
  }
  // Second look for $KARPLUS
  if (karpluspath.empty()) {
    env = getenv("KARPLUS");
    if (env==NULL) {
      mprinterr("Error: jcoupling: Karplus parameters $AMBERHOME/dat/Karplus.txt not found\n");
      mprinterr("       and KARPLUS environment variable not set.\n");
      return 1;
    }
    karpluspath.assign(env);
  }
  // Load Karplus parameters
  if (loadKarplus(karpluspath)) 
    return 1;

  mprintf("    J-COUPLING: Searching for dihedrals in mask [%s].\n",Mask1.MaskString());
  mprintf("                Using Karplus parameters in \"%s\"\n",karpluspath.c_str());
  mprintf("                %i parameters found for %i residues.\n",Nconstants,
          KarplusConstants.size());
  mprintf("                Writing output to %s\n",outfilename);

  // DEBUG - Open output
  outputfile.SetupFile(outfilename,WRITE,debug);
  outputfile.OpenFile();

  return 0;
}

// Jcoupling::setup()
/** Set up a j-coupling calculation for dihedrals defined by atoms within
  * the mask.
  */
int Jcoupling::setup() {
  std::string resName;
  std::vector<karplusConstant> *currentResList=NULL;
  int MaxResidues;
  jcouplingInfo JC;

  if ( currentParm->SetupCharMask(Mask1, activeReference) ) return 1;
  if (Mask1.None()) {
    mprinterr("    Error: Jcoupling::setup: Mask specifies no atoms.\n");
    return 1;
  }
  // If JcouplingInfo has already been set up, print a warning and reset for
  // new parm.
  if (JcouplingInfo.size() > 0) {
    mprintf("    Warning: Jcoupling has been set up for another parm.\n");
    mprintf("             Resetting jcoupling info for new parm %s\n",currentParm->parmName);
    JcouplingInfo.clear();
  }

  // For each residue, set up 1 jcoupling calc for each parameter defined in
  // KarplusConstants for this residue. Only set up the Jcoupling calc if all
  // atoms involved are present in the mask.
  MaxResidues = currentParm->FinalSoluteRes();
  for (int residue=0; residue < MaxResidues; residue++) {
    resName.assign(currentParm->ResidueName(residue));
    std::map<std::string,karplusConstantList>::iterator reslist = KarplusConstants.find(resName);
    // If list does not exist for residue, skip it.
    if (reslist == KarplusConstants.end() ) {
      mprintf("    Warning: Jcoupling::setup: Karplus parameters not found for residue [%i:%s]\n",
              residue+1, resName.c_str());
      continue;
    }
    currentResList = (*reslist).second;
    // For each parameter set in the list find the corresponding atoms.
    for (std::vector<karplusConstant>::iterator kc = currentResList->begin();
                                                kc != currentResList->end();
                                                kc++) 
    {
      // Init jcoupling info
      // Constants will point inside KarplusConstants
      JC.residue = residue;
      JC.atom[0] = -1;
      JC.atom[1] = -1;
      JC.atom[2] = -1;
      JC.atom[3] = -1;
      JC.C=(*kc).C;
      JC.type=(*kc).type;
      // For each atom in the dihedral specified in this Karplus constant, find
      // corresponding atoms in parm. 
      for (int idx=0; idx < 4; idx++) 
        JC.atom[idx] = currentParm->FindAtomInResidue(residue+(*kc).offset[idx],
                                                      (*kc).atomName[idx]       );
      // Check that all atoms were found
      for (int idx=0; idx < 4; idx++) {
        if (JC.atom[idx]==-1) {
          mprinterr("Error: jcoupling::setup: Atom %4s:%i not found for residue %i\n",
                    (*kc).atomName[idx], idx, residue+(*kc).offset[idx]);
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
  mprintf("    J-COUPLING: [%s] Will calculate J-coupling for %u dihedrals.\n",Mask1.MaskString(),
          JcouplingInfo.size());
  if (JcouplingInfo.size()==0) {
    mprintf("    Warning: No dihedrals found for J-coupling calculation!\n");
    mprintf("             Check that all atoms of dihedrals are included in mask [%s]\n",
            Mask1.MaskString());
    mprintf("             and/or that dihedrals are defined in Karplus parameter file.\n");
    return 1;
  }
  // DEBUG
  if (debug>0) {
    MaxResidues=0;
    for (std::vector<jcouplingInfo>::iterator jc = JcouplingInfo.begin();
                                              jc !=JcouplingInfo.end();
                                              jc++) 
    {
      mprintf("%8i [%i:%4s]",MaxResidues+1,(*jc).residue,currentParm->ResidueName((*jc).residue));
      mprintf(" %6i:%-4s",(*jc).atom[0],currentParm->AtomName((*jc).atom[0]));
      mprintf(" %6i:%-4s",(*jc).atom[1],currentParm->AtomName((*jc).atom[1]));
      mprintf(" %6i:%-4s",(*jc).atom[2],currentParm->AtomName((*jc).atom[2]));
      mprintf(" %6i:%-4s",(*jc).atom[3],currentParm->AtomName((*jc).atom[3]));
      mprintf(" %6.2lf%6.2lf%6.2lf%6.2lf %i\n",(*jc).C[0],(*jc).C[1],(*jc).C[2],(*jc).C[3],
              (*jc).type);
      MaxResidues++;
    }
  }    
  return 0;  
}

// Jcoupling::action()
/** For each dihedral defined in JcouplingInfo, perform the dihedral and
  * Jcoupling calculation.
  */
int Jcoupling::action() {
  double phi,J;
  double phitemp,C0,C1,C2,C3;
  int residue;

  outputfile.IO->Printf("#Frame %i\n",frameNum+OUTPUTFRAMESHIFT);

  for (std::vector<jcouplingInfo>::iterator jc = JcouplingInfo.begin();
                                            jc !=JcouplingInfo.end();
                                            jc++)
  {
    phi = currentFrame->DIHEDRAL( (*jc).atom[0], (*jc).atom[1], (*jc).atom[2], (*jc).atom[3] );
    C0 = (*jc).C[0];
    C1 = (*jc).C[1];
    C2 = (*jc).C[2];
    C3 = (*jc).C[3];
    if ((*jc).type==1) {
      //J = JcouplingC((*jc).C, phi);
      //phitemp = phi + C3; // Only necessary if offsets become used in perez-type calc
      J = C0 + (C1 * cos(phi)) + (C2 * cos(phi * 2.0)); 
    } else {
      //J = JcouplingABC((*jc).C, phi);
      phitemp = cos( phi + C3 );
      J = (C0 * phitemp * phitemp) + (C1 * phitemp) + C2;
    }

    residue = (*jc).residue;
    // DEBUG - output
    outputfile.IO->Printf("%5i %4s%4s%4s%4s%4s%12lf%12lf\n",
            residue+1,currentParm->ResidueName(residue),
            currentParm->AtomName((*jc).atom[0]),currentParm->AtomName((*jc).atom[1]),
            currentParm->AtomName((*jc).atom[2]),currentParm->AtomName((*jc).atom[3]),
            phi*RADDEG,J);
    //mprintf("%5i %4s",residue+1,P->resnames[residue]);
    //mprintf("%4s",P->names[(*jc).atom[0]]);
    //mprintf("%4s",P->names[(*jc).atom[1]]);
    //mprintf("%4s",P->names[(*jc).atom[2]]);
    //mprintf("%4s",P->names[(*jc).atom[3]]);
    //mprintf("%12lf%12lf\n",phi,J);
  }

  //fprintf(outfile,"%10i %10.4lf\n",frameNum,D);
  
  return 0;
} 
