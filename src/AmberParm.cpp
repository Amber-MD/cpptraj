/* AmberParm.cpp
 * Class that holds parameter information. Can be read in from Amber Topology
 * or PDB files.
 */
#include <cstdlib>
#include <cstring>
#include <cstdio> // For sscanf, sprintf
#include "AmberParm.h" // PtrajFile.h 
#include "ptrajmask.h"
#include "PDBfileRoutines.h"
#include "Mol2FileRoutines.h"
#include "CpptrajStdio.h"

// ================= PRIVATE FUNCTIONS =========================

/*
 * AmberParm::ResName()
 * Given a residue number, set buffer with residue name. Replace blanks with _
 */
void AmberParm::ResName(char *buffer, int res) {
  if (res<0 || res>nres) return;
  strcpy(buffer, resnames[res]);
  if (buffer[3]==' ') buffer[3]='_';
}

/*
 * AmberParm::SetFormat()
 * Given a fortran-type format string, set the format type, data type, the 
 * expected number of columns, width of data, and printf format string.
 * Also calculate the buffer size necessary to read or write 
 * the data in Amber Parm Format given the number of data elements N.
 */
void AmberParm::SetFormat(char *Format, int N) {
  int bufferLines;

  fFormat=UNKNOWN_FFORMAT;
  fType=UNKNOWN_FTYPE;
  width=0;
  numCols=0;
  FormatString="";

  if        ( strncmp(Format,"%FORMAT(10I8)"  ,13)==0 ) {
    fFormat=F10I8;
    fType=FINT;
    width=8;
    numCols=10;
    FormatString="%8i";
  } else if ( strncmp(Format,"%FORMAT(5E16.8)",15)==0 ) {
    fFormat=F5E16_8;
    fType=FDOUBLE;
    width=16;
    numCols=5;
    FormatString="%16.8lE";
  } else if ( strncmp(Format,"%FORMAT(20a4)"  ,13)==0 ) {
    fFormat=F20a4;
    fType=FCHAR;
    width=4;
    numCols=20;
    FormatString="%4s";
  } else if ( strncmp(Format,"%FORMAT(12I6)",  13)==0 ) {
    fFormat=F12I6;
    fType=FINT;
    width=6;
    numCols=12;
    FormatString="%6i";
  } else if ( strncmp(Format,"%FORMAT(3I8)"   ,12)==0 ) {
    fFormat=F3I8;
    fType=FINT;
    width=8;
    numCols=3;
    FormatString="%8i";
  } else {
    rprintf("Error: Unrecognized format in parm file %s: %s\n",parmName,Format);
    return;
  }

  BufferSize= N * width;
  bufferLines = N / numCols;
  if ((N % numCols)!=0) bufferLines++;
  // If DOS file there are CR before Newlines
  if (File.isDos) bufferLines*=2;
  BufferSize+=bufferLines;
  //if (debug>0) 
  //  fprintf(stdout,"*** Buffer size is %i including %i newlines.\n",BufferSize,bufferLines);
}

/* 
 * AmberParm::getFlagFileValues()
 * Search for the FLAG specified by Key and return the values. The values will 
 * be put into an array type according to the FORMAT string in the top file, 
 * but it is necessary to explictly type the returned array. maxval is used to 
 * allocate memory for the return array - only maxval values will be read.
 */
void *AmberParm::getFlagFileValues(const char *Key, int maxval){
  int i; 
  char lineBuffer[BUFFER_SIZE]; // Hold flag/format line from parmfile
  char value[83];      // Hold Key from Flag line
  char temp[17];       // Hold data element
  char *buffer,*ptr;
  char **C;
  int *I;
  double *D;

  if (debug>0) {
    mprintf("Reading %s\n",Key);
    mprintf("DEBUG: maxval= %i\n",maxval);
  }

  C=NULL; D=NULL; I=NULL;
  
  // First, rewind the input file.
  File.IO->Rewind();

  // Next, search for the required FLAG
  while ( File.IO->Gets(lineBuffer,BUFFER_SIZE) == 0) {
    if ( strncmp(lineBuffer,"%FLAG",5)==0 ) {
      sscanf(lineBuffer,"%*s %s",value);
      if (strcmp(value,Key)==0) {
        if (debug>0) mprintf("DEBUG: Found Flag Key [%s]\n",value);
        // Read format line 
        File.IO->Gets(lineBuffer,BUFFER_SIZE);
        if (debug>0) mprintf("DEBUG: Format line [%s]\n",lineBuffer);
        // Set format
        this->SetFormat(lineBuffer,maxval);
        if (debug>0) mprintf("DEBUG: Format type %i\n",fFormat);
        if (fFormat == UNKNOWN_FFORMAT) return NULL;
        // Allocate memory based on data type
        switch (fType) {
          case UNKNOWN_FTYPE : return NULL;
          case FINT   : I=(int*) malloc(maxval*sizeof(int)); break;
          case FDOUBLE: D=(double*) malloc(maxval*sizeof(double)); break;
          case FCHAR  : 
            C=(char**) malloc(maxval*sizeof(char*));
            for (i=0; i<maxval; i++) C[i]=(char*) malloc(5*sizeof(char));
            break;
        }
        // Allocate memory to read in entire section
        buffer=(char*) calloc(BufferSize,sizeof(char));
        if ( File.IO->Read(buffer,sizeof(char),BufferSize)==-1 ) {
          rprintf("ERROR in read of prmtop section %s\n",Key);
          free(buffer);
          break; // Send us outside the while loop
        }
        if (debug>3) mprintf("DEBUG: Buffer [%s]\n",buffer);

        // Convert values in buffer to their type
        ptr=buffer;
        temp[width]='\0';          
        for (i=0; i<maxval; i++) {
          // Advance past newlines - DOS newline is CR
          //fprintf(stdout,"0 %i: %c %i\n",i,ptr[0],ptr[0]);
          while (ptr[0]=='\n' || ptr[0]=='\r') { ptr++; }
          //fprintf(stdout,"1 %i: %c %i\n",i,ptr[0],ptr[0]);
          strncpy(temp,ptr,width);
          if (debug>3) mprintf("DEBUG:   %8i buffer %s\n",i,temp);
          // Sanity check: If we have hit another FLAG before we've read maxval
          // values this is bad.
          // NOTE: Could do addtional type checking too.
          if ( strncmp(ptr,"%FLAG",5)==0 ) {
            rprintf("Error: #values read (%i) < # expected values (%i).\n",i,maxval);
            if (I!=NULL) free(I);
            if (D!=NULL) free(D);
            if (C!=NULL) {for(i=0;i<maxval;i++) free(C[i]); free(C);}
            free(buffer);
            return NULL;
          }
          /* Now convert temp to appropriate type */
          if (I!=NULL) I[i]=atoi(temp);
          if (D!=NULL) D[i]=atof(temp);
          if (C!=NULL) strcpy(C[i],temp);
          ptr+=width;
        }
        /* All done! */
        free(buffer);
        if (I!=NULL) return I;
        if (D!=NULL) return D;
        if (C!=NULL) return C;
      } // End if (strcmp(value,Key)==0)
    } // End if ( strncmp(lineBuffer,"%FLAG",5)==0 )
  } // End While loop
  // If we have scanned through the input file and have not found Key, bad! 
  rprintf("Error: Could not find key %s in file.\n",Key);
  if (I!=NULL) free(I);
  if (D!=NULL) free(D);
  if (C!=NULL) {for(i=0;i<maxval;i++) free(C[i]); free(C);}

  return NULL;
}

// CONSTRUCTOR
AmberParm::AmberParm(int debugIn) {
  debug=debugIn;
  fFormat = UNKNOWN_FFORMAT; numCols=0; width=0; FormatString=""; BufferSize=0;
  fType = UNKNOWN_FTYPE;

  values=NULL; 
  resnums=NULL;  
  names=NULL;  
  resnames=NULL; 
  mass=NULL;    
  charge=NULL;
  bonds=NULL;  
  bondsh=NULL;   
  types=NULL;   
  atomsPerMol=NULL; 
  Box=NULL;    
  pindex=0;      
  parmFrames=0; 
  outFrame=0;
  finalSoluteRes=0; 
  molecules=0;       
  firstSolvMol=0;
  solventMask=NULL; 
  solventMolecules=0; 
  solventAtoms=0;
  solventMoleculeStart=NULL; 
  solventMoleculeStop=NULL;
  natom=0;      
  nres=0; 
  ifbox=0; 
  parmName=NULL; 
  SurfaceInfo=NULL;
}

// DESTRUCTOR
AmberParm::~AmberParm() {
  int i;

  if (names!=NULL) {
    for (i=0; i < natom; i++) free(names[i]);
    free(names);
  }
  if (resnames!=NULL) {
    for (i=0; i < nres; i++) free(resnames[i]);
    free(resnames);
  }
  if (types!=NULL) {
    for (i=0; i < natom; i++) free(types[i]);
    free(types);
  }
  if (mass!=NULL) free(mass);
  if (charge!=NULL) free(charge);
  if (values!=NULL) free(values);
  if (resnums!=NULL) free(resnums);
  if (bonds!=NULL) free(bonds);
  if (bondsh!=NULL) free(bondsh);
  if (atomsPerMol!=NULL) free(atomsPerMol);
  if (Box!=NULL) free(Box);
  if (solventMoleculeStart!=NULL) free(solventMoleculeStart);
  if (solventMoleculeStop!=NULL) free(solventMoleculeStop);
  if (solventMask!=NULL) free(solventMask);
  if (parmName!=NULL) free(parmName);
  if (SurfaceInfo!=NULL) free(SurfaceInfo);
}

/*
 * ---------========= ROUTINES FOR ACCESSING VALUES ARRAY =========---------
 */
int AmberParm::NbondsWithH() { 
  if (values!=NULL) 
    return values[NBONH]; 
  else
    return -1;
}
int AmberParm::NbondsWithoutH() { 
  if (values!=NULL)
    return values[MBONA];
  else
    return -1;
}

/* 
 * AmberParm::OpenParm()
 * Attempt to open file and read in parameters.
 */
int AmberParm::OpenParm(char *filename) {
  if ( File.SetupFile(filename,READ,UNKNOWN_FORMAT, UNKNOWN_TYPE,debug) ) return 1;

  // Copy parm filename to parmName. Separate from File.filename in case of stripped parm
  parmName=(char*) malloc( (strlen(File.basefilename)+1) * sizeof(char));
  strcpy(parmName,File.basefilename);

  if ( File.OpenFile() ) return 1;

  switch (File.fileFormat) {
    case AMBERPARM : if (ReadParmAmber()) return 1; break;
    case PDBFILE   : if (ReadParmPDB()  ) return 1; break;
    case MOL2FILE  : if (ReadParmMol2() ) return 1; break;
    default: 
      rprintf("Unknown parameter file type: %s\n",File.filename);
      return 1;
  }

  File.CloseFile();

  // Create a last dummy residue in resnums that holds natom, which would be
  // the atom number of the next residue if it existed. Shift the number by 1
  // to be consistent with the rest of the array. Do this to be
  // consistent with ptrajmask selection behavior - saves an if-then stmt
  resnums=(int*) realloc(resnums,(nres+1)*sizeof(int));
  resnums[nres]=natom+1;
  // DEBUG
  //fprintf(stdout,"==== DEBUG ==== Resnums for %s:\n",File.filename);
  //for (err=0; err<nres; err++) 
  //  fprintf(stdout,"    %i: %i\n",err,resnums[err]);

  // Set up solvent information
  SetSolventInfo();

  if (debug>0) {
    mprintf("  Number of atoms= %i\n",natom);
    mprintf("  Number of residues= %i\n",nres);
  }

  return 0;
}

/*
 * AmberParm::AssignLCPO()
 * Assign parameters for LCPO method. All radii are incremented by 1.4 Ang.
 */
void AmberParm::AssignLCPO(SurfInfo *S, double vdwradii, double P1, double P2, 
                           double P3, double P4) {
  S->vdwradii = vdwradii + 1.4;
  S->P1 = P1;
  S->P2 = P2;
  S->P3 = P3;
  S->P4 = P4;
}

/*
 * WarnLCPO()
 * Called when the number of bonds to the atom of type atype is not 
 * usual.
 */
void WarnLCPO(char *atype, int atom, int numBonds) {
  mprintf("Warning: Unusual number of bonds for atom %i (%i), type %-2s.\n",
          atom, numBonds, atype);
  mprintf("Using default atom parameters.\n");
}

/* 
 * AmberParm::SetSurfaceInfo()
 * Set up parameters only used in surface area calcs.
 * LCPO method from:
 *   J. Weiser, P.S. Shenkin, and W.C. Still,
 *   "Approximate atomic surfaces from linear combinations of pairwise
 *   overlaps (LCPO)", J. Comp. Chem. 20:217 (1999).
 * Adapted from gbsa=1 method in SANDER, mdread.f
 */
int AmberParm::SetSurfaceInfo() {
  int *numBonds; // # of bonded neighbors each atom has (LCPO only?)
  int i,atom1,atom2;
  char atype[2];
 
  // If surface info already set up exit 
  if (SurfaceInfo!=NULL) return 0;
 
  // If no bond information exit
  if (bonds==NULL) {
    mprintf("Error: SetSurfaceInfo(): Parm %s does not contain bond info.\n",parmName);
    return 1;
  } 

  // If no atom type information exit
  if (types==NULL) {
    mprintf("Error: SetSurfaceInfo(): Parm %s does not contain atom type info.\n",
            parmName);
    return 1;
  }
 
  // Get the number of bonded neighbors for each atom
  numBonds = (int*) malloc(natom * sizeof(int));
  memset(numBonds, 0, natom * sizeof(int));
  // NOTE: Why ignore hydrogens?
  for (i = 0; i < values[MBONA]*3; i+=3) {
     atom1 = bonds[i  ] / 3;
     atom2 = bonds[i+1] / 3;
     numBonds[atom1]++;
     numBonds[atom2]++;
  }

  // DEBUG
  //for (i=0; i<natom; i++)
  //  fprintf(stdout,"DEBUG:    Atom %6i_%4s: %2i bonds.\n",i,names[i],numBonds[i]);

  // Set vdw radii and LCPO parameters
  SurfaceInfo = (SurfInfo*) malloc(natom * sizeof(SurfInfo));
  for (i=0; i < natom; i++) {
    atype[0] = types[i][0];
    atype[1] = types[i][1];

    if (atype[0]=='C' && atype[1]=='T') {
      switch ( numBonds[i] ) {
        case 1: AssignLCPO(SurfaceInfo+i, 1.70, 0.77887, -0.28063, -0.0012968, 0.00039328); break;
        case 2: AssignLCPO(SurfaceInfo+i, 1.70, 0.56482, -0.19608, -0.0010219, 0.0002658);  break;
        case 3: AssignLCPO(SurfaceInfo+i, 1.70, 0.23348, -0.072627, -0.00020079, 0.00007967); break;
        case 4: AssignLCPO(SurfaceInfo+i, 1.70, 0.00000, 0.00000, 0.00000, 0.00000); break;
        default: WarnLCPO(atype,i,numBonds[i]);
                AssignLCPO(SurfaceInfo+i, 1.70, 0.77887, -0.28063, -0.0012968, 0.00039328);
      }
    } else if (atype[0]=='C' || atype[0]=='c') {
      switch ( numBonds[i] ) {
        case 2: AssignLCPO(SurfaceInfo+i, 1.70, 0.51245, -0.15966, -0.00019781, 0.00016392); break;
        case 3: AssignLCPO(SurfaceInfo+i, 1.70, 0.070344, -0.019015, -0.000022009, 0.000016875); break;
        default: WarnLCPO(atype,i,numBonds[i]);
                AssignLCPO(SurfaceInfo+i, 1.70, 0.77887, -0.28063, -0.0012968, 0.00039328);
      }
    } else if (atype[0]=='O' && atype[1]==' ') {
      AssignLCPO(SurfaceInfo+i, 1.60, 0.68563, -0.1868, -0.00135573, 0.00023743);
    } else if (atype[0]=='O' && atype[1]=='2') {
      AssignLCPO(SurfaceInfo+i, 1.60, 0.88857, -0.33421, -0.0018683, 0.00049372);
    } else if (atype[0]=='O' || atype[0]=='o') {
      switch (numBonds[i]) {
        case 1: AssignLCPO(SurfaceInfo+i, 1.60, 0.77914, -0.25262, -0.0016056, 0.00035071); break;
        case 2: AssignLCPO(SurfaceInfo+i, 1.60, 0.49392, -0.16038, -0.00015512, 0.00016453); break;
        default: WarnLCPO(atype,i,numBonds[i]);
                AssignLCPO(SurfaceInfo+i, 1.60, 0.77914, -0.25262, -0.0016056, 0.00035071);
      }
    } else if (atype[0]=='N' && atype[1]=='3') {
      switch (numBonds[i]) {
        case 1: AssignLCPO(SurfaceInfo+i, 1.65, 0.078602, -0.29198, -0.0006537, 0.00036247); break;
        case 2: AssignLCPO(SurfaceInfo+i, 1.65, 0.22599, -0.036648, -0.0012297, 0.000080038); break;
        case 3: AssignLCPO(SurfaceInfo+i, 1.65, 0.051481, -0.012603, -0.00032006, 0.000024774); break;
        default: WarnLCPO(atype,i,numBonds[i]);
                AssignLCPO(SurfaceInfo+i, 1.65, 0.078602, -0.29198, -0.0006537, 0.00036247);
      }
    } else if (atype[0]=='N' || atype[0]=='n') {
      switch (numBonds[i]) {
        case 1: AssignLCPO(SurfaceInfo+i, 1.65, 0.73511, -0.22116, -0.00089148, 0.0002523); break;
        case 2: AssignLCPO(SurfaceInfo+i, 1.65, 0.41102, -0.12254, -0.000075448, 0.00011804); break;
        case 3: AssignLCPO(SurfaceInfo+i, 1.65, 0.062577, -0.017874, -0.00008312, 0.000019849); break;
        default: WarnLCPO(atype,i,numBonds[i]);
                AssignLCPO(SurfaceInfo+i, 1.65, 0.078602, -0.29198, -0.0006537, 0.00036247);
      }
    } else if (atype[0]=='S' && atype[1]=='H') {
      AssignLCPO(SurfaceInfo+i, 1.90, 0.7722, -0.26393, 0.0010629, 0.0002179);
    } else if (atype[0]=='S' || atype[0]=='s') {
      AssignLCPO(SurfaceInfo+i, 1.90, 0.54581, -0.19477, -0.0012873, 0.00029247);
    } else if (atype[0]=='P' || atype[1]=='p') {
      switch (numBonds[i]) {
        case 3: AssignLCPO(SurfaceInfo+i, 1.90, 0.3865, -0.18249, -0.0036598, 0.0004264); break;
        case 4: AssignLCPO(SurfaceInfo+i, 1.90, 0.03873, -0.0089339, 0.0000083582, 0.0000030381); break;
        default: WarnLCPO(atype,i,numBonds[i]);
          AssignLCPO(SurfaceInfo+i, 1.90, 0.3865, -0.18249, -0.0036598, 0.0004264);
      }
    } else if (atype[0]=='Z') {
      AssignLCPO(SurfaceInfo+i, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000);
    } else if (atype[0]=='H' || atype[0]=='h') {
      AssignLCPO(SurfaceInfo+i, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000);
    } else if (atype[0]=='M' && atype[1]=='G') {
      //  Mg radius = 0.99A: ref. 21 in J. Chem. Phys. 1997, 107, 5422
      //  Mg radius = 1.18A: ref. 30 in J. Chem. Phys. 1997, 107, 5422
      //  Mg radius = 1.45A: Aqvist 1992
      //  The following P1-4 values were taken from O.sp3 with two bonded 
      //  neighbors -> O has the smallest van der Waals radius 
      //  compared to all other elements which had been parametrized
      AssignLCPO(SurfaceInfo+i, 1.18, 0.49392, -0.16038, -0.00015512, 0.00016453);
    } else {
      mprintf("Warning: Using carbon SA parms for unknown atom type %i %2s\n",i,atype);
      AssignLCPO(SurfaceInfo+i, 1.70, 0.51245, -0.15966, -0.00019781, 0.00016392);
    }
  } // END LOOP OVER natom

  // DEBUG
  /*
  for (i=0; i<natom; i++) {
    fprintf(stdout,"%6i %4s: %6.2lf %lf %lf %lf %lf\n",i+1,types[i],SurfaceInfo[i].vdwradii,
    fprintf(stdout,"%6i%6.2lf%12.8lf%12.8lf%12.8lf%12.8lf\n",i+1,SurfaceInfo[i].vdwradii,
            SurfaceInfo[i].P1,SurfaceInfo[i].P2,SurfaceInfo[i].P3,SurfaceInfo[i].P4);
  }
  */
  free(numBonds);
  return 0;
}

/*
 * AmberParm::SetSolventInfo()
 * Assuming atomsPerMol has been read in, set solvent information.
 */
int AmberParm::SetSolventInfo() {
  int mol, molAtom, maskAtom; // i, j, k

  // Allocate memory
  solventMask=(char*) malloc(natom * sizeof(char));
  for (maskAtom=0; maskAtom<natom; maskAtom++) solventMask[maskAtom]='F';
  solventMoleculeStart=(int*) malloc(natom * sizeof(int));
  solventMoleculeStop=(int*) malloc(natom * sizeof(int));
  solventMolecules=0;
  solventAtoms=0;

  // Treat all the molecules starting with firstSolvMol (nspsol) as solvent
  if (atomsPerMol!=NULL) {
    molAtom = 0;
    for (mol=0; mol < molecules; mol++) {
      if (mol+1 >= firstSolvMol) {
        // Add this molecule to the solvent list
        solventAtoms += atomsPerMol[mol];
        for (maskAtom=molAtom; maskAtom < molAtom+atomsPerMol[mol]; maskAtom++)
          solventMask[maskAtom] = 'T';
        solventMoleculeStart[solventMolecules] = molAtom;
        solventMoleculeStop[ solventMolecules] = molAtom+atomsPerMol[mol];
        solventMolecules++;
      }
      molAtom += atomsPerMol[mol];
    }

  // Treat all residues named WAT as solvent
  } else if (resnums!=NULL) {
    for (mol=0; mol < nres; mol++) { 
      if ( strcmp("WAT ", resnames[mol]) == 0 ||
           strcmp(" WAT", resnames[mol]) == 0    ) {
        // Add this residue to the list of solvent 
        molAtom = resnums[mol+1] - resnums[mol];
        solventAtoms += molAtom;
        solventMoleculeStart[solventMolecules] = resnums[mol] - 1;
        solventMoleculeStop[ solventMolecules] = resnums[mol+1] - 1;
        solventMolecules++;
        for (maskAtom=resnums[mol] - 1; maskAtom < resnums[mol+1] - 1; maskAtom++)
          solventMask[maskAtom] = 'T';
      }
    }
  }

  // NOTE: Deallocate memory if no solvent atoms?
  if (debug>0)
    mprintf("    %i solvent molecules, %i solvent atoms.\n",
            solventMolecules, solventAtoms);

  return 0; 
}
     
/* 
 * AmberParm::ReadParmAmber() 
 * Read parameters from Amber Topology file
 */
int AmberParm::ReadParmAmber() {
  int err, atom;
  int *solvent_pointer;

  if (debug>0) mprintf("Reading Amber Topology file %s\n",parmName);

  values=(int*) getFlagFileValues("POINTERS",31);
  if (values==NULL) {
    mprintf("Could not get values from topfile\n");
    return 1;
  }

  natom=values[NATOM];
  nres=values[NRES];
  ifbox=values[IFBOX];
  if (debug>0)
    mprintf("    Amber top contains %i atoms, %i residues.\n",natom,nres);

  err=0;
  names=(char**) getFlagFileValues("ATOM_NAME",natom);
  if (names==NULL) {mprintf("Error in atom names.\n"); err++;}
  types=(char**) getFlagFileValues("AMBER_ATOM_TYPE",natom);
  if (types==NULL) {mprintf("Error in atom types.\n"); err++;}
  resnames=(char**) getFlagFileValues("RESIDUE_LABEL",nres);
  if (resnames==NULL) {mprintf("Error in residue names.\n"); err++;}
  resnums=(int*) getFlagFileValues("RESIDUE_POINTER",nres);
  if (resnums==NULL) {mprintf("Error in residue numbers.\n"); err++;}
  mass=(double*) getFlagFileValues("MASS",natom);
  if (mass==NULL) {mprintf("Error in masses.\n"); err++;}
  charge=(double*) getFlagFileValues("CHARGE",natom);
  if (charge==NULL) {mprintf("Error in charges.\n"); err++;}
  // Convert charges to units of electron charge
  for (atom=0; atom < natom; atom++)
    charge[atom] /= 18.2223;
  bonds=(int*) getFlagFileValues("BONDS_WITHOUT_HYDROGEN",values[MBONA]*3);
  if (bonds==NULL) {mprintf("Error in bonds w/o H.\n"); err++;}
  bondsh=(int*) getFlagFileValues("BONDS_INC_HYDROGEN",values[NBONH]*3);
  if (bondsh==NULL) {mprintf("Error in bonds inc H.\n"); err++;}
  if (ifbox>0) {
    solvent_pointer=(int*) getFlagFileValues("SOLVENT_POINTERS",3);
    if (solvent_pointer==NULL) {
      mprintf("Error in solvent pointers.\n"); 
      err++;
    } else {
      finalSoluteRes=solvent_pointer[0];
      molecules=solvent_pointer[1];
      firstSolvMol=solvent_pointer[2];
      free(solvent_pointer);
    }
    atomsPerMol=(int*) getFlagFileValues("ATOMS_PER_MOLECULE",molecules);
    if (atomsPerMol==NULL) {mprintf("Error in atoms per molecule.\n"); err++;}
    Box=(double*) getFlagFileValues("BOX_DIMENSIONS",4);
    if (Box==NULL) {mprintf("Error in Box information.\n"); err++;}
    if (debug>0) {
      mprintf("    %s contains box info: %i mols, first solvent mol is %i\n",
              parmName, molecules, firstSolvMol);
      mprintf("    BOX: %lf %lf %lf %lf\n",Box[0],Box[1],Box[2],Box[3]);
    }
  }

  if ( err>0 ) {
    rprintf("Error reading topfile\n");
    return 1;
  }

  return 0;
}

/* 
 * AmberParm::ReadParmPDB()
 * Open the PDB file specified by filename and set up topology data.
 * Mask selection requires natom, nres, names, resnames, resnums.
 */
int AmberParm::ReadParmPDB() {
  char buffer[256];
  int bufferLen;  
  int currResnum;

  mprintf("    Reading PDB file %s as topology file.\n",parmName);
  currResnum=-1;
  memset(buffer,' ',256);

  while ( File.IO->Gets(buffer,256)==0 ) {
    // If ENDMDL or END is reached stop reading
    if ( strncmp(buffer,"END",3)==0) break;
    // If TER increment number of molecules and continue
    if ( strncmp(buffer,"TER",3)==0) {
      molecules++;
      continue;
    }
    // Skip all other non-ATOM records
    if (strncmp(buffer,"ATOM",4)!=0 &&
        strncmp(buffer,"HETATM",6)!=0 ) continue;

    // Detect and remove trailing newline
    bufferLen = strlen(buffer);
    if (buffer[bufferLen-1] == '\n') buffer[bufferLen-1]='\0';

    // Allocate memory for atom name
    names=(char**) realloc(names, (natom+1) * sizeof(char*));
    names[natom]=pdb_name(buffer);

    // If this residue number is different than the last, allocate mem for new res
    if (currResnum!=pdb_resnum(buffer)) {
      resnames=(char**) realloc(resnames, (nres+1) * sizeof(char*));
      resnames[nres]=pdb_resname(buffer);
      resnums=(int*) realloc(resnums, (nres+1) * sizeof(int));
      resnums[nres]=natom+1; // +1 since in Amber Top atoms start from 1
      currResnum=pdb_resnum(buffer);
      nres++;
    }
    // Clear the buffer
    memset(buffer,' ',256);

    natom++;
  }
  if (debug>0) 
    mprintf("    PDB contains %i atoms, %i residues, %i molecules.\n",
            natom,nres,molecules);
  // If no atoms, probably issue with PDB file
  if (natom<=0) {
    mprintf("Error: No atoms in PDB file.\n");
    return 1;
  }

  return 0;
}

/*
 * AmberParm::ReadParmMol2()
 * Read file as a Tripos Mol2 file.
 */
int AmberParm::ReadParmMol2() {
  char buffer[MOL2BUFFERSIZE];
  int mol2bonds, atom;
  int resnum, currentResnum;
  int numbonds, numbondsh;
  int numbonds3, numbondsh3;
  char resName[5];

  currentResnum=-1;
  mprintf("    Reading Mol2 file %s as topology file.\n",parmName);
  // Get @<TRIPOS>MOLECULE information
  if (Mol2ScanTo(&File, MOLECULE)) return 1;
  //   Scan title
  if ( File.IO->Gets(buffer,MOL2BUFFERSIZE) ) return 1;
  mprintf("      Mol2 Title: [%s]\n",buffer);
  //   Scan # atoms and bonds
  // num_atoms [num_bonds [num_subst [num_feat [num_sets]]]]
  if ( File.IO->Gets(buffer,MOL2BUFFERSIZE) ) return 1;
  mol2bonds=0;
  sscanf(buffer,"%i %i",&natom, &mol2bonds);
  mprintf("      Mol2 #atoms: %i\n",natom);
  mprintf("      Mol2 #bonds: %i\n",mol2bonds);

  // Allocate memory for atom names and types 
  names = (char**) malloc( natom * sizeof(char*));
  types = (char**) malloc( natom * sizeof(char*));
  charge = (double*) malloc( natom * sizeof(double));

  // Get @<TRIPOS>ATOM information
  if (Mol2ScanTo(&File, ATOM)) return 1;
  for (atom=0; atom < natom; atom++) {
    if ( File.IO->Gets(buffer,MOL2BUFFERSIZE) ) return 1;
    names[atom] = (char*) malloc( 5 * sizeof(char));
    types[atom] = (char*) malloc( 5 * sizeof(char));
    // atom_id atom_name x y z atom_type [subst_id [subst_name [charge [status_bit]]]]
    sscanf(buffer,"%*i %s %*f %*f %*f %s %i %s %lf", names[atom], types[atom],
           &resnum,resName, charge+atom);
    //mprintf("      %i %s %s %i %s %lf\n",atom,names[atom],types[atom],resnum,resName,charge[atom]);
    // Check if residue number has changed - if so record it
    if (resnum != currentResnum) {
      resnames = (char**) realloc(resnames, (nres+1) * sizeof(char*));
      resnames[nres] = (char*) malloc( 5 * sizeof(char));
      resnums=(int*) realloc(resnums, (nres+1) * sizeof(int));
      resnums[nres]=atom+1; // +1 since in Amber Top atoms start from 1
      strcpy(resnames[nres], resName);
      currentResnum = resnum;
      nres++;
    }
  }

  // Get @<TRIPOS>BOND information [optional]
  numbonds=0;
  numbondsh=0;
  numbonds3=0;
  numbondsh3=0;
  if (Mol2ScanTo(&File, BOND)==0) {
    for (atom=0; atom < mol2bonds; atom++) {
      if ( File.IO->Gets(buffer,MOL2BUFFERSIZE) ) return 1;
      // bond_id origin_atom_id target_atom_id bond_type [status_bits]
      //         resnum         currentResnum
      sscanf(buffer,"%*i %i %i\n",&resnum,&currentResnum);
      // mol2 atom #s start from 1
      if ( names[resnum-1][0]=='H' || names[currentResnum-1][0]=='H' ) {
        bondsh = (int*) realloc(bondsh, (numbondsh+1)*3*sizeof(int));
        bondsh[numbondsh3  ]=(resnum-1)*3;
        bondsh[numbondsh3+1]=(currentResnum-1)*3;
        bondsh[numbondsh3+2]=0; // Need to assign some force constant eventually
        //mprintf("      Bond to Hydrogen %s-%s %i-%i\n",names[resnum-1],names[currentResnum-1],
        //        resnum, currentResnum);
        numbondsh++;
        numbondsh3+=3;
      } else {
        bonds = (int*) realloc(bonds, (numbonds+1)*3*sizeof(int));
        bonds[numbonds3  ]=(resnum-1)*3;
        bonds[numbonds3+1]=(currentResnum-1)*3;
        bonds[numbonds3+2]=0; // Need to assign some force constant eventually
        //mprintf("      Bond %s-%s %i-%i\n",names[resnum-1],names[currentResnum-1],
        //        resnum, currentResnum);
        numbonds++;
        numbonds3+=3;
      }
    }
    
  } else {
    mprintf("      Mol2 file does not contain bond information.\n");
  }

  // Set up AMBER Pointers array
  values = (int*) calloc(31, sizeof(int));
  values[NATOM] = natom;
  values[NRES] = nres;
  values[NBONH] = numbondsh;
  values[MBONA] = numbonds;

  mprintf("    Mol2 contains %i atoms, %i residues,\n", natom,nres);
  mprintf("    %i bonds to H, %i other bonds.\n", numbondsh,numbonds);

  return 0;
}

/*
 * AmberParm::ParmInfo()
 * Print parm information for atoms in mask.
 */
void AmberParm::ParmInfo(char *maskstr) {
  char *mask;
  int atom,res;
  if (maskstr==NULL) {
    mprintf("Error: AmberParm::ParmInfo: No mask given.\n");
    return;
  }
  mask = parseMaskString(maskstr, natom, nres, names, resnames, resnums, NULL, 'd',debug);
  if (mask==NULL) {
    mprintf("Warning: AmberParm::ParmInfo: Atom mask is NULL.\n");
    return;
  }
  for (atom=0; atom < natom; atom++) {
    if (mask[atom]=='F') continue;
    res = atomToResidue(atom);
    mprintf("Atom %i:%4s Res %i:%4s Mol %i",atom+1,names[atom],res+1,resnames[res],
            atomToMolecule(atom)+1);
    mprintf(" Type=%4s",types[atom]);
    mprintf(" Charge=%lf",charge[atom]);
    mprintf(" Mass=%lf\n",mass[atom]);
  }
  free(mask);
}

/*
 * AmberParm::Info()
 * Print information about this parm to buffer.
 */
void AmberParm::Info(char *buffer) {

  sprintf(buffer,"%i atoms, %i res, box %i, %i mol, %i solvent mol, %i frames",
          natom,nres,ifbox,molecules,solventMolecules,parmFrames);
}
  

/* 
 * AmberParm::mask
 * Wrapper for the old ptraj mask parser from ptrajmask.h
 */
char *AmberParm::mask(char *maskstr) {
  
  // Should never be called with NULL string
  if (maskstr==NULL) return NULL;
  // NOTE: Args 7-8 are for distance criteria selection. Arg 7 is coords,
  // arg 8 should be f for float (NOT IMPLEMENTED) or d for double.
  // Last arg is debug level.
  return parseMaskString(maskstr, natom, nres, names, resnames, resnums, NULL, 'd', debug);
}

/*
 * AmberParm::mask
 * Wrapper for old ptraj mask parser with coordinates
 */
char *AmberParm::mask(char *maskstr, double *X) {
  if (maskstr==NULL) return NULL;
  return parseMaskString(maskstr, natom, nres, names, resnames, resnums, X, 'd', debug);
}

// NOTE: The following atomToX functions do not do any memory checks!
/* 
 * AmberParm::atomToResidue()
 * Given an atom number, return corresponding residue number.
 */
int AmberParm::atomToResidue(int atom) {
  int i, atom1;

  atom1 = atom + 1; // Since in resnums atom numbers start from 1
  for (i = 0; i < nres; i++)
    if ( atom1>=resnums[i] && atom1<resnums[i+1] )
      return i;

  return -1;
}

/*
 * AmberParm::atomToMolecule()
 * Given an atom number, return corresponding molecule number.
 */
int AmberParm::atomToMolecule(int atom) {
  int i, a, atom1;

  atom1 = atom + 1; // Since in atomsPerMol numbers start from 1
  a = 0;
  for (i = 0; i < molecules; i++) {
    a += atomsPerMol[i];
    if (atom1 <= a)
      return i;
  }
  return -1;
}

/*
 * AmberParm::atomToSolventMolecule()
 * Given an atom number, return corresponding solvent molecule
 * NOTE: Could this be achieved with atomToMolecule and solventMask?
 */
int AmberParm::atomToSolventMolecule(int atom) {
  int i, atom1;

  atom1 = atom + 1; 
  for (i = 0; i < molecules; i++) {
    if (atom1 <= solventMoleculeStart[i])
      return -1;
    else if (atom1>solventMoleculeStart[i] && atom1<=solventMoleculeStop[i])
      return i;
  }

  return -1;
}

/*
 * SetupBondArray()
 * Given an atom map and new parm, set up bond array
 * NOTE: Set up atom map to be atom*3??
 */
int *SetupBondArray(int *atomMap, int oldN3, int *oldBonds, int *newN) {
  int *bonds;
  int N3, i, atom1, atom2;

  bonds=NULL;
  N3=0;
  // Go through Bonds with/without H, use atomMap to determine what goes into newParm
  for (i=0; i < oldN3; i+=3) {
    // Check that atom1 and atom2 exist in newParm
    // In the bond arrays atom nums are multiplied by 3
    atom1 = atomMap[ oldBonds[i]/3   ];
    atom2 = atomMap[ oldBonds[i+1]/3 ];
    if ( atom1!=-1 && atom2!=-1 ) {
      // Put new atom 1 and new atom 2 in newParm array
      bonds=(int*) realloc(bonds, (N3+3) * sizeof(int)); 
      bonds[N3]   = atom1 * 3;
      bonds[N3+1] = atom2 * 3;
      bonds[N3+2] = oldBonds[i+2];
      N3+=3;
    }
  }
  
  *newN = N3 / 3;
  return bonds;
}

/*
 * AmberParm::modifyStateByMask()
 * Adapted from ptraj
 *  The goal of this routine is to create a new AmberParm (newParm)
 *  based on the current AmberParm (this), deleting atoms that are
 *  not in the Selected array.
 * NOTE: Make all solvent/box related info dependent on IFBOX only?
 */
AmberParm *AmberParm::modifyStateByMask(int *Selected, int Nselected) {
  AmberParm *newParm;
  int selected;
  int i, ires, imol; 
  int j, jres, jmol;
  int curres, curmol; 
  int *atomMap; // Convert atom # in oldParm to newParm; -1 if atom is not in newParm

  // Allocate space for the new state
  newParm = new AmberParm(debug); 

  // Allocate space for arrays and perform initialization
  newParm->values = (int*) calloc(31, sizeof(int));
  atomMap = (int*) malloc( this->natom * sizeof(int));
  for (i=0; i<this->natom; i++) atomMap[i]=-1;
  newParm->names    = (char**)  malloc( this->natom   * sizeof(char*) );
  for (i=0; i<this->natom; i++) 
    newParm->names[i]=(char*) malloc(5 * sizeof(char));
  newParm->types    = (char**)  malloc( this->natom   * sizeof(char*) );
  for (i=0; i<this->natom; i++)
    newParm->types[i]=(char*) malloc(5 * sizeof(char));
  newParm->charge   = (double*) malloc( this->natom   * sizeof(double));
  newParm->mass     = (double*) malloc( this->natom   * sizeof(double));
  newParm->resnames = (char**)  malloc( this->nres    * sizeof(char*) );
  for (i=0; i<this->nres; i++)
    newParm->resnames[i]=(char*) malloc(5 * sizeof(char));
  newParm->resnums  = (int*)    malloc((this->nres+1) * sizeof(int   ));

  if (this->molecules>0) 
    newParm->atomsPerMol = (int*) malloc(this->molecules * sizeof(int));

  if (this->solventMolecules>0) {
    // Set first solvent molecule to -1 for now
    newParm->firstSolvMol=-1;
  }

  j = 0; 
  jres = -1; jmol = -1;
  ires = -1; imol = -1;

  // Loop over Selected atoms and set up information for the newstate if the atom is 
  // not to be deleted...
  for (selected=0; selected < Nselected; selected++) {
    // i = old atom #, j = new atom number
    i = Selected[selected];          // Atom to be kept from oldParm
    curres = this->atomToResidue(i); // Residue number of atom in oldParm
    atomMap[i]=j;                    // Store this atom in the atom map
    // Copy over atom information
    strcpy(newParm->names[j], this->names[i]);
    strcpy(newParm->types[j], this->types[i]);
    newParm->charge[j]      = this->charge[i];
    newParm->mass[j]        = this->mass[i];

    // Check to see if we are in the same residue or not and copy relevant information
    if (ires == -1 || ires != curres) {
      jres++;
      strcpy(newParm->resnames[jres], this->resnames[curres]);
      newParm->resnums[jres] = j+1;
      ires = curres;
    }

    // Check to see if we are in the same molecule or not and increment #atoms in molecule
    if (this->molecules>0) {
      curmol = this->atomToMolecule(i);
      if (imol == -1 || imol != curmol) {
        jmol++;
        newParm->atomsPerMol[jmol]=1;
        imol = curmol;
      } else {
        newParm->atomsPerMol[jmol]++;
      }
    }

    // If we are keeping this atom and it belongs to a solvent molecule and 
    // the first solvent atom has not been set, set it.
    if (this->solventMolecules>0 && this->solventMask[i]=='T' && newParm->firstSolvMol<0) {
      newParm->firstSolvMol = jmol + 1;
      newParm->finalSoluteRes = jres;
    }

    // Increment the new atom counter
    j++;

  } // End loop over oldParm Selected atoms 

  // Set up bond arrays
  newParm->bondsh = SetupBondArray(atomMap, this->values[NBONH]*3, this->bondsh, 
                                   &(newParm->values[NBONH]));
  newParm->bonds  = SetupBondArray(atomMap, this->values[MBONA]*3, this->bonds, 
                                   &(newParm->values[MBONA]));
  free(atomMap);

  // Fix up IPRES
  newParm->resnums[jres+1] = j+1;

  // Set up new parm information
  newParm->natom = j;
  newParm->nres = jres+1; 
  newParm->ifbox = this->ifbox;
  newParm->parmFrames = this->parmFrames;
  if (this->molecules>0) 
    newParm->molecules = jmol+1;

  // Give stripped parm the same pindex as original parm
  newParm->pindex = this->pindex;
  
  // Reallocate memory 
  newParm->charge=(double*) realloc(newParm->charge, newParm->natom * sizeof(double));
  newParm->mass=(double*) realloc(newParm->mass, newParm->natom * sizeof(double));
  newParm->resnums=(int*) realloc(newParm->resnums, (newParm->nres+1) * sizeof(int));
  for (i=newParm->natom; i<this->natom; i++) free(newParm->names[i]);
  newParm->names=(char**) realloc(newParm->names, newParm->natom * sizeof(char*));
  for (i=newParm->natom; i<this->natom; i++) free(newParm->types[i]);
  newParm->types=(char**) realloc(newParm->types, newParm->natom * sizeof(char*));
  for (i=newParm->nres; i<this->nres; i++) free(newParm->resnames[i]);
  newParm->resnames=(char**) realloc(newParm->resnames, (newParm->nres+1) * sizeof(char*));
  if (newParm->molecules>0)
    newParm->atomsPerMol=(int*) realloc(newParm->atomsPerMol, newParm->molecules * sizeof(int));

  // Set up solvent info if necessary
  if (newParm->firstSolvMol < 0) {
    // No solvent in stripped parmtop
    newParm->solventMolecules=0;
  } else {
    // Set up new solvent info based on new resnums and firstSolvMol
    mprintf("           New parm: First solvent molecule is %i\n",newParm->firstSolvMol);
    newParm->SetSolventInfo();
  }
  
  // Copy box information
  if (this->Box!=NULL) {
    newParm->Box=(double*) malloc(4*sizeof(double));
    for (i=0; i<4; i++)
      newParm->Box[i] = this->Box[i];
  }

  // Set values up
  // NOTE: Eventually set all pointers up?
  newParm->values[NATOM] = newParm->natom;
  newParm->values[NRES] = newParm->nres;
  newParm->values[IFBOX] = newParm->ifbox;

  mprintf("           New parmtop contains %i atoms.\n",newParm->natom);
  mprintf("                                %i residues.\n",newParm->nres);
  mprintf("                                %i molecules.\n",newParm->molecules);
  mprintf("                                %i solvent molcules.\n",newParm->solventMolecules);

  return newParm;
}

/* 
 * AmberParm::DataToBuffer()
 * Return char buffer containing N data elements stored in I, D, or C with 
 * given fortran format.
 * NOTE: Needs error checking
 */
char *AmberParm::DataToBuffer(char *bufferIn, const char *format, 
                              int *I, double *D, char **C, int N) {
  int coord;
  char *ptr, *buffer;

  this->SetFormat((char*)format,N);

  if (bufferIn==NULL) 
    buffer=(char*) malloc( (BufferSize+1) * sizeof(char));
  else
    buffer=(char*) realloc(bufferIn, (BufferSize+1) * sizeof(char));
  
  if (buffer==NULL) return NULL;

  //fprintf(stdout,"*** Called DataToBuffer: N=%i, width=%i, numCols=%i, format=%s\n",
  //        N, width, numCols, FormatString);
  //fprintf(stdout,"*** Buffer address is %p\n",buffer);
  
  ptr=buffer;
  for (coord=0; coord<N; coord++) {
    if      (I!=NULL) sprintf(ptr,FormatString,I[coord]);
    else if (D!=NULL) sprintf(ptr,FormatString,D[coord]);
    else if (C!=NULL) sprintf(ptr,FormatString,C[coord]);
    ptr+=width;
    if ( ((coord+1)%numCols)==0 ) {
      sprintf(ptr,"\n");
      ptr++;
    }
  } 
  // If the coord record didnt end on a newline, print one
  if ( (coord%numCols)!=0 ) {
    sprintf(ptr,"\n");
    //ptr++; // Only needed if more will be written
  } 
  
  //fprintf(stdout,"*** Ptr ended on %p\n",ptr);
  return buffer;
}

/*
 * PrintFlagFormat()
 * Given an Amber Top file FLAG and fortran FORMAT, print to outfile.
 */
void PrintFlagFormat(PtrajFile *outfile, const char *Flag, const char *Format) {
  outfile->IO->Printf("%-80s\n",Flag);
  outfile->IO->Printf("%-80s\n",Format);
}

/*
 * AmberParm::WriteAmberParm()
 * Write out information from current AmberParm to an Amber parm file
 */
int AmberParm::WriteAmberParm() {
  PtrajFile outfile;
  char *buffer;
  int solvent_pointer[3];
  int atom;

  if (parmName==NULL) return 1;

  if ( outfile.SetupFile(parmName, WRITE, AMBERPARM, STANDARD, debug) )
    return 1;

  buffer=NULL;
  if (outfile.OpenFile()) return 1;

  mprintf("    Writing out amber topology file %s.\n",parmName);

  // HEADER AND TITLE - Eventually use actual date and time
  outfile.IO->Printf("%-80s\n","%VERSION  VERSION_STAMP = V0001.000  DATE = 12/03/01  13:16:16");
  PrintFlagFormat(&outfile, "%FLAG TITLE", "%FORMAT(20a4)");
  outfile.IO->Printf("%-80s\n","");

  // POINTERS
  PrintFlagFormat(&outfile, "%FLAG POINTERS", "%FORMAT(10I8)");
  buffer = DataToBuffer(buffer,"%FORMAT(10I8)", values, NULL, NULL, 31);
  outfile.IO->Write(buffer, sizeof(char), BufferSize);

  // ATOM NAMES
  PrintFlagFormat(&outfile, "%FLAG ATOM_NAME", "%FORMAT(20a4)");
  buffer = DataToBuffer(buffer,"%FORMAT(20a4)", NULL, NULL, names, natom);
  outfile.IO->Write(buffer, sizeof(char), BufferSize);

  // CHARGE
  // Convert charges to AMBER charge units
  for (atom=0; atom<natom; atom++)
    charge[atom] *= 18.2223;
  PrintFlagFormat(&outfile, "%FLAG CHARGE", "%FORMAT(5E16.8)");
  buffer = DataToBuffer(buffer,"%FORMAT(5E16.8)", NULL, charge, NULL, natom);
  outfile.IO->Write(buffer, sizeof(char), BufferSize);

  // MASS 
  PrintFlagFormat(&outfile, "%FLAG MASS", "%FORMAT(5E16.8)");
  buffer = DataToBuffer(buffer,"%FORMAT(5E16.8)", NULL, mass, NULL, natom);
  outfile.IO->Write(buffer, sizeof(char), BufferSize);

  // RESIDUE LABEL - resnames
  PrintFlagFormat(&outfile, "%FLAG RESIDUE_LABEL", "%FORMAT(20a4)");
  buffer = DataToBuffer(buffer,"%FORMAT(20a4)", NULL, NULL, resnames, nres);
  outfile.IO->Write(buffer, sizeof(char), BufferSize);

  // RESIDUE POINTER - resnums, IPRES
  PrintFlagFormat(&outfile, "%FLAG RESIDUE_POINTER", "%FORMAT(10I8)");
  buffer = DataToBuffer(buffer,"%FORMAT(10I8)", resnums, NULL, NULL, nres);
  outfile.IO->Write(buffer, sizeof(char), BufferSize);

  // AMBER ATOM TYPE - might be null if read from pdb
  if (types!=NULL) {
    PrintFlagFormat(&outfile, "%FLAG AMBER_ATOM_TYPE", "%FORMAT(20a4)");
    buffer = DataToBuffer(buffer,"%FORMAT(20a4)", NULL, NULL, types, natom);
    outfile.IO->Write(buffer, sizeof(char), BufferSize);
  }

  // BONDS INCLUDING HYDROGEN - might be null if read from pdb
  if (bondsh != NULL) {
    PrintFlagFormat(&outfile, "%FLAG BONDS_INC_HYDROGEN", "%FORMAT(10I8)");
    buffer = DataToBuffer(buffer,"%FORMAT(10I8)", bondsh, NULL, NULL, values[NBONH]*3);
    outfile.IO->Write(buffer, sizeof(char), BufferSize);
  }

  // BONDS WITHOUT HYDROGEN - might be null if read from pdb
  if (bonds!=NULL) {
    PrintFlagFormat(&outfile, "%FLAG BONDS_WITHOUT_HYDROGEN", "%FORMAT(10I8)");
    buffer = DataToBuffer(buffer,"%FORMAT(10I8)", bonds, NULL, NULL, values[MBONA]*3);
    outfile.IO->Write(buffer, sizeof(char), BufferSize);
  }

  // SOLVENT POINTERS
  if (ifbox>0) {
    PrintFlagFormat(&outfile, "%FLAG SOLVENT_POINTERS", "%FORMAT(3I8)");
    solvent_pointer[0]=finalSoluteRes;
    solvent_pointer[1]=molecules;
    solvent_pointer[2]=firstSolvMol;
    buffer = DataToBuffer(buffer,"%FORMAT(10I8)", solvent_pointer, NULL, NULL, 3);
    outfile.IO->Write(buffer, sizeof(char), BufferSize);

    // ATOMS PER MOLECULE
    PrintFlagFormat(&outfile, "%FLAG ATOMS_PER_MOLECULE", "%FORMAT(10I8)");
    buffer = DataToBuffer(buffer,"%FORMAT(10I8)", atomsPerMol, NULL, NULL, molecules);
    outfile.IO->Write(buffer, sizeof(char), BufferSize);

    // BOX DIMENSIONS
    PrintFlagFormat(&outfile, "%FLAG BOX_DIMENSIONS", "%FORMAT(5E16.8)");
    buffer = DataToBuffer(buffer,"%FORMAT(5E16.8)", NULL, Box, NULL, 4);
    outfile.IO->Write(buffer, sizeof(char), BufferSize);
  }

  free(buffer);
  outfile.CloseFile();

  return 0;
}
