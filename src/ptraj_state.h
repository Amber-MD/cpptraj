#ifndef INC_PTRAJ_STATE_H
#define INC_PTRAJ_STATE_H
#include <stdio.h> // FILE
#include "Name.h" // NAME in ptrajState

// GLOBAL PTRAJ DEBUG LEVEL
// Defined here since both Action_PtrajAction and Action_PtrajAnalysis include
// this file
#ifdef PTRAJ_STATE_MODULE
int prnlev = 0;
#else
extern int prnlev;
#endif
void SetPrnlev(int);

// ---------- ptrajMode --------------------------------------------------------
// originally from ptraj_local.h
typedef enum _ptrajMode {
  PTRAJ_NOOP,  PTRAJ_ACTION, PTRAJ_FIRSTPASS, PTRAJ_SECONDPASS,
  PTRAJ_SETUP, PTRAJ_STATUS, PTRAJ_PRINT,     PTRAJ_CLEANUP
} ptrajMode;


// ---------- Ptraj State ------------------------------------------------------
// originally from ptraj_local.h
// Redefine Name to match definition in Ptraj
typedef NAME Name;
typedef struct _ptrajState {
  double box[6];             // box lengths and angles 
  double *masses;            // atom masses 
  double *charges;           // atom charges 
  int atoms;                 // number of atoms 
  int residues;              // number of residues 
  int *ipres;                // atoms in each residue (atom #s start from 1)
  int *ipres_mask;           // ipres for mask parser, atom #s start from 0 
  int IFBOX;                 // is there box information? 
  int boxfixed[6];           // equal to 1 if this box coordinate is fixed 
  int molecules;             // total number of molecules 
  int *moleculeInfo;         // number of atoms in each molecule 
  int *solventMask;          // atoms in the solvent 
  int solventMolecules;      // number of solvent molecules 
  int *solventMoleculeStart; // pointer into solventMask for first atom of each solvent 
  int *solventMoleculeStop;  // pointer into solventMask for last atom of each solvent 
  int solventAtoms;          // number of solvent atoms 
  Name *atomName;            // atom names 
  Name *residueName;         // residue names 
  int maxFrames;             // number of useful frames in 1 pass of trajin's 
  double temp0;              // DAN TEST: for writing out netcdf temp trajs 
} ptrajState;

#define INITIALIZE_ptrajState( _p_ ) \
  _p_->box[0] = 0.0; _p_->box[1] = 0.0; _p_->box[2] = 0.0; \
  _p_->box[3] = 90.0; _p_->box[4] = 90.0; _p_->box[5] = 90.0; \
  _p_->masses = NULL; \
  _p_->charges = NULL; \
  _p_->atoms = 0; \
  _p_->residues = 0; \
  _p_->ipres = NULL; \
  _p_->IFBOX = 0; \
  _p_->boxfixed[0] = 0; _p_->boxfixed[1] = 0; _p_->boxfixed[2] = 0; \
  _p_->boxfixed[3] = 0; _p_->boxfixed[4] = 0; _p_->boxfixed[5] = 0; \
  _p_->molecules = 0; \
  _p_->moleculeInfo = NULL; \
  _p_->solventMask = NULL; \
  _p_->solventMolecules = 0; \
  _p_->solventMoleculeStart = NULL; \
  _p_->solventMoleculeStop = NULL; \
  _p_->solventAtoms = 0; \
  _p_->atomName = NULL; \
  _p_->residueName = NULL; \
  _p_->maxFrames = 0; \
  _p_->temp0 = 0.0

ptrajState *ptrajCopyState(ptrajState **);
int atomToResidue(int, int, int *);

int isActiveDetailed(int, int, int *, int, int, Name *, Name *, int *);
void printAtomMaskDetailed(FILE *, int *, int, int, Name *, Name *, int *);
void printAtomMask(FILE *, int *, ptrajState *);
void printAtomCompact2(FILE *, int, ptrajState *);
int *processAtomMask( char *, ptrajState * );


#endif
