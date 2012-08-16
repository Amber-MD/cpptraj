#include <cstdio> //sscanf
#include <cstring> //strncmp
#include <cmath> // sqrt, tanh
#include "ModesInfo.h"
#include "CpptrajStdio.h"
#include "CpptrajFile.h"

const size_t ModesInfo::BUFSIZE_ = 1024;

// CONSTRUCTOR
ModesInfo::ModesInfo() :
  type_(MT_UNKNOWN),
  source_(MS_UNKNOWN),
  navgelem_(0),
  avg_(0),
  nvect_(0),
  nvectelem_(0),
  freq_(0),
  evec_(0)
{
  dType_ = MODES;
}

// CONSTRUCTOR
ModesInfo::ModesInfo(modesType tIn, modesSource sIn, std::string& nameIn) :
  type_(tIn),
  source_(sIn),
  navgelem_(0),
  avg_(0),
  nvect_(0),
  nvectelem_(0),
  freq_(0),
  evec_(0)
{
  // DataSet
  name_ = nameIn;
  dType_ = MODES;
}

// Constants used in the CalcXXX routines.
const double ModesInfo::CONSQ = 2.39805E-3;
const double ModesInfo::TKBC2 = 0.46105E-34;
const double ModesInfo::AVO   = 6.023E23;
const double ModesInfo::CNST  = TKBC2 * AVO;
const double ModesInfo::CMTOA = 1.000E8;
const double ModesInfo::TWOPI = 6.2832;
const double ModesInfo::CONT  = CMTOA / TWOPI;

// ModesInfo::SetNavgElem()
int ModesInfo::SetNavgElem(int mask1tot) {
  // TODO: This is already calcd in MatrixType Nelt.
  if (type_ == MT_DIST || type_ == MT_IDEA || type_ == MT_IRED)
    navgelem_ = mask1tot;
  else if (type_ == MT_DISTCOVAR)
    navgelem_ = mask1tot * (mask1tot - 1) / 2;
  else // CORREL, COVAR, MWCOVAR
    navgelem_ = 3 * mask1tot;
  //mprintf("CDBG: mask1tot=%i ModesInfo::navgelem = %i\n",mask1tot, navgelem_);
  return navgelem_;
}

// DESTRUCTOR
ModesInfo::~ModesInfo() {
  if (avg_!=0) delete[] avg_;
  if (freq_!=0) delete[] freq_;
  if (evec_!=0) delete[] evec_;
}

// ModesInfo::ReadEvecFile()
int ModesInfo::ReadEvecFile(std::string& modesfile, int ibeg, int iend) {
  CpptrajFile infile;
  char buffer[BUFSIZE_];

  if (infile.SetupRead( modesfile.c_str(), 0 )) return 1;
  if (infile.OpenFile()) return 1;
  // Read title line
  if (infile.Gets(buffer, BUFSIZE_)!=0) {
    mprinterr("Error: ReadEvecFile(): error while reading title (%s)\n",infile.FullFileStr());
    return 1;
  }
  source_ = MS_FILE;

  std::string title(buffer);
  if (title.find("DISTCOVAR")!=std::string::npos)
    type_ = MT_DISTCOVAR;
  else if (title.find("MWCOVAR")!=std::string::npos)
    type_ = MT_MWCOVAR;
  else if (title.find("COVAR")!=std::string::npos)
    type_ = MT_COVAR;
  else if (title.find("DIST")!=std::string::npos)
    type_ = MT_DIST;
  else if (title.find("CORREL")!=std::string::npos)
    type_ = MT_CORREL;
  else if (title.find("IDEA")!=std::string::npos)
    type_ = MT_IDEA;
  else if (title.find("IRED")!=std::string::npos)
    type_ = MT_IRED;
  else {
    // For compatibility with quasih and nmode output
    mprintf("Warning: ReadEvecFile(): Unrecognized type [%s]\n", title.c_str());
    mprintf("         Assuming MWCOVAR.\n");
    type_ = MT_MWCOVAR;
  }

  // Read number of coords for avg and evec
  if ( infile.Gets(buffer, BUFSIZE_)!=0) {
    mprinterr("Error: ReadEvecFile(): error while reading number of atoms (%s)\n",
              infile.FullFileStr());
    return 1;
  }
  int nvals = sscanf(buffer, "%i %i", &navgelem_, &nvectelem_);
  if (nvals == 0) {
    mprinterr("Error: ReadEvecFile(): sscanf on coords failed (%s)\n",infile.FullFileStr());
    return 1;
  } else if (nvals == 1) {
    mprintf("Warning: ReadEvecFile(): No value for nvectelem found in %s,\n", 
            infile.FullFileStr());
    mprintf("         assuming it is navgelem (%i)\n",navgelem_);
    nvectelem_ = navgelem_;
  }

  // Allocate memory for avg, freq, evec
  if (avg_!=0) delete[] avg_;
  if (freq_!=0) delete[] freq_;
  if (evec_!=0) delete[] evec_;
  avg_ = new double[ navgelem_ ];
  // TODO: Check bounds
  freq_ = new double[ iend - ibeg + 1 ];
  evec_ = new double[ (iend - ibeg + 1) * nvectelem_ ];

  // Read avg coordinates
  int ncoords = navgelem_;
  int nlines = ncoords / 7;
  if (ncoords > 0 && (ncoords % 7) != 0)
    ++nlines;
  int nent = 0;
  double tmpval[7];
  for (int i = 0; i < nlines; ++i) {
    if (infile.Gets(buffer, BUFSIZE_)!=0) {
      mprinterr("Error: ReadEvecFile(): error while reading avg coords (%s)\n",
                infile.FullFileStr());
      return 1;
    }
    nvals = sscanf(buffer, "%lf %lf %lf %lf %lf %lf %lf", tmpval, tmpval+1, tmpval+2,
                   tmpval+3, tmpval+4, tmpval+5, tmpval+6);
    for (int j = 0; j < nvals; j++)
      avg_[nent++] = tmpval[j];
  }
  //mprintf("CDBG: ReadEvecFile(): nent=%i\n", nent);

  // Read eigenvectors
  ncoords = nvectelem_;
  nlines = ncoords / 7;
  if (ncoords > 0 && (ncoords % 7) != 0)
    ++nlines;
  nvect_ = 0;
  int nno = 0;
  while ( infile.Gets(buffer, BUFSIZE_)==0 ) { // This should read in ' ****'
    if (strncmp(buffer," ****", 5)!=0) {
      mprinterr("Error: ReadEvecFile(): When reading eigenvector %i, expected ' ****',\n",nvect_);
      mprinterr("       got %s [%s]\n", buffer, infile.FullFileStr());
      return 1;
    }
    // Read number and freq
    if (infile.Gets(buffer, BUFSIZE_)!=0) {
      mprinterr("Error: ReadEvecFile(): error while reading number and freq (%s)\n",
                infile.FullFileStr());
      return 1;
    }
    if (sscanf(buffer, "%i%lf", &nno, tmpval) != 2) {
      mprinterr("Error: ReadEvecFile(): error while scanning number and freq (%s)\n",
                infile.FullFileStr());
      return 1;
    }
    //mprintf("CDBG:\tVec[%i]: #%i frequency=%lf\n", nvect_, nno, tmpval[0]);
    if (nno >= ibeg && nno <= iend)
      freq_[nvect_] = tmpval[0];
    else if (nno > iend)
      break;
    // Read coords
    nent = 0;
    for (int i = 0; i < nlines; ++i) {
      if (infile.Gets(buffer, BUFSIZE_)!=0) {
        mprinterr("Error: ReadEvecFile(): error while reading evec coords (%s)\n",
                  infile.FullFileStr());
        return 1;
      }
      nvals = sscanf(buffer, "%lf %lf %lf %lf %lf %lf %lf", tmpval, tmpval+1, tmpval+2,
                     tmpval+3, tmpval+4, tmpval+5, tmpval+6);
      for (int j = 0; j < nvals; j++) {
        evec_[ncoords * nvect_ + nent] = tmpval[j];
        ++nent;
      }
    }
    if (nno >= ibeg && nno <= iend)
      ++nvect_;
  }

  if (nvect_ != (iend - ibeg + 1)) {
    mprintf("Warning: Number of read evecs is %i, number of requested evecs is %i\n",
            nvect_, iend - ibeg + 1);
  }
 
  infile.CloseFile();
  // DEBUG
  //mprintf("CDBG: ModesInfo read from %s:",modesfile.c_str());
  //mprintf("navgelem=%i nvect=%i nvectelem=%i\n", navgelem_, nvect_, nvectelem_);
  return 0;
}

/** Calc spectral density (JACS 2002, 124, 4522; eq. A24) */
double ModesInfo::calc_spectral_density(double *taum, int i, double omega) {
  double J = 0.0;
  for(int j = 0 ; j < nvect_; j++){
    J += (freq_[j] * (evec_[j * nvectelem_ + i] * evec_[j * nvectelem_ + i])) * 2.0 * taum[j] /
         ( 1.0 + omega*omega * taum[j]*taum[j] ); //check order XXXX
  }                                                                                                  
  return J;
}

/** Calc rms atomic fluctuations */
double* ModesInfo::CalcRMSfluct(int ibeg, int iend, bool ibose) {
  int natoms = navgelem_ / 3;
  double* results = new double[ natoms * 4 ];
  memset(results, 0, natoms*4*sizeof(double));

  for (int i = 0; i < natoms; ++i) {
    double sumx = 0;
    double sumy = 0;
    double sumz = 0;
    for (int j = 0; j < nvect_; ++j) {
      if (j+1 >= ibeg && j+1 <= iend) { // TODO: Pre calc the shift
        double frq = freq_[j];
        if (frq >= 0.5) {
          // Don't use eigenvectors associated with zero or negative eigenvalues
          int idx = nvectelem_*j + 3*i;
          double distx = evec_[idx  ] * evec_[idx  ];
          double disty = evec_[idx+1] * evec_[idx+1];
          double distz = evec_[idx+2] * evec_[idx+2];
          double fi = 1.0 / (frq*frq);
          if (ibose) {
            double argq = CONSQ * frq;
            fi *= (argq / tanh(argq));
          }
          sumx += distx * fi;
          sumy += disty * fi;
          sumz += distz * fi;
        }
      }
    }

    sumx *= CNST;
    sumy *= CNST;
    sumz *= CNST;

    int ind = 4*i;
    results[ind  ] = sqrt(sumx) * CONT;
    results[ind+1] = sqrt(sumy) * CONT;
    results[ind+2] = sqrt(sumz) * CONT;
    results[ind+3] = sqrt(sumx + sumy + sumz) * CONT;
  }

  return results;
}

/** Calc displacement of coordinates along normal mode directions */
double* ModesInfo::CalcDisplacement(int ibeg, int iend, bool ibose, double factor) {
  int natoms = navgelem_ / 3;
  double* results = new double[ navgelem_ ];
  memset(results, 0, navgelem_*sizeof(double));
  double sqrtcnst = sqrt(CNST) * CONT * factor;

  for (int i = 0; i < nvect_; ++i) {
    if (i+1 >= ibeg && i+1 <= iend) { // TODO: Pre calc the shift
      double frq = freq_[i];
      if (frq >= 0.5) {
        // Don't use eigenvectors associated with zero or negative eigenvalues
        double fi = 1.0 / frq;
        if (ibose) {
          double argq = CONSQ * frq;
          fi *= (fi * argq / tanh(argq));
          fi = sqrt(fi);
        }
        fi *= sqrtcnst; // * CONT * factor
    
        int nvi = nvectelem_ * i;
        for (int j = 0; j < natoms; ++j) {
          int ind = j * 3;
          results[ind  ] += evec_[nvi + ind  ] * fi;
          results[ind+1] += evec_[nvi + ind+1] * fi;
          results[ind+2] += evec_[nvi + ind+2] * fi;
        }
      }
    }
  }

  return results;
}

/** Calc dipole-dipole correlation functions */
double* ModesInfo::CalcDipoleCorr(int ibeg, int iend, bool ibose, 
                                  modestackType const& Apairs)
{
  double vec[3], del[3][3];
  double qcorr = 1.0; // For when ibose is false.
  int rsize = Apairs.size() * (iend - ibeg + 1);
  double* results = new double[ rsize ];
  memset(results, 0, rsize * sizeof(double));

  int ncnt = 0; // Index into results
  for (modestack_it apair = Apairs.begin(); apair != Apairs.end(); ++apair)
  {
    // Convert atom #s to coordinate #s
    int at1 = (*apair).first * 3;
    int at2 = (*apair).second * 3;
    // Calc unit vector along at2->at1 bond
    double dnorm = 0.0;
    for (int i = 0; i < 3; ++i) {
      vec[i] = avg_[at1+i] - avg_[at2+i];
      dnorm += (vec[i] * vec[i]);
    }
    dnorm = sqrt(dnorm);
    vec[0] /= dnorm;
    vec[1] /= dnorm;
    vec[2] /= dnorm;
    // Loop over desired modes
    for (int i = 0; i < nvect_; ++i) {
      if (i+1 >= ibeg && i+1 <= iend) { // TODO: Pre calc the shift
        if (freq_[i] >= 0.5) {
          // Don't use eigenvectors associated with zero or negative eigenvalues
          // NOTE: Where is 11791.79 from? Should it be a const?
          double frq = freq_[i] * freq_[i] / 11791.79;
          if (ibose) {
            double argq = CONSQ * freq_[i];
            qcorr = argq / tanh(argq);
          }
          /*
           *  Calc the correlation matrix for delta
           *    as in eq. 7.16 of lamm and szabo, J Chem Phys 1986, 85, 7334.
           *  Note that the rhs of this eq. should be multiplied by kT
           */
          double qcorrf = qcorr / frq;
          int nvi = nvectelem_ * i;
          for (int j = 0; j < 3; ++j) {
            int idxi1 = at1 + j;
            int idxi2 = at2 + j;
            for (int k = 0; k < 3; ++k) {
              int idxj1 = at1 + k;
              int idxj2 = at2 + k;
              del[j][k] = 0.6 * qcorrf * (evec_[nvi + idxi1] - evec_[nvi + idxi2])
                                       * (evec_[nvi + idxj1] - evec_[nvi + idxj2]);
            }
          }
          // Correlation in length, eq. 10.2 of lamm and szabo
          // NOTE: Commented out in PTRAJ
          /*****
          rtr0 = 0.0;
          for(j = 0; j < 3; j++)
            for(k = 0; k < 3; k++)
              rtr0 += e[j] * e[k] * del[j][k];
          *****/
          // Librational correlation function, using eq. 7.12 of lamm and szabo
          // (w/o beta on the lhs).
          double val = 0.0;
          for (int j = 0; j < 3; ++j) {
            val -= del[j][j];
            for (int k = 0; k < 3; ++k) 
              val += vec[j] * vec[k] * del[j][k];
          }
          val *= (3.0 / (dnorm * dnorm));

          results[ncnt] = val;
          ncnt++;
        } // END if positive definite eigenvalue
      } // END if in between ibeg and iend
    } // END loop over nvect
  } // END loop over atom pairs

  return results;
}

// ModesInfo::ProjectCovar()
void ModesInfo::ProjectCovar(CpptrajFile& outfile, Frame& currentFrame, 
                             AtomMask& mask, std::vector<double> const& sqrtMasses)
{
  double XYZ[3];
  int idx1 = 0;
  for (int i = 0; i < nvect_; ++i) {
    double proj = 0;
    int idx2 = 0;
    std::vector<double>::const_iterator sqrtmass = sqrtMasses.begin();
    for (AtomMask::const_iterator atom = mask.begin(); atom != mask.end(); ++atom)
    {
      currentFrame.GetAtomXYZ( XYZ, *atom );
      double mass = *sqrtmass;
      proj += (XYZ[0] - avg_[idx2  ]) * mass * evec_[idx1  ];
      proj += (XYZ[1] - avg_[idx2+1]) * mass * evec_[idx1+1];
      proj += (XYZ[2] - avg_[idx2+2]) * mass * evec_[idx1+2];
      idx1 += 3;
      idx2 += 3;
      ++sqrtmass;
    }
    // Output projection
    outfile.Printf(" %9.3f", proj);
  }
  outfile.Printf("\n");
}

// ModesInfo::ProjectIDEA()
void ModesInfo::ProjectIDEA(CpptrajFile& outfile, Frame& currentFrame,
                             AtomMask& mask)
{
  double XYZ[3];
  int idx1 = 0;
  for (int i = 0; i < nvect_; ++i) {
    double proj1 = 0;
    double proj2 = 0;
    double proj3 = 0;
    for (AtomMask::const_iterator atom = mask.begin(); atom != mask.end(); ++atom)
    {
      proj1 += XYZ[0] * evec_[idx1];
      proj2 += XYZ[1] * evec_[idx1];
      proj3 += XYZ[2] * evec_[idx1];
      ++idx1;
    }
    // Output projection
    outfile.Printf(" %9.3f %9.3f %9.3f %9.3f", proj1, proj2, proj3,
                   sqrt(proj1*proj1 + proj2*proj2 + proj3*proj3) );
  }
  outfile.Printf("\n");
}

