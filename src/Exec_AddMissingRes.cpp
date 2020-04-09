#include "Exec_AddMissingRes.h"
#include "BufferedLine.h"
#include "CharMask.h"
#include "CpptrajStdio.h"
#include "DataSet_Coords_CRD.h"
#include "ParmFile.h"
#include "StringRoutines.h"
#include "Trajin_Single.h"
#include "Trajout_Single.h"
#include <cstdlib> // atoi
#include <cstring> //strncmp

/** Get gap info from PDB */
int Exec_AddMissingRes::FindGaps(Garray& Gaps, CpptrajFile& outfile, std::string const& pdbname)
const
{
  BufferedLine infile;
  if (infile.OpenFileRead( pdbname )) {
    mprinterr("Error: Could not open '%s' for reading.\n", pdbname.c_str());
    return 1;
  }
  const char* linePtr = infile.Line();
  int inMissing = 0;
  int nmissing = 0;
  bool firstGap = true;
  bool atTheEnd = false;
  std::string lastname, lastchain;
  int lastres = 0;
  std::string currentname, currentchain;
  int currentres;
  while (linePtr != 0) {
    if (strncmp(linePtr, "REMARK", 6) == 0)
    {
      ArgList line(linePtr);
      if (line.Nargs() > 2) {
        if (inMissing == 0) {
          // MISSING section not yet encountered.
          if (line[0] == "REMARK" && line[2] == "MISSING") {
            inMissing = 1;
          }
        } else if (inMissing == 1) {
          // In MISSING, looking for start of missing residues
          if (line[0] == "REMARK" && line[2] == "M") {
            inMissing = 2;
            nmissing = 0;
            firstGap = true;
            atTheEnd = false;
          }
        } else if (inMissing == 2) {
          // Reading MISSING residues
          if (line[1] != "465") {
            atTheEnd = true;
          } else
            nmissing++;
          if (firstGap) {
            // The very first gap TODO check Gaps.empty
            Gaps.push_back( Gap(line[2], atoi(line[4].c_str()), line[3]) );
            lastname = Gaps.back().LastName();
            lastres = Gaps.back().StartRes();
            lastchain = Gaps.back().Chain();
            firstGap = false;
          } else {
            currentname = line[2];
            currentres = atoi(line[4].c_str());
            currentchain = line[3];
            if (atTheEnd || currentres - lastres > 1 || currentchain != lastchain) {
              // New sequence starting or end. Finish current.
              Gaps.back().SetStopRes(lastres);
              /*
              mprintf("  Gap %c %4s %6i to %4s %6i %6u\n",
                      Gaps.back().Chain(),
                      Gaps.back().FirstName().c_str(), Gaps.back().StartRes(),
                      Gaps.back().LastName().c_str(), Gaps.back().StopRes(),
                      Gaps.back().Nres());*/
              if (atTheEnd) {
                break;
              }
              Gaps.push_back( Gap(currentres, currentchain) );
            }
            // Continue the current sequence
            Gaps.back().AddGapRes(currentname);
            lastname = currentname;
            lastres = currentres;
            lastchain = currentchain;
          }
        } // END inMissing == 2
      } // END nargs > 2
    } // END REMARK
    linePtr = infile.Line();
  } // END while linePtr != 0

  // Printout
  for (Garray::const_iterator it = Gaps.begin(); it != Gaps.end(); ++it) {
    outfile.Printf("  Gap %c %4s %6i to %4s %6i %6u\n",
                   it->Chain(),
                   it->FirstName().c_str(), it->StartRes(),
                   it->LastName().c_str(), it->StopRes(),
                   it->Nres());
    // Print residues
    unsigned int col = 1;
    for (Gap::name_iterator name = it->nameBegin(); name != it->nameEnd(); ++name) {
      outfile.Printf("%c", Residue::ConvertResName(*name));
      col++;
      if (col > 80) {
        outfile.Printf("\n");
        col = 1;
      }
    }
    if (col > 1)
      outfile.Printf("\n");
  }
  outfile.Printf("%i missing residues.\n", nmissing);
  if (Gaps.empty()) {
    mprintf("Warning: No gaps found.\n");
  }
  return 0;
} 

/** Try to minimize using steepest descent. */
int Exec_AddMissingRes::Minimize(Topology const& topIn, Frame& frameIn, CharMask const& maskIn)
const
{
  double min_tol = 1.0E-5;
  int max_iteration = 1000;

  // Output trajectory
  int iteration = 0;
  Trajout_Single trajOut;
  if (trajOut.InitTrajWrite("min.nc", ArgList(), DataSetList(), TrajectoryFile::AMBERNETCDF))
    return 1;
  if (trajOut.SetupTrajWrite((Topology*)&topIn, CoordinateInfo(), 0))
    return 1;
  if (trajOut.WriteSingle(iteration, frameIn)) return 1;

  // Selected bonds
  BondArray activeBonds;
  for (BondArray::const_iterator bnd = topIn.Bonds().begin(); bnd != topIn.Bonds().end(); ++bnd)
  {
    if (maskIn.AtomInCharMask( bnd->A1() ) ||
        maskIn.AtomInCharMask( bnd->A2() ))
    {
      mprintf("DEBUG: Bond %i to %i\n", bnd->A1()+1, bnd->A2()+1);
      activeBonds.push_back( *bnd );
    }
  }
  
  // Forces
  std::vector<Vec3> Farray(topIn.Natom(), Vec3(0.0));
  // Coordinates
  //std::vector<Vec3> Xarray;
  //Xarray.reserve( topIn.Natom() );
  //for (int at = 0; at < topIn.Natom(); at++)
  //  Xarray.push_back( Vec3(frameIn.XYZ(at)) );
  // Degrees of freedom
  double deg_of_freedom = 3 * maskIn.Nselected();
  double fnq = sqrt(deg_of_freedom);
  // Main loop for steepest descent
  //const double Rk = 1.0;
  const double dxstm = 1.0E-5;
  const double crits = 1.0E-6;
  double rms = 1.0;
  double dxst = 1.0;
  double last_e = 0.0;
  mprintf("          \t%8s %12s %12s\n", " ", "ENE", "RMS");
  while (rms > min_tol && iteration < max_iteration) {
    double e_total = 0.0;
    // Determine bond energy and forces
    for (BondArray::const_iterator bnd = activeBonds.begin(); bnd != activeBonds.end(); ++bnd)
    {
      BondParmType BP = topIn.BondParm()[ bnd->Idx() ];
      //Vec3 const& XYZ0 = Xarray[ bnd->A1() ];
      //Vec3 const& XYZ1 = Xarray[ bnd->A2() ];
      const double* XYZ0 = frameIn.XYZ( bnd->A1() );
      const double* XYZ1 = frameIn.XYZ( bnd->A2() );
      double rx = XYZ0[0] - XYZ1[0];
      double ry = XYZ0[1] - XYZ1[1];
      double rz = XYZ0[2] - XYZ1[2];
      double r2 = rx*rx + ry*ry + rz*rz;
      if (r2 > 0.0) {
        double r2inv = 1.0/r2;
        double r = sqrt(r2);
        //mprintf("DBG: %i A1=%i A2=%i R=%g\n", iteration, bnd->A1()+1, bnd->A2()+1, r);
        double rinv = r * r2inv;

        double db = r - BP.Req();
        double df = BP.Rk() * db;
        double e = df * db;
        e_total += e;

        df *= 2.0 * rinv;

        double dfx = df * rx;
        double dfy = df * ry;
        double dfz = df * rz;

        if (maskIn.AtomInCharMask(bnd->A1())) {
          Farray[bnd->A1()][0] -= dfx;
          Farray[bnd->A1()][1] -= dfy;
          Farray[bnd->A1()][2] -= dfz;
        }

        if (maskIn.AtomInCharMask(bnd->A2())) {
          Farray[bnd->A2()][0] += dfx;
          Farray[bnd->A2()][1] += dfy;
          Farray[bnd->A2()][2] += dfz;
        }
      }
    }

/*
    unsigned int idx = 0; // Index into FrameDistances
    for (unsigned int f1 = 0; f1 != nframes; f1++)
    {
      for (unsigned int f2 = f1 + 1; f2 != nframes; f2++)
      {
        //double Req = FrameDistances().GetCdist(f1, f2);
        Vec3 V1_2 = Xarray[f1] - Xarray[f2];
        double r2 = V1_2.Magnitude2();
        double s = sqrt(r2);
        double r = 2.0 / s;
        double db = s - FrameDistances().GetElement(idx++);
        double df = Rk * db;
        double e = df * db;
        e_total += e;
        df *= r;
        // Apply force
        V1_2 *= df;
        Farray[f1] -= V1_2;
        Farray[f2] += V1_2;
      }
    }
*/
    // Calculate the magnitude of the force vector.
    double sum = 0.0;
    for (std::vector<Vec3>::const_iterator FV = Farray.begin(); FV != Farray.end(); ++FV)
      sum += FV->Magnitude2();
    rms = sqrt( sum ) / fnq;
    // Adjust search step size
    if (dxst < crits) dxst = dxstm;
    dxst = dxst / 2.0;
    if (e_total < last_e) dxst = dxst * 2.4;
    double dxsth = dxst / sqrt( sum );
    last_e = e_total;
    // Update positions and reset force array.
    double* Xptr = frameIn.xAddress();
    for (int idx = 0; idx != topIn.Natom(); idx++, Xptr += 3)
    //std::vector<Vec3>::iterator FV = Farray.begin();
    //for (std::vector<Vec3>::iterator XV = Xarray.begin();
    //                                 XV != Xarray.end(); ++XV, ++FV)
    {
      //mprintf("xyz0= %g %g %g  Fxyz= %g %g %g\n", Xptr[0], Xptr[1], Xptr[2], Farray[idx][0], Farray[idx][1], Farray[idx][2]);
      Xptr[0] += Farray[idx][0] * dxsth;
      Xptr[1] += Farray[idx][1] * dxsth;
      Xptr[2] += Farray[idx][2] * dxsth;
      Farray[idx] = 0.0;
      //mprintf("xyz1= %g %g %g\n", Xptr[0], Xptr[1], Xptr[2]);
      //*XV += (*FV * dxsth);
      //*FV = 0.0;
    }
    // Write out current E.
    mprintf("Iteration:\t%8i %12.4E %12.4E\n", iteration, e_total, rms);
    iteration++;
    if (trajOut.WriteSingle(iteration, frameIn)) return 1;
  }
  // RMS error
  double sumdiff2 = 0.0;
  for (BondArray::const_iterator bnd = activeBonds.begin(); bnd != activeBonds.end(); ++bnd)
  {
    BondParmType BP = topIn.BondParm()[ bnd->Idx() ];
    //Vec3 const& XYZ0 = Xarray[ bnd->A1() ];
    //Vec3 const& XYZ1 = Xarray[ bnd->A2() ];
    Vec3 XYZ0( frameIn.XYZ( bnd->A1() ) );
    Vec3 XYZ1( frameIn.XYZ( bnd->A2() ) );
    Vec3 V1_2 = XYZ0 - XYZ1;
    double r1_2 = sqrt( V1_2.Magnitude2() );
    double diff = r1_2 - BP.Req();
    sumdiff2 += (diff * diff);
    //if (debug_ > 0)
        mprintf("\t\t%u to %u: D= %g  Eq= %g  Delta= %g\n",
                bnd->A1()+1, bnd->A2()+1, r1_2, BP.Req(), fabs(diff));
  }
/* 
  unsigned int idx = 0; // Index into FrameDistances
  double sumdiff2 = 0.0;
  for (unsigned int f1 = 0; f1 != nframes; f1++)
  {
    for (unsigned int f2 = f1 + 1; f2 != nframes; f2++)
    {
      Vec3 V1_2 = Xarray[f1] - Xarray[f2];
      double r1_2 = sqrt( V1_2.Magnitude2() );
      double Req = FrameDistances().GetElement(idx);
      double diff = r1_2 - Req;
      sumdiff2 += (diff * diff);
      if (debug_ > 0)
        mprintf("\t\t%u to %u: D= %g  Eq= %g  Delta= %g\n",
                f1+1, f2+1, r1_2, Req, fabs(diff));
      ++idx;
    }
  }
*/
  double rms_err = sqrt( sumdiff2 / (double)activeBonds.size() );
  mprintf("\tRMS error of final graph positions: %g\n", rms_err);
/*
  // Write out final graph with cluster numbers.
  std::vector<int> Nums;
  Nums.reserve( nframes );
  if (cnumvtime != 0) {
    ClusterSieve::SievedFrames const& sievedFrames = FrameDistances().FramesToCluster();
    DataSet_1D const& CVT = static_cast<DataSet_1D const&>( *cnumvtime );
    for (unsigned int n = 0; n != nframes; n++)
      Nums.push_back( (int)CVT.Dval(sievedFrames[n]) );
  } else
    for (int n = 1; n <= (int)nframes; n++)
      Nums.push_back( n );
*/
  return 0;
}

/** Write topology and frame to pdb. */
int Exec_AddMissingRes::WriteStructure(std::string const& fname, Topology* newTop, Frame const& newFrame,
                                       TrajectoryFile::TrajFormatType typeOut)
const
{
  Trajout_Single trajOut;
  if (trajOut.InitTrajWrite(fname, ArgList(), DataSetList(), typeOut))
    return 1;
  if (trajOut.SetupTrajWrite(newTop, CoordinateInfo(), 1))
    return 1;
  if (trajOut.WriteSingle(0, newFrame)) return 1;
  trajOut.EndTraj();
  return 0;
}


/// Placeholder for Residues
class Pres {
  public:
    Pres() : oresnum_(0), tresnum_(-1), chain_(' ') {}
    /// CONSTRUCTOR - Take Residue
    Pres(Residue const& res, int resnum) :
      name_(res.Name()), oresnum_(res.OriginalResNum()), tresnum_(resnum), chain_(res.ChainID())
      {}
    /// CONSTRUCTOR - Take name, number, chain
    Pres(std::string const& name, int rnum, char chain) :
      name_(name), oresnum_(rnum), tresnum_(-1), chain_(chain)
      {}
    /// First sort by chain, then by original residue number
    bool operator<(const Pres& rhs) const {
      if (chain_ == rhs.chain_)
        return (oresnum_ < rhs.oresnum_);
      else
        return (chain_ < rhs.chain_);
    }

    NameType const& Name() const { return name_; }
    int OriginalResNum()   const { return oresnum_; }
    int TopResNum()        const { return tresnum_; }
    char ChainID()         const { return chain_; }
  private:
    NameType name_;
    int oresnum_;   ///< Original (PDB) residue number.
    int tresnum_;   ///< Topology residue index; -1 if it was missing.
    char chain_;    ///< Original (PDB) chain ID.
};

/** Try to add in missing residues. */
int Exec_AddMissingRes::AddMissingResidues(DataSet_Coords_CRD* dataOut,
                                           Topology const& topIn,
                                           Frame const& frameIn,
                                           Garray const& Gaps)
{
  typedef std::set<Pres> Pset;
  Pset AllResidues;
  // First add all existing residues
  for (int rnum = 0; rnum < topIn.Nres(); ++rnum) {
    std::pair<Pset::iterator, bool> ret = AllResidues.insert( Pres(topIn.Res(rnum), rnum) );
    if (!ret.second) {
      mprinterr("Internal Error: Somehow residue %s was duplicated.\n",
                topIn.TruncResNameNum(rnum).c_str());
      return 1;
    }
  }

  // Loop over gaps
  for (Garray::const_iterator gap = Gaps.begin(); gap != Gaps.end(); ++gap)
  {
    mprintf("\tGap %c %i to %i\n", gap->Chain(), gap->StartRes(), gap->StopRes());
    int currentRes = gap->StartRes();
    for (Gap::name_iterator it = gap->nameBegin(); it != gap->nameEnd(); ++it, ++currentRes) {
      mprintf("DEBUG: %s %i\n", it->c_str(), currentRes);
      std::pair<Pset::iterator, bool> ret = AllResidues.insert( Pres(*it, currentRes, gap->Chain()) );
      if (!ret.second) {
        mprinterr("Internal Error: Somehow residue %s %i in chain %c was duplicated.\n",
                  it->c_str(), currentRes, gap->Chain());
        return 1;
      }
    }
    /*
    // Start res connector mask
    std::string maskStr0("::" + std::string(1,gap->Chain()) + "&:;" + integerToString(gap->StartRes()-1));
    // Stop res connector mask
    std::string maskStr1("::" + std::string(1,gap->Chain()) + "&:;" + integerToString(gap->StopRes()+1));
    mprintf("\t  Mask0=[%s] Mask1=[%s]\n", maskStr0.c_str(), maskStr1.c_str());
    // Find start res connector
    AtomMask mask0( maskStr0 );
    if (topIn.SetupIntegerMask( mask0, coordsIn )) return 1;
    // Find stop res connector
    AtomMask mask1( maskStr1 );
    if (topIn.SetupIntegerMask( mask1, coordsIn )) return 1;
    mprintf("\t  Selected0=%i Selected1=%i\n", mask0.Nselected(), mask1.Nselected());
    */
  }

  // Print residues.
  // Count number of present atoms and missing residues.
  int nAtomsPresent = 0;
  int nResMissing = 0;
  for (Pset::const_iterator it = AllResidues.begin(); it != AllResidues.end(); ++it)
  {
    mprintf("\t  %6s %8i %8i %c\n", *(it->Name()), it->OriginalResNum(), it->TopResNum()+1, it->ChainID());
    if (it->TopResNum() < 0)
      nResMissing++;
    else
      nAtomsPresent += topIn.Res(it->TopResNum()).NumAtoms();
  }
  mprintf("\t%i atoms present, %i residues missing.\n", nAtomsPresent, nResMissing);

  // Create new Frame
  Frame newFrame(nAtomsPresent + nResMissing);
  newFrame.ClearAtoms();
  // Create a new topology with all residues. For missing residues, create a CA atom.
  Topology* newTop = dataOut->TopPtr();
  // Zero coord for new CA atoms
  Vec3 Zero(0.0);
  // Topology for CA atoms
  Topology CAtop;
  // Frame for CA atoms
  Frame CAframe;
  // Mask for missing CA atoms
  CharMask CAmissing;
  for (Pset::const_iterator it = AllResidues.begin(); it != AllResidues.end(); ++it)
  {
    int topResNum = it->TopResNum();
    if (topResNum < 0) {
      // This was a missing residue
      newTop->AddTopAtom( Atom("CA", "C "),
                          Residue(it->Name(), it->OriginalResNum(), ' ', it->ChainID()) );
      newFrame.AddVec3( Zero );
      // CA top
      CAtop.AddTopAtom( Atom("CA", "C "),
                        Residue(it->Name(), it->OriginalResNum(), ' ', it->ChainID()) );
      CAframe.AddVec3( Zero );
      CAmissing.AddAtom(true);
    } else {
      Residue const& topres = topIn.Res(it->TopResNum());
      Residue newres(topres.Name(), topres.OriginalResNum(), topres.Icode(), topres.ChainID());
      int caidx = -1;
      for (int at = topres.FirstAtom(); at < topres.LastAtom(); at++) {
        if (topIn[at].Name() == "CA") caidx = at;
        newTop->AddTopAtom( Atom(topIn[at].Name(), topIn[at].ElementName()), newres );
        frameIn.printAtomCoord(at);
        newFrame.AddXYZ( frameIn.XYZ(at) );
        newFrame.printAtomCoord(newTop->Natom()-1);
      }
      // CA top
      if (caidx == -1) {
        mprinterr("Error: No CA atom found for residue %s\n", topIn.TruncResNameNum(it->TopResNum()).c_str());
        return 1;
      }
      CAtop.AddTopAtom( Atom(topIn[caidx].Name(), topIn[caidx].ElementName()), newres );
      CAframe.AddXYZ( frameIn.XYZ(caidx) );
      CAmissing.AddAtom(false);
    }
  } // END loop over all residues
  newTop->SetParmName("newpdb", "temp.pdb");
  newTop->CommonSetup( false ); // No molecule search
  newTop->Summary();
  if (WriteStructure("temp.pdb", newTop, newFrame, TrajectoryFile::PDBFILE)) {
    mprinterr("Error: Write of temp.pdb failed.\n");
    return 1;
  }
  // CA top
  // Add pseudo bonds between adjacent CA atoms in the same chain.
  BondParmType CAbond(1.0, 3.8);
  for (int cares = 1; cares < CAtop.Nres(); cares++) {
    // Since only CA atoms, residue # is atom #
    Residue const& res0 = CAtop.Res(cares - 1);
    Residue const& res1 = CAtop.Res(cares);
    if (res0.ChainID() == res1.ChainID()) {
      CAtop.AddBond(cares-1, cares, CAbond);
    }
  }
  CAtop.SetParmName("capdb", "temp.ca.mol2");
  CAtop.CommonSetup(true); // molecule search
  CAtop.Summary();
  if (WriteStructure("temp.ca.mol2", &CAtop, CAframe, TrajectoryFile::MOL2FILE)) {
    mprinterr("Error: Write of temp.ca.mol2 failed.\n");
    return 1;
  }

  if (Minimize(CAtop, CAframe, CAmissing)) {
    mprinterr("Error: Minimization of CA atoms failed.\n");
    return 1;
  }
                           
  return 0;
}

// Exec_AddMissingRes::Help()
void Exec_AddMissingRes::Help() const
{
  mprintf("\tpdbname <pdbname> name <setname> [out <filename>]\n"
          "\t[parmargs <parm args>] [trajargs <trajin args>]\n");
}

// Exec_AddMissingRes::Execute()
Exec::RetType Exec_AddMissingRes::Execute(CpptrajState& State, ArgList& argIn)
{
  std::string pdbname = argIn.GetStringKey("pdbname");
  if (pdbname.empty()) {
    mprinterr("Error: provide PDB name.\n");
    return CpptrajState::ERR;
  }
  mprintf("\tPDB name: %s\n", pdbname.c_str());
  CpptrajFile* outfile = State.DFL().AddCpptrajFile(argIn.GetStringKey("out"),
                                                    "AddMissingRes", DataFileList::TEXT, true);
  if (outfile==0) {
    mprinterr("Internal Error: Unable to allocate 'out' file.\n");
    return CpptrajState::ERR;
  }
  mprintf("\tOutput file: %s\n", outfile->Filename().full());
  ArgList parmArgs;
  std::string parmArgStr = argIn.GetStringKey("parmargs");
  if (!parmArgStr.empty()) {
    parmArgs.SetList(parmArgStr, ",");
    mprintf("\tParm args: %s\n", parmArgStr.c_str());
  }
  ArgList trajArgs;
  std::string trajArgStr = argIn.GetStringKey("trajargs");
  if (!trajArgStr.empty()) {
    trajArgs.SetList(trajArgStr, ",");
    mprintf("\tTraj args: %s\n", trajArgStr.c_str());
  }
  std::string dsname = argIn.GetStringKey("name");
  if (dsname.empty())  {
    mprinterr("Error: Output set name must be specified with 'name'.\n");
    return CpptrajState::ERR;
  }
  DataSet_Coords_CRD* dataOut = (DataSet_Coords_CRD*)State.DSL().AddSet(DataSet::COORDS, dsname);
  if (dataOut == 0) {
    mprinterr("Error: Unable to allocate output coords data set.\n");
    return CpptrajState::ERR;
  }
  mprintf("\tOutput set: %s\n", dataOut->legend());

  // Find missing residues/gaps in the PDB
  Garray Gaps;
  if (FindGaps(Gaps, *outfile, pdbname)) {
    mprinterr("Error: Finding missing residues failed.\n");
    return CpptrajState::ERR;
  }
  mprintf("\tThere are %zu gaps in the PDB.\n", Gaps.size());

  // Read in topology
  ParmFile parmIn;
  Topology topIn;
  if (parmIn.ReadTopology(topIn, pdbname, parmArgs, State.Debug())) {
    mprinterr("Error: Read of topology from PDB failed.\n");
    return CpptrajState::ERR;
  }
  topIn.Summary();

  // Set up input trajectory
  Trajin_Single trajIn;
  if (trajIn.SetupTrajRead(pdbname, trajArgs, &topIn)) {
    mprinterr("Error: Setup of PDB for coordinates read failed.\n");
    return CpptrajState::ERR;
  }
  trajIn.PrintInfo(1);
  // Create input frame
  Frame frameIn;
  frameIn.SetupFrameV(topIn.Atoms(), trajIn.TrajCoordInfo());
  // Read input
  if (trajIn.BeginTraj()) {
    mprinterr("Error: Opening PDB for coordinates read failed.\n");
    return CpptrajState::ERR;
  }
  trajIn.GetNextFrame(frameIn);
  trajIn.EndTraj();

  // Try to add in missing residues
  if (AddMissingResidues(dataOut, topIn, frameIn, Gaps)) {
    mprinterr("Error: Attempt to add missing residues failed.\n");
    return CpptrajState::ERR;
  }

  return CpptrajState::OK;
}
