#include "Traj_SQM.h"
#include "CpptrajStdio.h"

// Traj_SQM::processWriteArgs()
int Traj_SQM::processWriteArgs(ArgList& argIn) {
  return 0;
}

int Traj_SQM::setupTrajout(std::string const& fname, Topology* trajParm,
                           int NframesToWrite, bool append)
{
  if (trajParm==0) return 1;
  if (append) {
    mprinterr("Error: Append not supported for SQM.\n");
    return 1;
  }
  if (outfile_.SetupWrite( fname, debug_ )) return 1;
  sqmParm_ = trajParm;
  if (NframesToWrite==1) singleWrite_ = true;
  // Set up title
  std::string outTitle = Title();
  if (outTitle.empty()) {
    outTitle.assign("Cpptraj generated SQM input");
  } else {
    if ( outTitle.size() > 80) {
      mprintf("Warning: Amber SQM title for '%s' too long: truncating.\n[%s]\n",
              outfile_.Filename().base(), outTitle.c_str());
      outTitle.resize(80);
    }
  }
  SetTitle( outTitle );

  return 0;
}

// Traj_SQM::writeFrame()
int Traj_SQM::writeFrame(int set, Frame const& frameOut) {
  // If just writing 1 frame dont modify output filename
  if (singleWrite_) {
    if (outfile_.OpenFile()) return 1;
  } else {
    if (outfile_.OpenWriteNumbered( set + 1 )) return 1;
  }
  const Topology& parm = static_cast<const Topology&>( *sqmParm_ );
  for (int atom = 0; atom < parm.Natom(); atom++) {
    const double* XYZ = frameOut.XYZ( atom );
    outfile_.Printf("%2d %-4s %12.7f %12.7f %12.7f\n", parm[atom].AtomicNumber(),
                    parm[atom].c_str(), XYZ[0], XYZ[1], XYZ[2]);
  }
  outfile_.CloseFile();
  return 0;
}

// Traj_SQM::Info()
void Traj_SQM::Info() {
  mprintf("is an SQM input file");
}
