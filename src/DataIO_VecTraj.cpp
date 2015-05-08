#include "DataIO_VecTraj.h"
#include "CpptrajStdio.h"
#include "DataSet_Vector.h"
#include "ParmFile.h"
#include "Trajout_Single.h"

// CONSTRUCTOR
DataIO_VecTraj::DataIO_VecTraj() : trajoutFmt_(TrajectoryFile::UNKNOWN_TRAJ) {
  SetValid( DataSet::VECTOR );
}

void DataIO_VecTraj::WriteHelp() {
  mprintf("\t[trajfmt <format>] [parmout <file>]\n");
}

int DataIO_VecTraj::processWriteArgs(ArgList& argIn) {
  trajoutFmt_ = TrajectoryFile::GetFormatFromString( argIn.GetStringKey("trajfmt") );
  parmoutName_ = argIn.GetStringKey("parmout");
  return 0;
}

int DataIO_VecTraj::WriteData(std::string const& fname, DataSetList const& SetList) {
  if (SetList.empty()) return 1;
  // Create pseudo-topology for all vectors.
  Topology pseudo;
  int nres = 1;
  int natom = 0;
  int vec_size = -1;
  for (DataSetList::const_iterator set = SetList.begin(); set != SetList.end(); ++set) {
    if ((*set)->Type() != DataSet::VECTOR)
      mprintf("Warning: Set '%s' is not a vector, skipping.\n", (*set)->legend());
    else {
      DataSet_Vector const& Vec = static_cast<DataSet_Vector const&>( *(*set) );
      if (vec_size == -1)
        vec_size = (int)Vec.Size();
      else if (vec_size != (int)Vec.Size()) {
        mprinterr("Error: Vector '%s' size %zu != first vector size %zu.\n"
                  "Error:   All vectors must have same size.\n",
                  Vec.legend(), Vec.Size(), vec_size);
        return 1;
      }
      pseudo.AddTopAtom(Atom("OXYZ", ' ', 0), nres, "VEC", 0);
      pseudo.AddTopAtom(Atom("VXYZ", ' ', 0), nres, "VEC", 0);
      pseudo.AddBond(natom, natom+1);
      natom += 2;
      ++nres;
    }
  }
  pseudo.CommonSetup(false);
  if (!parmoutName_.empty()) {
    ParmFile pfile;
    if (pfile.WriteTopology( pseudo, parmoutName_, ParmFile::UNKNOWN_PARM, 0 )) {
      mprinterr("Error: Could not write pseudo topology to '%s'\n", parmoutName_.c_str());
      return 1;
    }
  }
  Trajout_Single out;
  if (out.InitTrajWrite(fname, ArgList(), &pseudo, trajoutFmt_) == 0) {
    Frame outFrame(pseudo.Natom());
    for (int i = 0; i != vec_size; ++i) {
      outFrame.ClearAtoms();
      for (DataSetList::const_iterator set = SetList.begin(); set != SetList.end(); ++set) {
        if ((*set)->Type() == DataSet::VECTOR) {
          DataSet_Vector const& Vec = static_cast<DataSet_Vector const&>( *(*set) );
          Vec3 const& ovec = Vec.OXYZ(i);
          outFrame.AddVec3( ovec );
          outFrame.AddVec3( Vec[i] + ovec );
        }
      }
      if (out.WriteFrame(i, &pseudo, outFrame)) return 1;
    }
    out.EndTraj();
  } else {
    mprinterr("Error: Could not set up '%s' for write.\n", fname.c_str());
    return 1;
  }
  return 0;
} 
