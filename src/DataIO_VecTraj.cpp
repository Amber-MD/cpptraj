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

int DataIO_VecTraj::WriteData(FileName const& fname, DataSetList const& SetList) {
  if (SetList.empty()) return 1;
  // Create pseudo-topology for all vectors.
  Topology pseudo;
  BondArray bonds;
  // 1 pseudo bond type, Rk = 0.0, Req = 1.0 Ang.
  BondParmArray bParm(1, BondParmType(0.0, 1.0));
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
      Residue vec_res("VEC", nres, ' ', ' ');
      pseudo.AddTopAtom(Atom("OXYZ", 0), vec_res);
      pseudo.AddTopAtom(Atom("VXYZ", 0), vec_res);
      bonds.push_back( BondType(natom, natom+1, 0) ); // Bond parm index 0
      natom += 2;
      ++nres;
    }
  }
  pseudo.SetBondInfo( bonds, BondArray(), bParm );
  pseudo.CommonSetup();
  if (!parmoutName_.empty()) {
    ParmFile pfile;
    if (pfile.WriteTopology( pseudo, parmoutName_, ParmFile::UNKNOWN_PARM, 0 )) {
      mprinterr("Error: Could not write pseudo topology to '%s'\n", parmoutName_.c_str());
      return 1;
    }
  }
  Trajout_Single out;
  if (out.PrepareTrajWrite(fname, ArgList(), &pseudo, CoordinateInfo(),
                           vec_size, trajoutFmt_) == 0)
  {
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
      if (out.WriteSingle(i, outFrame)) return 1;
    }
    out.EndTraj();
  } else {
    mprinterr("Error: Could not set up '%s' for write.\n", fname.full());
    return 1;
  }
  return 0;
} 
