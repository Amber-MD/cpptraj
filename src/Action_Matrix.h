#ifndef INC_ACTION_MATRIX_H
#define INC_ACTION_MATRIX_H
#include "Action.h"
#include "DataSet_MatrixDbl.h"
#include "DataSet_Vector.h"
#include "Array1D.h"
#include "ActionFrameCounter.h"
/// Calculate various types of matrices
class Action_Matrix : public Action, ActionFrameCounter {
  public:
    Action_Matrix();
    DispatchObject* Alloc() const { return (DispatchObject*)new Action_Matrix(); }
    void Help() const;
  private:
    Action::RetType Init(ArgList&, ActionInit&, int);
    Action::RetType Setup(ActionSetup&);
    Action::RetType DoAction(int, ActionFrame&);
#   ifdef MPI
    int SyncAction();
    Parallel::Comm trajComm_;
#   endif
    void Print();

    typedef DataSet_MatrixDbl::Darray Darray;     ///< Mass/vector array type.
    typedef Darray::iterator          v_iterator; ///< Iterator over vector.
    typedef Darray::const_iterator    M_iterator; ///< Iterator over mass.

    DataSet_MatrixDbl* Mat_;      ///< Matrix used in calculations.
    DataSet_MatrixDbl* matByRes_; ///< Matrix averaged by residue. 
    DataFile* outfile_;           ///< DataFile for writing out matrix.
    CpptrajFile* byMaskOut_;      ///< For outputting matrix averaged over mask(s).
    AtomMask mask1_;
    AtomMask mask2_;
    enum OutputType { BYATOM=0, BYRESIDUE, BYMASK };
    OutputType outtype_;
    int debug_;
    // IRED only
    int order_;                                ///< Legendre order
    std::vector<DataSet_Vector*> IredVectors_; ///< IRED vectors
    // DIHCOVAR only
    Array1D DihedralSets_;
    // MWcovar only
    Darray mass1_; ///< Atom masses corresponding to mask1_.
    Darray mass2_; ///< Atom masses corresponding to mask2_.

    Darray vect2_; ///< Hold diagonal elements squared.
#   ifdef _OPENMP
    /// For OPENMP only, save coord indices (X-Y) for speed 
    std::vector<int> crd_indices_;
#   endif
    bool useMask2_;
    bool useMass_;

    // For ByResidue output.
    typedef std::vector<int> Iarray;
    struct matrix_res {
      Iarray maskIdxs_; ///< Residue atom indices into mask, matrix row/col.
      int resnum_;      ///< Topology residue number.
    };
    typedef std::vector<matrix_res> MatResArray;
    MatResArray residues1_;
    MatResArray residues2_;
    MatResArray MaskToMatResArray(Topology const&, AtomMask const&) const;

    Darray FillMassArray(Topology const&, AtomMask const&) const;
    void CalcIredMatrix(int);
    void CalcDistanceMatrix(Frame const&);
    inline void StoreVec(v_iterator&, v_iterator&, const double*) const;
    void CalcCovarianceMatrix(Frame const&);
    void StoreXY(v_iterator&, v_iterator&, const double*) const;
    void CalcDihedralCovariance( int );
    void CalcIdeaMatrix(Frame const&);
    void CalcCorrelationMatrix(Frame const&);
    void CalcDistanceCovarianceMatrix(Frame const&);
    void Vect2MinusVect();
    void FinishCovariance(size_t);
    inline void DotProdAndNorm(DataSet_MatrixDbl::iterator&, v_iterator&,
                               v_iterator&, v_iterator&, v_iterator&) const;
    void FinishCorrelation();
    void FinishDistanceCovariance();
    double ByMaskAverage(unsigned int, unsigned int) const;
};
#endif
