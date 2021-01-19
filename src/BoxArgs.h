#ifndef INC_BOXARGS_H
#define INC_BOXARGS_H
class ArgList;
class Box;
/// Hold XYZ ABG arguments for setting up a box.
class BoxArgs {
  public:
    BoxArgs();
    /// Keywords for setting XYZ ABG
    static const char* Keywords_XyzAbg();
    /// Keywords for setting truncated octohedron
    static const char* Keywords_TruncOct();
    /// Set box arguments from keywords
    int SetBoxArgs(ArgList&);
    /// Set all lengths to given value
    void SetLengths(double);
    /// Set all angles to given value
    void SetAngles(double);
    /// Set info not already set from given box
    int SetMissingInfo(Box const&);
    /// Print info that has been set to STDOUT
    void PrintXyzAbg() const;
    /// \return Pointer to XYZ ABG array
    const double* XyzAbg() const { return xyzabg_; }
  private:
    static int SetEmptyInfo(const char*, double&, const char*, double, const char*, double);

    double xyzabg_[6]; ///< Hold box information to be set.
    bool setVar_[6];   ///< If true, that part of the xyzabg_ array has been set.
};
#endif
