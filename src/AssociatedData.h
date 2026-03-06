#ifndef INC_ASSOCIATEDDATA_H
#define INC_ASSOCIATEDDATA_H
class ArgList;
class AssociatedData {
  public:
    /// Destructor. Virtual since this class is inherited.
    virtual ~AssociatedData() {}
    /// Associated data types
    enum AssociatedType { NOE = 0, CONNECT, RESID };
    /// CONSTRUCTOR - take associated data type
    AssociatedData(AssociatedType t) : type_(t) {}
    /// \return Associated data type
    AssociatedType Type() { return type_; }
    /// \return A copy of the associated data
    virtual AssociatedData* Copy() const = 0;
    /// Print associated data info to stdout
    virtual void Ainfo() const = 0;
    /// Process arguments for associated data
    virtual int ProcessAdataArgs(ArgList&) = 0;
  private:
    AssociatedType type_;
};
#endif
