#ifndef INC_ASSOCIATEDDATA_H
#define INC_ASSOCIATEDDATA_H
class AssociatedData {
  public:
    /// Destructor. Virtual since this class is inherited.
    virtual ~AssociatedData() {}
    enum AssociatedType { NOE = 0 };
    AssociatedData(AssociatedType t) : type_(t) {}
    AssociatedType Type() { return type_; }
    virtual AssociatedData* Copy() const = 0;
    virtual void Ainfo() const = 0;
  private:
    AssociatedType type_;
};
#endif
