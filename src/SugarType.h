#ifndef INC_SUGARTYPE_H
#define INC_SUGARTYPE_H
/// Hold sugar information
class SugarType {
  public:
    /** The alpha form has the CH2OH substituent (C5-C6 etc in Glycam) on the 
      * opposite side of the OH on the anomeric carbon (C1 in Glycam), while
      * in the beta form it is on the same side. Since the OH is not always 
      * there, cpptraj will look relative to the hydrogen on the anomeric 
      * carbon. So the alpha form has CH2OH on same side as H1 (torsion
      * between -90 and 90), beta form has it on the opposite side.
      */
    enum FormType { UNKNOWN=0, ALPHA, BETA };

    SugarType() : form_(UNKNOWN), resChar_(' '), linkage_(-1) {}

    FormType Form() const { return form_; }
    char ResChar()  const { return resChar_; }
    int Linkage()   const { return linkage_; }
  private:
    FormType form_;
    char resChar_;
    int linkage_;
};
#endif
