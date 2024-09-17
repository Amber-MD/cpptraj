#ifndef INC_KERNEL_LJSWITCH_H
#define INC_KERNEL_LJSWITCH_H
/** Switching function for Lennard-Jones. */
template <class T> class Kernel_LJswitch {
  public:
    /// CONSTRUCTOR
    Kernel_LJswitch() : cut2_0_(0), cut2_1_(0) {}
    /// Setup
    void setup_switch(T const& cut2_0in, T const& cut2_1in) {
      cut2_0_ = cut2_0in;
      cut2_1_ = cut2_1in;
    }
    /// Calculate switched value
    T switch_fn(T const& rij2)
    {
      if (rij2 <= cut2_0_)
        return 1.0;
      else if (rij2 > cut2_1_)
        return 0.0;
      else {
        T xoff_m_x = cut2_1_ - rij2;
        T fac = 1.0 / (cut2_1_ - cut2_0_);
        return (xoff_m_x*xoff_m_x) * (cut2_1_ + 2.0*rij2 - 3.0*cut2_0_) * (fac*fac*fac);
      }
    }
  private:
    T cut2_0_;
    T cut2_1_;
};
#endif
