#ifndef INC_TIMER_H
#define INC_TIMER_H
class Timer {
  public:
    Timer();
    void Start() { GetWallTime(start_sec_, start_ns_); }
    void Stop()  { GetWallTime(stop_sec_,  stop_ns_);  }
    double Total();
  private:
    void GetWallTime(int&, int&);
    int start_sec_;
    int start_ns_;
    int stop_sec_;
    int stop_ns_;
};
#endif
