#ifndef INC_TIMER_H
#define INC_TIMER_H
class Timer {
  public:
    Timer();
    void Start() { GetWallTime(start_sec_, start_ns_); }
    void Stop()  {
      int stop_sec, stop_ns;
      GetWallTime(stop_sec,  stop_ns);
      double seconds = (double)( stop_sec - start_sec_ );
      double nano = (double)( stop_ns - start_ns_);
      total_ += ( seconds + (nano / 1000000000) );
    }
    double Total() const { return total_; }
  private:
    void GetWallTime(int&, int&);
    int start_sec_;
    int start_ns_;
    double total_;
};
#endif
