#ifndef INC_DATABLOCK_H
/// An abstraction used for copying one or more elements of a DataSet
class DataBlock {
  public:
    DataBlock() : nelts_(0) {}
    DataBlock(void* pin, unsigned int n) : ptr0_(pin), nelts_(n) {}

    const void* Ptr0()   const { return ptr0_; }
    unsigned int Nelts() const { return nelts_; }
  private:
    void* ptr0_;         ///< Pointer to memory address to copy
    unsigned int nelts_; ///< Number of elements to copy
};
#endif
