#ifndef INC_ARRAYITERATOR_H
#define INC_ARRAYITERATOR_H
#include <iterator>
/// Template for iterator over C-type array.
template <class T>
class ArrayIterator : public std::iterator<std::forward_iterator_tag, T> {
  public:
    ArrayIterator() : ptr_(0), increment_(1) {}
    ArrayIterator(const ArrayIterator& rhs) : ptr_(rhs.ptr_), increment_(1) {}
    ArrayIterator(T* pin) : ptr_(pin), increment_(1) {}
    ArrayIterator& operator=(const ArrayIterator& rhs) {
      if (this == &rhs) return *this;
      ptr_ = rhs.ptr_;
      increment_ = rhs.increment_;
      return *this;
    }
    // Relations
    bool operator==(const ArrayIterator& rhs) const { return (ptr_==rhs.ptr_);}
    bool operator!=(const ArrayIterator& rhs) const { return (ptr_!=rhs.ptr_);}
    // Increment
    ArrayIterator& operator++() {
      ptr_ += increment_;
      return *this;
    }
    ArrayIterator operator++(int) {
      ArrayIterator tmp(*this);
      ++(*this);
      return tmp;
    }
    // Value
    T& operator*()  { return *ptr_; }
    // Address
    T* operator->() { return ptr_; }
    // Addition
    ArrayIterator& operator+=(int offset) {
      ptr_ += (offset * increment_);
      return *this;
    }
    ArrayIterator operator+(int offset) {
      ArrayIterator tmp(*this);
      tmp += offset;
      return tmp;
    }
  private:
    T* ptr_;
    int increment_; ///< Size of iteration increment
};
#endif
