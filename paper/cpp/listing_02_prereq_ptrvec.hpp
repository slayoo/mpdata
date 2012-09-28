#include <boost/ptr_container/ptr_vector.hpp>
template <class T>
struct vec_t : boost::ptr_vector<T>
{
  const T &operator[](const int i) const
  {   
    return this->at(
      (i + this->size()) % this->size()
    ); 
  }
};
