
//
// ALGH::Stack<T> class template
//
// Jaeil Choi
// Last modified in Sep, 2006
//
// An simple implementation of classic stack with
//   typical member functions -- push(), read(), pop(), resize()
// Note that it does not extend/shrink the array automatically.
//

#ifndef ALGH_STACK_HPP
#define ALGH_STACK_HPP

#include <iostream>

namespace ALGH {
  

template <class T_t>
class Stack {
private:
  T_t *array;
  int array_size;
  int entry_count;
public:
  Stack() : array(NULL), array_size(0), entry_count(0) {}
  Stack(int size) : array(NULL), array_size(0), entry_count(0) { resize(size); }
  ~Stack() { if (array) free(array); }
  void clear(void) { if (array) free(array); array = NULL; array_size = entry_count = 0; }
  void clearAndFree(void) {
    for (int i=0; array && i<entry_count; i++) if (array[i]) free(array[i]);
    clear();
  }
  inline int count()   { return entry_count; }
  inline int isEmpty() { return (entry_count==0); }
  inline int isFull()  { return (entry_count>=array_size); }
  
  bool resize(int size) {
    if (size < entry_count) return false;
    if (size == array_size) return true;
    array_size = size;
    array = (T_t*) realloc( array, array_size * sizeof(T_t) );
    return (array != NULL);
  }
  bool push(T_t value) { 
    if (entry_count >= array_size) return false;
    array[entry_count++] = value;
    return true;
  }
  T_t  pop(void) { 
    if (entry_count == 0) return (T_t)0;
    return array[--entry_count]; 
  }
  T_t  read(void) {	// read the value on top without poping it out
    if (entry_count == 0) return (T_t)0;
    return array[entry_count-1]; 
  }
  void empty(void) { entry_count = 0; }
};

}	// end of namespace ALGH

#endif  // ALGH_STACK_HPP
