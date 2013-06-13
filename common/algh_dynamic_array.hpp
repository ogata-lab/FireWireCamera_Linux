
//
// ALGH::DynamicArray<T> class template
//
// Jaeil Choi
// last modified in Mar, 2007
//
// This class is designed for:
//   - simple array with automatic and dynamic memory allocation
//   - fast sorting using QuickSort
//   - fast searching using BinarySearch
// Implementation details:
//   - it's similar to 'vector<>' in STL (only in basic array access)
//   - the array is expanded automatically by about 10%
//   - You can access its data directly using array notation 'a[]'.
//   - The allocated memory can be shinked down by calling 'setSize(n)'.
//


#ifndef ALGH_DYNAMIC_ARRAY_HPP
#define ALGH_DYNAMIC_ARRAY_HPP

#include <iostream>

namespace ALGH {
  
using namespace std;

template <class T>
class DynamicArray {
public:
  int	size;		// current size of the array
  T	*data;		// pointer to the allocated memory
private:
  int	reserved_size;	// size of the reserved memory
public:
  DynamicArray() : size(0), data(NULL), reserved_size(0) { reserve(20); }
  DynamicArray(int reserved_size) : size(0), data(NULL), reserved_size(0) { reserve(reserved_size); }
  ~DynamicArray() { clear(); }
  void clear(bool free_data=true) { if (data && free_data) free(data); data = NULL; size = reserved_size = 0; }
  void clearAndFree(void) { 
    if (data) for (int i=0; i<size; i++) if (data[i]) free(data[i]);
    clear();
  }
  
  void reserve(int new_size) {		// reserve memory space in advance
    if (new_size <= reserved_size) return;
    data = (T*) realloc (data, new_size * sizeof(T));
    reserved_size = new_size;
  }
  
  void setSize(int new_size) {
    if (new_size > reserved_size) {	// expand (10% more)
      reserved_size = new_size + (new_size>29 ? new_size/10 : 3);
      data = (T*) realloc (data, reserved_size * sizeof(T));
      memset( (data + size), 0, (reserved_size - size) * sizeof(T) );
    } else if (new_size < size) {	// shrink, only if specified
      if (new_size < 0) new_size = 0;
      reserved_size = new_size;
      data = (T*) realloc (data, reserved_size * sizeof(T));
    }
    size = new_size;
  }
  
  T &operator[](int index) {
    if (index >= size) setSize( index + 1 );
    return data[index];
  }
  
  int add(void) {
    setSize( size + 1 );  // increase memory, if necessary
    return size - 1;
  }
  
  int add(T value) { 
    setSize( size + 1 );  // increase memory, if necessary
    data[ size - 1 ] = value;
    return size - 1;
  }
  
  int find(T value) {
    for (int i = 0; i < size; i++) if (data[i] == value) return i;
    return -1;
  }
  
  bool removeByIndex(int idx, bool preserve_order=false) {
    if (idx < 0 || idx >= size) return false;
    if (preserve_order) for (; idx < size-1; idx++) data[idx] = data[idx+1];
    else                data[idx] = data[size-1];
    size--;
    return false;
  }
  
  bool removeByValue(T value) {
    int idx = find( value );
    if (idx < 0) return false;
    for (; idx < size-1; idx++) data[idx] = data[idx+1];
    size--;
    return true;
  }
  
  void copyFrom(DynamicArray<T> *da, bool preserve_old_values) {
    if (!da) return;
    if (preserve_old_values) reserve( da->size + this->size );
    else { this->size = 0;   reserve( da->size ); }
    for (int i = 0; i < da->size; i++) add( da->data[i] );
  }
  void swapWith(DynamicArray<T> *da) {
    DynamicArray<T> tmp = *this;
    *this = *da;
    *da   = tmp;
    tmp.data = NULL; 
    tmp.reserved_size = tmp.size = 0;
  }
  
  // -----------------------------------------------------------------
  // sorting / binary search
  // -----------------------------------------------------------------
private:
  int  (*cmp)(T entry1, T entry2);	// compare function
  
public:
  void sort(int(*compare)(T entry1, T entry2)=NULL) {
    // Quick sort using given compare function.
    //   compare() function will return 0, +1, or -1, given two entries.
    //   if compare==NULL, then generic '-' operator will be used to compare.
    this->cmp = compare;
    runQuickSort(0, size-1);
  }
  int search(T entry, int(*compare)(T entry1, T entry2)=NULL) {
    // Binary search using given compare function, and return the index. (-1 if failed)
    //   compare() function will return 0, +1, or -1, given two entries.
    //   if compare==NULL, then generic '-' operator will be used to compare.
    this->cmp = compare;
    return BinarySearch( entry, 0, size-1 );
  }
private:
  void runQuickSort(int sidx, int eidx) {
    // Sort the array 'data[]' in the index range [sidx,eidx]
    if (sidx >= eidx) return;
    int midx = runQuickSortPartition( sidx, eidx );
    runQuickSort( sidx, midx );
    runQuickSort( midx+1, eidx );
  }
  int  runQuickSortPartition(int sidx, int eidx) {
    T    x = data[sidx];
    int  i = sidx-1, j = eidx+1;
    if (cmp) {
      while (1) {
	do { j = j - 1; } while ( cmp(data[j], x) > 0 );
	do { i = i + 1; } while ( cmp(data[i], x) < 0 );
	if (i < j) { T tmp = data[i]; data[i] = data[j]; data[j] = tmp; }
	else return j;
      }
    } else {
      while (1) {
	do { j = j - 1; } while ( data[j] - x > 0 );
	do { i = i + 1; } while ( data[i] - x < 0 );
	if (i < j) { T tmp = data[i]; data[i] = data[j]; data[j] = tmp; }
	else return j;
      }
    }
  }
  int BinarySearch(T entry, int sidx, int eidx) {
    if (sidx > eidx) return -1;
    int ret, midx = (sidx + eidx) / 2;
    if (cmp) ret = cmp( entry, data[midx] );
    else     ret = entry - data[midx];
    if      (ret < 0) return BinarySearch( entry, sidx, midx-1 );
    else if (ret > 0) return BinarySearch( entry, midx+1, eidx );
    else return midx;
  }
  
};

// template <class T>
// static ostream &operator<< (ostream &os, const DynamicArray<T> &da) {
//   os << "[ ";
//   for (int i = 0; i < da.size; i++) os << da.data[i] << " ";
//   os << "]";
//   return os;
// }

}	// end of namespace ALGH

#endif // ALGH_DYNAMIC_ARRAY_HPP
