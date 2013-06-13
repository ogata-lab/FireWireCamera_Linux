
//
// ALGH::Set<T> class template
//
// Jaeil Choi
// last modified in Mar, 2007
//
// This class is designed for:
//   - maintaining a set of unique entries,
//     while adding up / subtracting the weight of each entry.
//   - providing basic set operations -- add(v) find(v) remove(v)
//   - ordered retreaval -- pop()
//   - finding the most frequent entry -- findMostFrequent()
//   - serialization -- convertToArray()
// Implementation detail:
//   - This class uses a binary tree internally.
//   - If you want multiple instances for an entry, then use ALGH::Group.
//   - If you need (key,data) pairs, then use ALGH::BinaryTree.
// Usage example:
//   Set<int> iset;
//   iset.add(5);  iset.add(3);  iset.add(8);  iset.add(5);  iset.add(7);
//   cout << "set : " << iset << "   size = << iset.count << endl;
//   iset.remove(5);   bool yn = iset.find(5);
//   iset.convertToArray( array );
//


#ifndef ALGH_SET_HPP
#define ALGH_SET_HPP

#include <iostream>

namespace ALGH {
  
using namespace std;


// ===================================================================
// SetNode  class template
// ===================================================================

template <class T>
class SetNode {
public:
  T		value;
  float		weight;
  SetNode<T>	*left, *right;
public:
  SetNode<T>(T v, float w=1.0f) : left(NULL), right(NULL) { value = v; weight = w; }
  ~SetNode<T>() { clear(); }
  void clear() {
    if (left)  left->clear();  if (right) right->clear();
    free( this ); // SetNode should be a dynamic variable.
  }
  void clearAndFree( ) {
    if (left)  left->clearAndFree();  if (right) right->clearAndFree();
    if ((void*)(this->value)) free((void*)(this->value));
    free( this ); // SetNode should be a dynamic variable.
  }

  SetNode<T>* add(T v, int (*compare)(T v1, T v2), bool &added, float w=1.0f) {
    int sign;
    if (compare) sign = compare(v, this->value);
    else sign = (v < this->value) ? -1 : ((v > this->value) ? +1 : 0);
    if        (sign < 0) {
      if (left) left->add( v, compare, added, w );
      else { left = new SetNode(v, w); added = true; }
    } else if (sign > 0) {
      if (right) right->add( v, compare, added, w );
      else { right = new SetNode(v, w); added = true; }
    } else { 
      this->weight += w;  added = false; 
    }
    return this;
  }
      
  SetNode<T>* remove(T value, int (*compare)(T v1, T v2), bool &removed, float *w) {
    int sign;
    if (compare) sign = compare(value, this->value);
    else sign = (value < this->value) ? -1 : ((value > this->value) ? +1 : 0);
    if      (sign < 0 && left)  left  = left->remove( value, compare, removed, w );
    else if (sign > 0 && right) right = right->remove( value, compare, removed, w );
    else {
      removed = true;
      if (w) *w = this->weight;
      if (right) {
	right = right->pop( &this->value, true, w );
      } else if (left) {
	SetNode<T> *temp = left;
	free( this );  return temp;
      } else {
	free( this );  return NULL;
      }
    }
    return this;
  }
  
  SetNode<T>*  pop(T *vp, bool smallest, float *w) {
    SetNode<T> *cp1, *cp2;
    if (smallest) { cp1 = left;  cp2 = right; }
    else          { cp1 = right; cp2 = left;  }
    if (cp1) {
      if (smallest) left  = cp1->pop(vp, smallest, w);
      else          right = cp1->pop(vp, smallest, w);
      return this;
    } else {
      *vp = value;
      if (w) *w = this->weight;
      free( this );
      return (cp2 ? cp2 : NULL);
    }
  }
  
  bool find(T value, int (*compare)(T v1, T v2), float *weight=NULL) {
    int sign;
    if (compare) sign = compare(value, this->value);
    else sign = (value < this->value) ? -1 : ((value > this->value) ? +1 : 0);
    if (sign == 0) { if (weight) *weight = this->weight; return true; }
    else if (sign <  0) { if (left)  return  left->find(value, compare, weight); }
    else if (sign >  0) { if (right) return right->find(value, compare, weight); }
    return false;
  }
  
  T findMostFrequent(float *w=NULL) {
    int lv, rv;  float lw, rw;
    if (left)  lv =  left->findMostFrequent( &lw ); else lw = lv = 0; 
    if (right) rv = right->findMostFrequent( &rw ); else rw = rv = 0; 
    if (this->weight >= lw) {
      if (this->weight >= rw) {
	if (w) *w = this->weight;  return this->value;
      } else        { if (w) *w = rw; return rv; }
    } else {
      if (lw >= rw) { if (w) *w = lw; return lv; }
      else          { if (w) *w = rw; return rv; }
    }
  }
  
  void convertToArray(T *array, int *pos) {
    if (left)  left->convertToArray(array, pos);
    array[(*pos)++] = value;
    if (right) right->convertToArray(array, pos);
  }
  
  void print(void (*print)(T v)) {
    if (left) left->print(print);
    if (print) { print(value);  cout << " "; }
    else cout << value << " ";
    if (right) right->print(print);
  }
  
};

// ===================================================================
// Set  class template
// ===================================================================

template <class T>
class Set {
public:
  int		count;
  SetNode<T>	*root;
  double	wsum;
  int  (*compare)(T v1, T v2);
  void (*print)(T v);
  
public:
  Set() : count(0), root(NULL), wsum(0), compare(NULL), print(NULL) {}
  Set(int (*compare)(T v1, T v2), void (*print)(T v)) 
    : count(0), root(NULL), wsum(0) { 
    this->compare = compare;  this->print = print;
  }
  ~Set() { clear(); }
  void clear() { if (root) root->clear(); root = NULL; wsum = count = 0; }
  void clearAndFree() { if (root) root->clearAndFree(); root = NULL; wsum = count = 0; }
  
  // -----------------------------------------------------------------
public:
  
  inline bool isEmpty(void) { return (root ? false : true); }
  inline int  size(void) { return count; }
  inline void setCompareFunction(int (*compare)(T v1, T v2)) { this->compare = compare; }
  inline void setPrintFunction(void (*print)(T v1)) { this->print = print; }
  
  bool add(T value, float weight=1.0f) {
    bool added = false;
    if (root) {
//       if (sizeof(T) == sizeof(void*) && compare && value == root->value)
// 	cerr << "Warning: (Set) Do not add a temporary pointer to Set." << endl;
      root = root->add( value, compare, added, weight );
      if (added) count++;
      wsum += weight;
    } else {
      root = new SetNode<T>(value, weight);
      count = 1;
      wsum  = weight;
      added = true;
    }
    return added;
  }
  
  bool remove(T value) { 
    // Remove 'value' from the set
    bool  removed = false;
    if (root) {
      float weight;
      root = root->remove( value, compare, removed, &weight );
      if (removed) { count--; wsum -= weight; }
    }
    return removed;
  }
  
  T    pop (bool smallest=true, float *weight=NULL) {
    T value;
    if (root) {
      float wdummy;
      if (weight == NULL) weight = &wdummy;
      root = root->pop( &value, smallest, weight );
      count--;  wsum -= *weight;
    } else value = (T)0;
    return value;
  }
  
  bool find (T value, float *weight=NULL) {
    // find the value in the set, and return the number of addition of it
    if (root) return root->find(value, compare, weight);
    else return false;
  }
  
  T findMostFrequent(float *weight=NULL) {
    return root->findMostFrequent( weight );
  }
  
  void convertToArray(T *array) {
    if (root == NULL) return;
    int pos = 0;
    root->convertToArray( array, &pos ); 
  }
  
};

template <class T>
static ostream &operator<< (ostream &os, const Set<T> &s) {
  cout << "{ ";
  if (s.root) s.root->print(s.print);
  cout << "}";
  return os;
}


}	// end of namespace ALGH

#endif  // ALGH_SET_HPP



// ===================================================================
// example code
// ===================================================================

#if 0

#include <iostream>
#include "algh_set.hpp"
using namespace std;

// int  compare(CNode *a, CNode *b) { return (a->key < b->key) ? -1 : ((a->key > b->key) ? +1 : 0); }
// void print(CNode *c) { cout << c->value; }
void print(int c) { cout << "[" << c << "]"; }

int main(void)
{
  ALGH::Set<int> iset(NULL, print);
  int value;
  
  iset.add(5);  iset.add(3);  iset.add(8);  iset.add(5);  iset.add(7);
  cout << "Set count : " << iset.count << "    " << iset << endl;
  
  iset.add(7);  iset.add(1);  iset.add(3);  iset.add(9);
  iset.add(4);  iset.add(3);  iset.add(2);  iset.add(8);
  
  cout << "Set count : " << iset.count << "    " << iset << endl;
  cout << "  find() ?  " 
       << "  1:" << iset.find(1) << "  9:" << iset.find(9)
       << "  0:" << iset.find(0) << "  6:" << iset.find(6) << endl;
  
  int *array = (int*)malloc( iset.count * sizeof(int) );
  iset.convertToArray( array );
  cout << "convertToArray : ";
  for (int i = 0; i < iset.count; i++) cout << array[i] << ", ";
  cout << endl;
  
  cout << "  pop()   : ";
  while(iset.count > 0) {
    value = iset.pop(false);
    cout << value << " " << flush;
  }
  cout << endl;
  
  iset.add(5);  iset.add(3);  iset.add(8);  iset.add(5);  iset.add(7);
  iset.clear();
}

#endif
