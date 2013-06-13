
//
// ALGH::PriorityQueue <T>   class template
// ALGH::PriorityQueue2<K,V> class template
//
// Jaeil Choi
// last modified in Nov, 2005
//
// common:
//  - balanced and compact tree using heap
//  - three children per node (tenary heap)
// PriorityQueue<T>:
//  - allow multiple entries in the queue with same value
//  - comparison function should be specified (NULL for basic data types)
// PriorityQueue2<K,V>:
//  - every entry has unique key,
//  - multiple 'pushing' for a key updates its value and rank in the queue.
//  - hash function should be specified (pre-defined: hash_integer, hash_pointer, hash_string)
//

#include <iostream>
#include <cmath>

#ifndef ALGH_PRIORITY_QUEUE_H
#define ALGH_PRIORITY_QUEUE_H

namespace ALGH {
  
using namespace std;
  
// ===================================================================
// PriorityQueue<K>		priority queue
// ===================================================================

template <class T>
class PriorityQueue {
private:
  int    len, max;
  T*     array;
  int    (*Cmp)(T data1, T data2);
  // ascending  : if ( < ) return -1; else if ( > ) return +1; else return 0;
  // descending : if ( < ) return +1; else if ( > ) return -1; else return 0;

public:
  PriorityQueue() : len(0), max(0), array(NULL), Cmp(NULL) { set(NULL, 64); }
  PriorityQueue(int (*cmp)(T, T), int init_size=64) : len(0), max(0), array(NULL) { set(cmp, init_size); }
  ~PriorityQueue()  { clear(false); }
  void  clear(bool free_data=false) { 
    if (array) free (array);   array = NULL;
    max = len = 0;
  }
  void  clearAndFree(void) { 
    if (array) {
      for (int i = 0; i < len; i++) if ((void*)array[i]) free((void*)array[i]);
      free (array); 
    }
    array = NULL;
    max = len = 0;
  }
  
  void  set(int (*cmp)(T, T), int init_size=64) {
    if (array) free(array);
    array = (T*) calloc (init_size, sizeof(T));
    max = init_size;  len = 0;  Cmp = cmp;
  }
  
  void  push(T data) {
    if (len >= max - 1) setStorage(0);	// double the size of storage
    BubbleUp(len++, data);
  }

  T pop(void) {
    if (len == 0) return NULL;
    T top    = array[0];
    T bottom = array[--len];
    TrickleDown(0, bottom);
    return top;
  }

  T     Top   (void)  { return (len <= 0 ? NULL : array[0]); }
  bool  isEmpty (void)  { return (len <= 0); }
  int   size  (void)  { return len; }
  void  setStorage(int size = 0) {
    if      (size == 0 )  max *= 2;
    else if (size > len)  max = size;
    else return;
    if (max < 16) max = 16;
    array = (T*) realloc(array, max * sizeof(T));
    if (array == NULL) {
      fprintf(stderr, "Error: cannot allocate %ld bytes memory for priority queue\n", max * sizeof(T));
      exit(EXIT_FAILURE);
    }
  }
  void  print (void (*print)(T), int count = 0) {
    int i, j = 0, k = 0;
    if  (print == NULL)  return;
    if  (count == 0) count = len;
    printf("Priority_queue size(%d/%d)\n", len, max);
    for (i = 0; i < len && i < count; i++) { 
      printf(" ");  print(array[i]);  k++;
      if (k % 3 == 0)  printf("  ");
      if (k == (int)pow(3.0, j)) { printf("\n");  k = 0;  j++; }
    }
    printf("\n");
  }
  
private:
  void  BubbleUp(int i, T data) {		// for  PUSH
    int parent = (i - 1) / 3;
    if ( i > 0 && Cmp(array[parent], data) > 0 ) {
      array[i] = array[parent];
      BubbleUp( parent, data );
    } else  array[i] = data;
  }
  void  TrickleDown(int i, T data) {		// for  POP
    int c1 = i+i+i+1, c2 = c1+1, c3 = c2+1, winner;
    if (c3 < len) {
      if ( Cmp(array[c1], array[c2]) > 0 ) {
	winner = ( Cmp(array[c2], array[c3]) > 0 ) ? c3 : c2;
      } else {
	winner = ( Cmp(array[c1], array[c3]) > 0 ) ? c3 : c1;
      }
    } else if (c2 < len) {
      winner = ( Cmp(array[c1], array[c2]) > 0 ) ? c2 : c1;
    } else winner = (c1 < len) ? c1 : i;
    if (winner != i && Cmp(data, array[winner]) > 0) {
      array[i] = array[winner];
      TrickleDown(winner, data);
    } else  array[i] = data;
  }

};

}	// end of namespace ALGH

// ===================================================================
// PriorityQueue2<K,V>		priority queue with unique entries
// ===================================================================

#include "algh_hash_table.hpp"

namespace ALGH {
  
template <class K, class V>
class PQUNode {
public:
  K	key;
  V	value;
  PQUNode() {}
  PQUNode(K key, V value) { this->key = key;  this->value = value; }
  ~PQUNode() {}
};

template <class K, class V>
class PriorityQueue2 {
private:
  ALGH::HashTable<K,int> htable;
  PQUNode<K,V>     *array;	// array for the heap
  int               len, max;
  bool              ascending;
  
public:
  PriorityQueue2() : array(NULL), len(0), max(0) {}
  PriorityQueue2(bool asc, int size, int (*hash_ftn)(K key)) 
    : array(NULL), len(0), max(0) { set(asc, size, hash_ftn); }
  ~PriorityQueue2() { clear(false); }
  void  clear(bool free_data=false) { 
    if (array) free (array);   array = NULL;
    max = len = 0;
    htable.clear();
  }
  void  clearAndFree(void) { 
    if (array) {
      for (int i = 0; i < len; i++) if ((void*)array[i].key) free((void*)array[i].key);
      free (array); 
    }
    array = NULL;
    max = len = 0;
    htable.clear();
  }
  
  void set(bool ascending, int size, int (*hash_ftn)(K key)) {
    // ascending: 
    // size     : the size of heap array and hash table.
    //           (set it to the maximum number of entries in the queue)
    // hash_ftn : hash function (pre-defined: hash_integer, hash_pointer, hash_string)
    if (array) free(array);
    array = (PQUNode<K,V>*) calloc (size, sizeof(PQUNode<K,V>));
    max = size;  len = 0;
    this->ascending = ascending;
    htable.create( size, hash_ftn );
  }

  bool push(K key, V value) {
    // push a (key, value) tuple, updating the value.
    // Same key can be 'pushed' multiple times, without creating multiple entries.
    int pos;
    PQUNode<K,V> data(key, value);
    if (htable.lookup( key, &pos )) {
      int c = Compare(data, array[pos]);
      //cout << "key " << key << " found at " << pos << "  c=" << c << endl;
      if      (c < 0) BubbleUp(pos, data);
      else if (c > 0) TrickleDown(pos, data);
      return false;	// existing entry
    } else {
      if (len >= max - 1) setStorage(0);	// double the size of storage
      BubbleUp(len++, data);
      return true;	// new entry
    }
  }

  bool pop(K *key, V *value) {
    // pop a (key, value) tuple from the queue
    if (len == 0) return false;
    PQUNode<K,V> top    = array[0];
    PQUNode<K,V> bottom = array[--len];
    TrickleDown(0, bottom);
    *key   = top.key;
    *value = top.value;
    return true;
  }
  
  bool Top(K *key, V *value) { 
    // Query the (key, value) tuple at the top of the queue
    if (len <= 0) return false;
    *key   = array[0].key;
    *value = array[0].value;
    return true;
  }
  
  bool find(K key, V *value) { 
    int pos;
    if (htable.lookup( key, &pos )) {
      *value = array[pos].value;
      return true;
    } else return false;
  }
  
  bool  isEmpty (void)  { return (len <= 0); }
  int   Size  (void)  { return len; }
  void  setStorage(int size = 0) {
    if      (size == 0 )  max *= 2;
    else if (size > len)  max = size;
    else return;
    if (max < 16) max = 16;
    array = (PQUNode<K,V>*) realloc(array, max * sizeof(PQUNode<K,V>));
    if (array == NULL) {
      fprintf(stderr, "Error: cannot allocate %ld bytes memory for priority queue\n", max * sizeof(PQUNode<K,V>));
      exit(EXIT_FAILURE);
    }
  }
  void  printInfo (int count = 0) {
    int i, j = 0, k = 0;
    if  (count == 0) count = len;
    cout << "Priority_queue size (" << len << "/" << max << ")" << endl;
    for (i = 0; i < len && i < count; i++) { 
      cout << " (" << array[i].key << ", " << array[i].value << ")";
      if (++k % 3 == 0)  printf("  ");
      if (k == (int)pow(3.0, j)) { printf("\n");  k = 0;  j++; }
    }
    if (!isEmpty()) printf("\n");
    htable.printInfo(false);
  }
  
private:
  int  Compare(PQUNode<K,V> d1, PQUNode<K,V> d2) {
    int n = ((d1.value < d2.value) ? -1 : ((d1.value > d2.value) ? +1 : 0));
    return (ascending ? n : -n);
  }
  void BubbleUp(int i, PQUNode<K,V> data) {	// for  PUSH
    int parent = (i - 1) / 3;
    if ( i > 0 && Compare(array[parent], data) > 0 ) {
      array[i] = array[parent];
      BubbleUp( parent, data );
    } else {
      array[i] = data;
    }
    htable.insert( array[i].key, i, false );	// update the position
  }
  void TrickleDown(int i, PQUNode<K,V> data) {	// for  POP
    int c1 = i+i+i+1, c2 = c1+1, c3 = c2+1, winner;
    if (c3 < len) {
      if ( Compare(array[c1], array[c2]) > 0 ) {
	winner = ( Compare(array[c2], array[c3]) > 0 ) ? c3 : c2;
      } else {
	winner = ( Compare(array[c1], array[c3]) > 0 ) ? c3 : c1;
      }
    } else if (c2 < len) {
      winner = ( Compare(array[c1], array[c2]) > 0 ) ? c2 : c1;
    } else winner = (c1 < len) ? c1 : i;
    if (winner != i && Compare(data, array[winner]) > 0) {
      array[i] = array[winner];
      TrickleDown(winner, data);
    } else {
      array[i] = data;
    }
    htable.insert( array[i].key, i, false );	// update the position
  }
  
};


}	// end of namespace ALGH

#endif  // ALGH_PRIORITY_QUEUE_H


// ===================================================================
// #if 0    // Example Code for PriorityQueue<T>

// #include <iostream>
// #include "algh_priority_queue.hpp"
// using namespace std;

// class CData {
// public:
//   int   i;
//   float dist;
// public:
//   CData(int ii, float dd) { i = ii; dist = dd; }
//   ~CData() {}
// };

// static int compare_ascending (CData *d1, CData* d2)
// {
//   if      (d1->dist < d2->dist) return -1;
//   else if (d1->dist > d2->dist) return +1;
//   else return 0;
// }

// static int compare_descending (CData *d1, CData* d2)
// {
//   if      (d1->dist < d2->dist) return +1;
//   else if (d1->dist > d2->dist) return -1;
//   else return 0;
// }

// void print(CData* p)
// {
//   printf("%.4f", p->dist);
// }

// int main(void)
// {
//   ALGH::PriorityQueue<CData*> q (compare_descending, 400);
//   CData *px;
  
//   q.push( new CData(1, 0.4) );
//   q.push( new CData(0, 0.3) );
//   q.push( new CData(0, 0.1) );
//   q.push( new CData(0, 0.2) );
//   q.push( new CData(0, 0.7) );
//   q.push( new CData(0, 0.5) );
//   q.push( new CData(1, 0.42) );
//   q.push( new CData(0, 0.32) );
//   q.push( new CData(0, 0.12) );
//   q.push( new CData(0, 0.22) );
//   q.push( new CData(0, 0.72) );
//   q.push( new CData(0, 0.52) );
//   px = (CData*) q.pop();  printf("(%2d, %7.4f)\n", px->i, px->dist);  delete(px);
//   px = (CData*) q.pop();  printf("(%2d, %7.4f)\n", px->i, px->dist);  delete(px);
//   px = (CData*) q.pop();  printf("(%2d, %7.4f)\n", px->i, px->dist);  delete(px);
//   q.push( new CData(0, 0.11) );
//   q.push( new CData(0, 0.21) );
//   q.push( new CData(0, 0.71) );
//   q.print(print);
//   while (!q.isEmpty()) {
//     px = q.pop();  printf("(%2d, %7.4f)\n", px->i, px->dist);  delete(px);
//   }
//   return 0;
// }

// #endif
// ===================================================================


// ===================================================================
// #if 0    // Example Code for PriorityQueue2<K,V>

// #include <iostream>
// #include "algh_priority_queue_unique.hpp"
// using namespace std;

// int main(void)
// {
//   ALGH::PriorityQueue2<int, float> q (true, 10, hash_integer);
//   int key;
//   float value;
  
//   q.push( 1, 0.1 );
//   q.push( 2, 0.2 );
//   q.push( 3, 0.3 );
//   q.push( 4, 0.4 );
//   q.push( 5, 0.5 );
//   q.push( 9, 0.9 );
//   q.push( 8, 0.8 );
//   q.push( 7, 0.7 );
//   q.push( 6, 0.6 );
//   q.printInfo();
//   q.push( 1, 0.91 );
//   q.printInfo();
//   q.push( 2, 0.92 );
//   q.push( 3, 0.93 );
//   q.printInfo();
//   q.pop( &key, &value );  printf("(%2d, %7.4f)\n", key, value);
//   q.pop( &key, &value );  printf("(%2d, %7.4f)\n", key, value);
//   q.pop( &key, &value );  printf("(%2d, %7.4f)\n", key, value);
//   q.push( 0, 0.00 );
//   q.printInfo();
//   while (!q.isEmpty()) {
//     q.pop( &key, &value );  printf("(%2d, %7.4f)\n", key, value);
//   }
//   q.printInfo();
  
//   return 0;
// }

// #endif
// ===================================================================
