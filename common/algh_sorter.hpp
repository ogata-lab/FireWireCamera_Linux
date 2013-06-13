
//
// ALGH::Sorter <T>   class template 
// Sorter2<T,K> class template for (entry,value) tuples
//
// Jaeil Choi
// last modified in Nov, 2005
//
//

#ifndef ALGH_SORTER_HPP
#define ALGH_SORTER_HPP

#include <iostream>

namespace ALGH {
  
using namespace std;

#define SWAP(a,b) do { temp = a; a = b; b = temp; } while(0)
#define SWAP3(a,b,t) do { t = a; a = b; b = t; } while(0)

// ===================================================================
// Sorter<T>			Sorting singleton entries
// ===================================================================

template <class T>
class Sorter {
  
public:
  Sorter()  {}
  ~Sorter() {}
private:
  int (*compare_ftn)(T a, T b);	// return the sign of (a - b) for ascending order
  
  // -----------------------------------------------------------------
  // QuickSort
  // -----------------------------------------------------------------
  
public:
  void QuickSort(T array[], int size, bool asc) {
    if (!array) return;
    QuickSortRecursive( array, 0, size-1 );
    if (!asc) {
      T temp;
      for (int i = 0; i < size/2; i++) SWAP( array[i], array[size-i-1] );
    }
  }
  void QuickSort(T array[], int size, int (*compare)(T, T)) {
    if (!array) return;
    this->compare_ftn = compare;
    QuickSortRecursive2( array, 0, size-1 );
  }
  
private:  
  void QuickSortRecursive(T array[], int p, int r) {
    if (p >= r) return;
    int q = QuickSortPartition(array,p,r);
    QuickSortRecursive( array, p, q );
    QuickSortRecursive( array, q+1, r );
  }
  void QuickSortRecursive2(T array[], int p, int r) {
    if (p >= r) return;
    int q = QuickSortPartition2(array,p,r);
    QuickSortRecursive2( array, p, q );
    QuickSortRecursive2( array, q+1, r );
  }
  int  QuickSortPartition(T array[], int p, int r) {
    T    x = array[p], temp;
    int  i = p-1, j = r+1;
    while (1) {
      do { j = j - 1; } while (array[j] > x);
      do { i = i + 1; } while (array[i] < x);
      if (i < j)  SWAP( array[i], array[j] );
      else return j;
    }
  }
  int  QuickSortPartition2(T array[], int p, int r) {
    T    x = array[p], temp;
    int  i = p-1, j = r+1;
    while (1) {
      do { j = j - 1; } while ( compare_ftn( array[j], x ) > 0 );
      do { i = i + 1; } while ( compare_ftn( array[i], x ) < 0 );
      if (i < j)  SWAP( array[i], array[j] );
      else return j;
    }
  }
  
  // -----------------------------------------------------------------
  // Shell Sort
  // -----------------------------------------------------------------
public:
  void ShellSort(T array[], int size, bool asc) {
    int i, n;
    if (asc) {
      for (n=size/2; n>=1; n /= 2)
	for (i=0; i<n; i++) InsertionSortAsc(array+i, size-i, n);
    } else {
      for (n=size/2; n>=1; n /= 2)
	for (i=0; i<n; i++) InsertionSortDesc(array+i, size-i, n);
    }
  }
  void InsertionSortAsc(T array[], int size, int step) {
    int  i, j, k;
    for (i=step; i<size; i+=step) {
      T val = array[i];
      for (j=0; j<i; j+=step) if (array[j] > val) break; // find insertion point
      if  (j==i) continue;
      for (k=i; k>j; k-=step) array[k] = array[k-step];  // insert 'val' at [j]
      array[j] = val;
    }
  }
  void InsertionSortDesc(T array[], int size, int step) {
    int  i, j, k;
    for (i=step; i<size; i+=step) {
      T val = array[i];
      for (j=0; j<i; j+=step) if (array[j] < val) break; // find insertion point
      if  (j==i) continue;
      for (k=i; k>j; k-=step) array[k] = array[k-step];  // insert 'val' at [j]
      array[j] = val;
    }
  }
  
  // -----------------------------------------------------------------
  // 
  // -----------------------------------------------------------------
public:
  void findSmallerValues(T array[], int total, int n) {
    // Find 'n' smaller values from all the values in 'array[]',
    //   by placing them in the first 'n' slots of the array.
    int  i, j;  T cv, tv;
    for (i=n; i<total ; i++) 
      for (j=0,cv=array[i]; j<n; j++) 
	if (cv<array[j]) { tv=array[j]; array[j]=cv; cv=tv; }
  }
  
};  
	    
	    
// ===================================================================
// Sorter2<T,K>		Sorting entry-value tuples
// ===================================================================

template <class T, class K>
class Sorter2 {
public:
  Sorter2()  {}
  ~Sorter2() {}
  // -----------------------------------------------------------------
  // QuickSort
  // -----------------------------------------------------------------
public:
  void QuickSort(T array[], K value[], int size, bool asc) {
    // sort entries in 'array[]' based on their values in 'value[]'
    QuickSortRecursive( array, value, 0, size-1 );
    if (asc) return;
    T tmpT;  K tmpK;
    for (int i = 0; i < size/2; i++) {
      SWAP3( array[i], array[size-i-1], tmpT );
      SWAP3( value[i], value[size-i-1], tmpK );
    }
  }
private:  
  void QuickSortRecursive(T array[], K value[], int p, int r) {
    if (p >= r) return;
    int q = QuickSortPartition(array, value, p, r);
    QuickSortRecursive( array, value, p, q );
    QuickSortRecursive( array, value, q+1, r );
  }
  int  QuickSortPartition(T array[], K value[], int p, int r) {
    K    x = value[p], tmpK;  T tmpT;
    int  i = p-1, j = r+1;
    while (1) {
      do { j = j - 1; } while (value[j] > x);
      do { i = i + 1; } while (value[i] < x);
      if (i < j) {
	SWAP3( array[i], array[j], tmpT );
	SWAP3( value[i], value[j], tmpK );
      }
      else return j;
    }
  }
};  


}	// end of namespace ALGH

#endif // ALGH_SORTER_HPP

	    
// #if 0  // Example ====================================================

// #include <iostream>
// #include "algh_sorter.hpp"
// using namespace std;

// class CNode {
// public:
//   int value;
// public:
//   CNode(int value) { this->value = value; }
//   ~CNode(void) {}
// };

// static int compare(CNode* n1, CNode* n2) 
// {
//   // This function must return ( n1 - n2 ) for ascending order
//   if      (n1->value < n2->value) return -1;
//   else if (n1->value > n2->value) return +1;   // THIS is important!
//   else return 0;
// }

// int main(void)
// {
//   ALGH::Sorter<int> sorter1;
//   int a1[] = { 0, 4, 3, 2, 1, 8, 6, 7, 5, 9, 3, 10, 9 };
//   int size = 13;
//   for (int i=0; i<size; i++) printf("%d ", a1[i]);  printf("\n");
//   sorter1.QuickSort(a1, size, true);
//   for (int i=0; i<size; i++) printf("%d ", a1[i]);  printf("quicksort\n");
    
//   ALGH::Sorter<CNode*> sorter2;
//   CNode* a2[10];
//   for (int i = 0; i<size; i++) a2[i] = new CNode(a1[i]);
//   sorter2.QuickSort(a2, size, compare);
//   for (int i=0; i<size; i++) printf("%d ", a2[i]->value);  printf("quicksort of CNode\n");
    
//   return EXIT_SUCCESS;
// }

// #endif  // ===========================================================

