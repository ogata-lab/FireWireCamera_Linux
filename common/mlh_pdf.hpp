
//
// MLH::PDF<T> 
// probability distribution function class for machine learning
//
// Jaeil Choi
// Last modified in July, 2009
//


#ifndef MLH_PDF_HPP
#define MLH_PDF_HPP

#include <iostream>
#include <cfloat>

namespace MLH {
  
template <class T>
class PDF {
public:
  int size;
  T   *data;
private:
  bool allocated;
  
public:
  PDF() : size(0), data(NULL), allocated(false) {}
  PDF(int size, char init_mode=' ') 
    : size(0), data(NULL), allocated(false) { setPDF(size, init_mode); }
  PDF(int size, T *data, char init_mode=' ') 
    : size(0), data(NULL), allocated(false) { setPDFWithBuffer(size, data, init_mode); }
  ~PDF(void) { clear(); }
  void clear(void) { 
    if (data && allocated) free(data); 
    size = 0;  data = NULL;  allocated = false; 
  }
  
  void setPDF(int size, char init_mode=' ') {
    clear();
    this->data = (T*)malloc(size*sizeof(T));
    this->size = size;
    allocated = true;
    switch (init_mode) {
    case '0': memset( this->data, 0, size*sizeof(T) );   break;  // zero
    case '1': clearValues(1.0f);                         break;  // 1.0
    case 'U': clearValues((size>0 ? 1.0f/size : 1.0f )); break;  // uniform
    case 'N': clearValues((size>0 ? 1.0f/size : 1.0f )); break;  // normalize
    default : break;  // leave it uninitialized
    }
  }
  void setPDFWithBuffer(int size, T *data, char init_mode=' ') {
    clear();
    this->data = data;
    this->size = size;
    allocated = false;
    switch (init_mode) {
    case '0': memset( this->data, 0, size*sizeof(T) );   break;  // zero
    case '1': clearValues(1.0f);                         break;  // 1.0
    case 'U': clearValues((size>0 ? 1.0f/size : 1.0f )); break;  // uniform
    case 'N': normalize();                               break;  // normalize
    default : break;  // leave it uninitialized
    }
  }
  T& operator[](int idx) { return data[idx]; }
  void clearValues(T v=0) { for (int i=0; i<size; i++) data[i] = v; }
  
  // -----------------------------------------------------------------
  
  T getSum(void) {
    T sum = 0;
    for (int i = 0; i < size; i++) sum += data[i];
    return sum;
  }
  bool normalize(T npdf[]=NULL) {
    T sum = getSum();
    if (sum == 0) return false; 
    if (npdf) for (int i=0; i<size; i++) npdf[i] = data[i] / sum;
    else      for (int i=0; i<size; i++) data[i] /= sum;
    return true;
  }
  
  int getRandomSample(void) {
    // Take a random sample among the entries, assumming normalized PDF.
    double sum=0, pr = rand() / (double)RAND_MAX;
    for (int i=0; i<size; i++) {
      sum += data[i];
      if (pr <= sum) return i;
    }
    return (size-1);
  }
  
  // -----------------------------------------------------------------
  
  void printInfo(char *cmmt=NULL) {
    T   max=0;  int i;  char c;
    for (i=0; i<size; i++) if (data[i]>max) max = data[i];
    if (max <= 0) {
      printf("%s(%d): Max=%.4f < All Zero >\n", (cmmt ? cmmt:"PDF"), size, max);
    } else if (size>5000) {
      printf("%s(%d): Max=%.4f < TooManyToShow >\n", (cmmt ? cmmt:"PDF"), size, max);
    } else {
      printf("%s(%d): Max=%.4f <", (cmmt ? cmmt:"PDF"), size, max);
      for (i=0; i<size; i++) {
	T   value = data[i] / max;
	if      (value == 0.00) c = ' ';
	else if (value == 1.00) c = 'M';
	else if (value  < 0.00) c = '-';  // invalid value
	else if (value >= 0.95) c = '*';
	else if (value <= 0.05) c = '.';
	else                   c = '0' + (int)(value*10+0.5);
	printf("%c", c);
	if (i%10==9) printf(" ");
      }
      printf(">\n");
    }
  }
  
};  
	    
	    
}	// end of namespace MLH

#endif  // MLH_PDF_HPP


