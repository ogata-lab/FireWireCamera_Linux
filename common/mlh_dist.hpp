
//
// MLH::Dist<T> 
// distribution function class for machine learning
//
// Jaeil Choi
// Last modified in Apr, 2005
//


#ifndef MLH_DIST_HPP
#define MLH_DIST_HPP

#include <iostream>
#include <cfloat>
#include "mth_probability.hpp"

namespace MLH {
  
template <class T>
class Dist {
public:
  int size[3], total;
  T   *data;
public:
  Dist()  { clear(false); }
  Dist(int size, T value = 0)  { clear(false); create(size); if (value != 0) initializeValues(value); }
  Dist(int size, bool uniform) { clear(false); create(size); initializeProb(uniform); }
  Dist(int row, int col, T value = 0)  { clear(false); create(row, col); if (value != 0) initializeValues(value); }
  Dist(int row, int col, bool uniform) { clear(false); create(row, col); initializeProb(uniform); }
  Dist(int row, int col, int hei, T value = 0)  { clear(false); create(row, col, hei); if (value != 0) initializeValues(value); }
  Dist(int row, int col, int hei, bool uniform) { clear(false); create(row, col, hei); initializeProb(uniform); }
  ~Dist() { clear(true);  }
  void clear(bool bFree=true) {
    if (bFree && data) free(data);
    data = NULL;
    total = size[0] = size[1] = size[2] = 0;
  }
  void create(int row, int col=1, int hei=1) {
    clear();
    this->size[0] = row;
    this->size[1] = col;
    this->size[2] = hei;
    this->total = row * col * hei;
    data = (T*) calloc( total, sizeof(T) );
  }
  void create(Dist<T> *d) { create(d->size[0], d->size[1], d->size[2]); }
  void initializeValues(T value) { 
    if (!data) return;
    // if 'value' < 0, then it is set to (rand() * (-value))
    for (int i = 0; i < total; i++)
      data[i] = ( value >= 0 ? value : -value * rand() / (T)RAND_MAX );
  }
  void initializeProb(bool uniform) { 
    int i;
    if (!data) return;  
    if (uniform) for (i = 0; i < total; i++) data[i] = 1.0/total;
    else {
      for (i = 0; i < total; i++) data[i] = rand() / (T)RAND_MAX;
      normalize();
    }
  }
  inline T& operator() (int idx) { return data[idx]; }
  inline T& operator() (int row, int col) { return data[ Index(row, col) ]; }
  inline T& operator() (int row, int col, int hei) { return data[ Index(row, col, hei) ]; }
  inline int Index(int row, int col) { return row * size[1] + col; }
  inline int Index(int row, int col, int hei) { return row * size[1] * size[2] + col * size[2] + hei; }
  inline int countDimensions(void) { int n;  for (n=3; n>0; n--) if (size[n-1] > 1) break;  return n; }

  
  // -----------------------------------------------------------------
  // Basic Algebra
  // -----------------------------------------------------------------
public:
  void add(T value)      { for (int i=0; i<total; i++) data[i] += value; }
  void add(Dist *d) { for (int i=0; i<total; i++) data[i] += d->data[i]; }
  void add(Dist *d1, Dist *d2) { for (int i=0; i<total; i++) data[i] = d1->data[i] + d2->data[i]; }
  void add(Dist *d1, Dist *d2, T s) { for (int i=0; i<total; i++) data[i] = d1->data[i] + d2->data[i] * s; }
  void add(Dist *d1, Dist *d2, T s1, T s2) { for (int i=0; i<total; i++) data[i] = d1->data[i] * s1 + d2->data[i] * s2; }
  void mult(T value)      { for (int i=0; i<total; i++) data[i] *= value; }
  void mult(Dist *d) { for (int i=0; i<total; i++) data[i] *= d->data[i]; }
  void mult(Dist *d1, Dist *d2) { for (int i=0; i<total; i++) data[i] = d1->data[i] * d2->data[i]; }
  
  T    Sum(void) { T sum = 0;  for (int i=0; i<total; i++) sum += data[i]; return sum; }
  T    Max(void) { T max = -FLT_MAX;  for (int i=0; i<total; i++) if (data[i] > max) max = data[i]; return max; }
  void Log(Dist *d) { create(d); for (int i=0; i<total; i++) data[i] = log(d->data[i]); }
  void Log(void) { for (int i=0; i<total; i++) data[i] = log(data[i]); }
  void Exp(void) { for (int i=0; i<total; i++) data[i] = exp(data[i]); }
  void setMin(T min) { for (int i=0; i<total; i++) if (data[i]<min) data[i]=min; }
  void setMax(T max) { for (int i=0; i<total; i++) if (data[i]>max) data[i]=max; }
  void setMinMax(T min, T max) { for (int i=0; i<total; i++) { if (data[i] < min) data[i] = min; if (data[i] > max) data[i] = max; } }
  Dist<T>& operator= (Dist<T>& dist) {
    create(dist.size[0], dist.size[1], dist.size[2]);
    memcpy (data, dist.data, total * sizeof(T));
    return *this;
  }
  
  // -----------------------------------------------------------------
  // Distribution as probability values
  // -----------------------------------------------------------------
public:
  void normalize(T lower_bound=0) {
    int i;  T sum = 0;
    for (i = 0; i < total; i++) {
      data[i] += lower_bound;
      sum += data[i];
    }
    if (sum > 0) for (i = 0; i < total; i++) data[i] /= sum;
  }
  void SumUpDimension( Dist<T> *d, int dim ) {
    int i, j, k, n = countDimensions();  
    T sum;
    if (d->size[dim] < 2) { *this = *d; return; }
    switch (dim) {
    case 0:
      create(d->size[1], d->size[2]);
      for (i = 0; i < size[0]; i++)
	for (j = 0; j < size[1]; j++) {
	  for (sum = k = 0; k < d->size[0]; k++) sum += (*d)(k,i,j);
	  data[ Index(i,j) ] = sum;
	}
      break;
    case 1:
      create(d->size[0], d->size[2]);
      for (i = 0; i < size[0]; i++)
	for (j = 0; j < size[1]; j++) {
	  for (sum = k = 0; k < d->size[1]; k++) sum += (*d)(i,k,j);
	  data[ Index(i,j) ] = sum;
	}
      break;
    case 2:
      create(d->size[0], d->size[1]);
      for (i = 0; i < size[0]; i++)
	for (j = 0; j < size[1]; j++) {
	  for (sum = k = 0; k < d->size[2]; k++) sum += (*d)(i,j,k);
	  data[ Index(i,j) ] = sum;
	}
      break;
    }
  }
  
  // -----------------------------------------------------------------
  // Distribution as log-probability values
  // -----------------------------------------------------------------
  
  void normalizeLog(T lower_bound=0) {
    // normalize the probability distribution in log-domain
    int i;
    T sum = 0;
    if (lower_bound > 0) {
      T lb = log( lower_bound );
      for (i = 0; i < total; i++) if (data[i] < lb) data[i] = lb;
    }
    add( - Max() );
    for (sum = i = 0; i < total; i++) sum += exp(data[i]);
    add( -log(sum) );
  }
  void SubtractMaxOverRow(void) {
    // for each row, find max and subtract
    int i, p, p0, p1, count = size[1] * size[2];
    for (i = 0; i < size[0]; i++) {
      p0 = i * count;  p1 = (i+1) * count;
      T maxv = -FLT_MAX;
      for (p = p0; p < p1; p++) if (data[p] > maxv) maxv = data[p];
      for (p = p0; p < p1; p++) data[p] -= maxv;
    }
  }
  void normalizeLogOverRow(T lower_bound=0) {
    // normalize over the first dimension
    if (countDimensions() != 2) { std::cerr << "Error: in normalizeLog" << std::endl; return; }
    int i, p, p0, p1, count = size[1] * size[2];
    T   sum, maxv;
    for (i = 0; i < size[0]; i++) {
      p0 = i * count;  p1 = (i+1) * count;
      sum = 0;
      for (p = p0; p < p1; p++) 	// convert into probability
	sum += data[p] = exp(data[p]) + lower_bound;
      for (p = p0; p < p1; p++) {
	data[p] /= sum;			// normalize
	data[p]  = log(data[p]);	// convert into log-probability
      }
    }
  }
  
  // -----------------------------------------------------------------
  // Machine Learning
  // -----------------------------------------------------------------
public:
  void LLGaussian( Dist<T> *x, Dist<T> *m, Dist<T> *s2, Dist<T> *ln_s2=NULL ) {
    if (size[0] != x->size[0] || size[1] != x->size[1] || size[2] != x->size[2]) 
      create( x->size[0], x->size[1], x->size[2] );
    MLH::Probability<T> prob;
    int i;
    if (ln_s2)
      for (i = 0; i < total; i++) 
	data[i] = prob.LogGaussian2( x->data[i], m->data[i], s2->data[i], ln_s2->data[i] );
    else
      for (i = 0; i < total; i++) 
	data[i] = prob.LogGaussian2( x->data[i], m->data[i], s2->data[i] );
  }
  void LL2Probability(T max = 0) {
    // convert log-likelihood values to probability values
    int i;
    T sum = 0;
    for (i = 0; i < total; i++) {
      data[i] = exp( data[i] - max );
      sum    += data[i];
    }
    if (sum > 0) for (i = 0; i < total; i++) data[i] /= sum;
  }
  void MutuallyExclusiveMax( Dist<T> *d, T lower_bound=0 ) {
    // Assuming 'this' and 'd' are probability distributions of 
    // mutually exclusive masks, set matching values of max of 'd' to zero.
    T max = d->Max();
    for (int i = 0; i < total; i++) if (d->data[i] == max) data[i] = 0;
    normalize(lower_bound);
  }
  int MaxOnly(T lower_bound=0) {
    int i, count = 0, max_idx, indices[20];
    T   max = Max();
    for (i = 0; i < total; i++) 
      if (data[i] == max && count<20) indices[count++] = i;
    i = (int)(count * (rand() / (T)RAND_MAX));
    if (i == count) i = count - 1;
    max_idx = indices[i];
    for (i = 0; i < total; i++) data[i] = (i == max_idx ? 1 : 0);
    normalize(lower_bound);
    return max_idx;
  }
  
  // -----------------------------------------------------------------
  // Etc
  // -----------------------------------------------------------------
  
  void printInfo(char *comment = NULL, int indent=0) {
    int i, j, k, ndim = countDimensions();
    char  format[40], blank[80];
    for (i = 0; i < indent; i++) blank[i] = ' ';  blank[i] = '\0';
    if (ndim == 1) {
      printf("%s%s[ ", blank, (comment ? comment : ""));
      for (i = 0; i < size[0]; i++) printf("%.2f ", data[i]);
      printf("]\n");
    } else if (ndim == 2) {
      if (comment) printf("%s%s\n", blank, comment);
      else         printf("%sDist<> (%dx%dx%d, %d dimensional)\n", blank, size[0], size[1], size[2], countDimensions());
      for (i = 0; i < size[0]; i++) {
	printf("%s[ ", blank);
	for (j = 0; j < size[1]; j++) printf("%.2f ", data[ Index(i,j) ]);
	printf("]\n");
      }
    }
  }
  
  bool checkProbilityValidity(T tolerance = 0.0001) {
    int i;  T sum;
    for (sum = i = 0; i < total; i++) {
      if (!(data[i] >= 0.0 && data[i] <= 1.0)) {
	printf("Debug (MLH::Dist): invalid prob value %g at %d\n", data[i], i); 
	return false;
      }
      sum += data[i];
    }
    if (!(sum >= 1.0 - tolerance && sum <= 1.0 + tolerance)) {
      printf("Debug (MLH::Dist): invalid prob sum %g at %d\n", sum, i);  
      return false;
    }
    return true;
  }
  
};  
	    
	    
}	// end of namespace MLH

#endif  // MLH_DIST_HPP

