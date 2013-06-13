
//
// MTX::SMatrix<> class template
//
// Jaeil Choi
// last modified in Nov, 2006
//
// This sparse matrix implementation is based on 
// the linked list of entries with row/col index.
// 


#ifndef MTX_SMATRIX_NEW_HPP
#define MTX_SMATRIX_NEW_HPP

#include <iostream>
#include <iomanip>
#include "mtx_matrix.hpp"

namespace MTX {
  
using namespace std;

//static bool debug_sparse_matrix = false;

template <class T>
class SMatrixEntry {
public:
  SMatrixEntry<T> *next;	// pointer to the next entry in the row
  int row, col;
  T   value;
public:
  SMatrixEntry() { row=col=0; value=0; next=NULL; }
  SMatrixEntry(int r, int c, T v=0) { row=r; col=c; value=v; next=NULL; }
  ~SMatrixEntry() {}
};

// ===================================================================
// SMatrix
// ===================================================================

template <class T>
class SMatrix {
public:
  int   nRow, nCol;		// size of the matrix
  SMatrixEntry<T> **Row;
  int   total;			// total number of nonzero entries
private:
  T     dummy;

public:
  SMatrix() {  clear(false); }
  SMatrix(int row, int col) {  clear(false); set(row, col); }
  SMatrix(int row) {  clear(false); setVector(row); }
  SMatrix(Matrix<T> &m) { convertFromMatrix(m); }
  ~SMatrix() { clear(); }
  
  void clear(bool bFree=true) {
    nRow = nCol = total = 0;   
    if (bFree && Row) {
      SMatrixEntry<T> *ep, *next;
      for (int row = 0; row < nRow; row++)
	for (ep = Row[row]; ep; ep = next) {
	  next = ep->next;
	  delete(ep);
	}
      free(Row);
    }
    Row = NULL;
  }
  
  // Creating the matrix ---------------------------------------------
public:
  void set(int row, int col) {
    clear();
    Row = (SMatrixEntry<T>**) calloc( row, sizeof(SMatrixEntry<T>*) );
    nRow = row;  nCol = col;
    total = 0;
  }
  
  void setVector(int row, int col=1) {
    clear();
    int  i, j;
    if        (col == 1) {
      Row = (SMatrixEntry<T>**) calloc( row, sizeof(SMatrixEntry<T>*) );
      for (i = 0; i < row; i++) Row[i] = new SMatrixEntry<T>(row, 0, 0.0); 
      nRow = row;  nCol = col;  total = nRow * nCol;
    } else if (row == 1) {
      Row = (SMatrixEntry<T>**) calloc( 1, sizeof(SMatrixEntry<T>*) );
      for (j = 0; j < col; j++) {
	SMatrixEntry<T> *ep = new SMatrixEntry<T>(0, col-1-j, 0.0); 
	ep->next = Row[0];
	Row[0] = ep;
      }
      nRow = row;  nCol = col;  total = nRow * nCol;
    } else {
      cerr << "Error: SMatrix::setVector() is only for vectors" << endl;
    }
  }
  
  // entry -----------------------------------------------------------
public:
  T& operator()(int row, int col) { return getEntry(row, col)->value; }
  
  T& operator()(int idx) {
    if      (nCol == 1) return (*this)(idx, 0);  // column vector
    else if (nRow == 1) return (*this)(0, idx);  // row vector
    else return dummy;
  }
  
private:
  SMatrixEntry<T>* getEntry(int row, int col) {
    // If not found, create the entry.
    SMatrixEntry<T> *ep, *lastp=NULL;
    ep = findEntry(row, col, &lastp);
    if (ep == NULL) {	// create the entry
      ep = new SMatrixEntry<T>(row, col);
      if (lastp) {
	ep->next = lastp->next;
	lastp->next = ep;
      } else {
	ep->next = Row[row];
	Row[row] = ep;
      }
      total++;
    }
    return ep;
  }
  SMatrixEntry<T>* findEntry(int row, int col, SMatrixEntry<T> **lastp=NULL) {
    // If not found, return NULL.
    if (row < 0 || row >= nRow || col < 0 || col >= nCol) return NULL;
    SMatrixEntry<T> *ep = NULL, *lp = NULL;
    if (lastp == NULL) lastp = &lp;
    *lastp = NULL;
    for (ep = Row[row]; ep; ep = ep->next) {
      if      (ep->col == col) return ep;
      else if (ep->col  > col) return NULL;
      else *lastp = ep;
    }
    return NULL;
  }
  
  // Using the matrix ------------------------------------------------
public:
  void assign(int row, int col, ... ) {
    // IMPORTANT!! - All values must be given as double types
    //   Incorrect! : M.assign(3, 3,  11,12,0, 21,22,0,  0,0,1);
    //   Correct    : M.assign(3, 3,  11.0, 12.0, 0.0,  21.0, 22.0, 0.0,  0.0, 0.0, 1.0); 
    set(row, col);
    va_list valist;
    va_start( valist, col );
    for (int i = 0; i < row*col; i++) {
      T data = (T)va_arg(valist, double);
      if (data != 0) getEntry(i/col, i%col)->value = data;
    }
    va_end( valist );
  }
  
  void assignZero(void) {
    // remove all the entries
    SMatrixEntry<T> *ep, *next;
    for (int row = 0; row < nRow; row++) {
      for (ep = Row[row]; ep; ep = next) {
	next = ep->next;
	delete(ep);
      }
      Row[row] = NULL;
    }
    total = 0;
  }
  
  bool isEmpty(int i, int j) {
    // check if an entry exists, without creating it.
    return ( findEntry(i, j) ? false : true );
  }
  
  void removeZeros(void) {
    SMatrixEntry<T> *ep, *prev, *next;
    int  removed=0;
    for (int row = 0; row < nRow; row++) {
      prev = NULL;
      for (ep = Row[row]; ep; ep = next) {
	next = ep->next;
	if (ep->value == 0) {
	  if (prev) { prev->next = next; delete(ep); }
	  else      { Row[row]   = next; delete(ep); }
	  removed++;
	} else prev = ep;
      }
    }
    total -= removed;
  }
  
  void transpose(void) {
    SMatrix<T> m(nCol, nRow);
    // export all the entries to the temporary matrix
    SMatrixEntry<T> *ep, *last=NULL, *next=NULL;
    int  row, temp;
    for (row = 0; row < nRow; row++) {
      for (ep = Row[row]; ep; ep = next) {
	next = ep->next; 
	// take the entry off
	temp = ep->row;  ep->row = ep->col;  ep->col = temp;
	ep->next = NULL;
	// put the entry in new matrix
	m.findEntry(ep->row, ep->col, &last);
	if (last) { ep->next = last->next; last->next = ep; }
	else { ep->next = m.Row[ep->row]; m.Row[ep->row] = ep; }
      }
      Row[row] = NULL;
    }
    free(Row);
    // import all the entries back
    Row   = m.Row;   m.Row  = NULL;
    nRow  = m.nRow;  m.nRow = 0;
    nCol  = m.nCol;  m.nCol = 0;
  }
  
  bool isSymmetric(void) {
    SMatrixEntry<T> *ep;
    for (int row = 0; row < nRow; row++)
      for (ep = Row[row]; ep; ep = ep->next)
	if ( findEntry(ep->col, ep->row) == NULL) return false;
    return true;
  }
  
  bool isPositiveDefinite(void) {
    //// NOT IMPLEMENTED YET
    return false;
  }
  
  void operator= (SMatrix<T> &m) {
    set( m.nRow, m.nCol );
    SMatrixEntry<T> *ep;
    for (int row = 0; row < m.nRow; row++)
      for (ep = m.Row[row]; ep; ep = ep->next)
	(*this)(ep->row, ep->col) = ep->value;
    total = m.total;
  }
  
  bool operator== (SMatrix<T> &m) {
    if (nRow != m.nRow || nCol != m.nCol || total != m.total) return false;
    SMatrixEntry<T> *ep1, *ep2;
    for (int row = 0; row < nRow; row++) {
      ep1 = Row[row];
      ep2 = m.Row[row];
      while (ep1 || ep2) {
	if (!ep1 || !ep2 || ep1->col != ep2->col || ep1->value != ep2->value) return false;
	ep1 = ep1->next;
	ep2 = ep2->next;
      }
    }
    return true;
  }
  
  bool operator== (Matrix<T> &m) { 
    if (nRow != m.nRow || nCol != m.nCol) return false;
    for (int row = 0; row < nRow; row++)
      for (int col = 0; col < nCol; col++)
	if (isEmpty(row, col)) {
	  if (m(row, col) != 0) return false;
	} else {
	  if ((*this)(row, col) != m(row, col)) return false;
	}
    return true;
  }
  
  // conversion -----------------------------------------------------
public:
  void convertFromMatrix (Matrix<T>& m) {
    set( m.nRow, m.nCol );
    SMatrixEntry<T> *ep;
    int row, col, pos;
    for (row = 0; row < nRow; row++) {
      pos = row * nCol;
      for (col = nCol-1; col >= 0; col--) {
	if (m(row, col) == 0) continue;
	ep = new SMatrixEntry<T>(row, col);
	ep->value = m.data[ pos + col ];
	ep->next = Row[row];
	Row[row] = ep;
	total++;
      }
    }
    //printInfo(true);
  }
  
  void convertToMatrix (Matrix<T>& m) {
    m.set(nRow, nCol);
    SMatrixEntry<T> *ep;
    for (int row = 0; row < nRow; row++)
      for (ep = Row[row]; ep; ep = ep->next) 
	m(ep->row, ep->col) = ep->value;
  }
  
  void copyToMatrix (Matrix<T>& m, int row, int col) {
    T value;
    for (int i = 0; i < nRow; i++) {
      SMatrixEntry<T> *ep = Row[i];
      for (int j = 0; j < nCol; j++) {
	if (ep && ep->col == j) { value = ep->value;  ep = ep->next; }
	else  value = 0;
	m(row + i, col + j) = value;
      }
    }
  }

  // ----------------------------------------------------------------
  // algebraic operations
  // ----------------------------------------------------------------
public:
  inline SMatrix<T>& add( SMatrix<T> &A, SMatrix<T> &B) {
    return addWithScale( A, (T)(+1.0), B );
  }
  
  inline SMatrix<T>& sub( SMatrix<T> &A, SMatrix<T> &B) {
    return addWithScale( A, (T)(-1.0), B );
  }
  
  SMatrix<T>& addWithScale( SMatrix<T> &A, T s, SMatrix<T> &B) {
    assert( A.nRow == B.nRow && A.nCol == B.nCol );
    SMatrixEntry<T> *ea, *eb;
    int  row, col, type;
    T    value;
    if (this == &A) {		// as in "A.addWithScale( A, 0.1, B);"
      for (row = 0; row < nRow; row++)
	for (eb = B.Row[row]; eb; eb = eb->next)
	  getEntry(row, eb->col)->value += s * eb->value;
    } else if (this == &B) {	// as in "A.addWithScale( B, 0.1, A);"
      multValue(s);
      for (row = 0; row < nRow; row++)
	for (ea = A.Row[row]; ea; ea = ea->next)
	  getEntry(row, ea->col)->value += ea->value;
    } else {			// as in "A.addWithScale( B, 0.1, C);"
      set( A.nRow, A.nCol );
      for (row = 0; row < nRow; row++) {
	ea = A.Row[row];
	eb = B.Row[row];
	while (ea || eb) {
	  if      (ea == NULL) type = 2;
	  else if (eb == NULL) type = 1;
	  else if (ea->col < eb->col) type = 1;
	  else if (ea->col > eb->col) type = 2;
	  else type = 0;
	  switch (type) {
	  case 0:  col = ea->col;  value = ea->value + s * eb->value;  ea = ea->next;  eb = eb->next; break;
	  case 1:  col = ea->col;  value = ea->value;  ea = ea->next;  break;
	  case 2:  col = eb->col;  value = s * eb->value;  eb = eb->next;  break;
	  }
	  if (value != 0) getEntry(row, col)->value = value;
	}
      }
    }
    return *this;
  }
  
  SMatrix<T>& multValue(T v) {
    SMatrixEntry<T> *ep;
    for (int row = 0; row < nRow; row++)
      for (ep = Row[row]; ep; ep = ep->next)  ep->value *= v;
    return *this;
  }
  
  SMatrix<T>& mult( SMatrix<T> &m1, SMatrix<T> &m2 ) {
    if (m1.nCol != m2.nRow) { cerr << "Error (SMatrix::Mult): Invalid size" << endl; return *this; }
    set(m1.nRow, m2.nCol);
    int row, col;
    T   sum;
    SMatrixEntry<T> *ep1, *ep2;
    for (row = 0; row < nRow; row++)
      for (col = 0; col < nCol; col++) {
	// multiply the row of m1 and the column of m2
	sum = 0;
	for (ep1 = m1.Row[row]; ep1; ep1 = ep1->next) {
	  if ((ep2 = m2.findEntry( ep1->col, col )) == NULL) continue;
	  if (ep2) sum += ep1->value * ep2->value;
	}
	if (sum != 0) (*this)(row, col) = sum;
      }
    return *this;
  }
  
  SMatrix<T>& mult( Matrix<T> &m1, SMatrix<T> &m2 ) {
    if (m1.nCol != m2.nRow) { cerr << "Error (SMatrix::Mult): Invalid size" << endl; return *this; }
    set(m1.nRow, m2.Col);
    int row, col, j;
    T   sum;
    SMatrixEntry<T> *ep2;
    for (row = 0; row < nRow; row++)
      for (col = 0; col < nCol; col++) {
	// multiply the row of m1 and the column of m2
	sum = 0;
	for (j = 0; j < m1.nCol; j++) {
	  for (ep2 = m2.Row[j]; ep2 && ep2->col <= col; ep2 = ep2->next)
	    if (ep2->col == col) sum += m1(row, j) * ep2->value;
	}
	if (sum != 0) (*this)(row, col) = sum;
      }
    return *this;
  }
  
  void MultVector( T v[], T r[] ) {	// r[] = (*this) * v[]
    // multiply nCol x 1 vector 'v[]' to the matrix, and save the result in 'r[]'
    // Note that 'v[]' and 'r[]' should be different arrays.
    SMatrixEntry<T> *ep;
    for (int row = 0; row < nRow; row++) {
      r[row] = 0;
      for (ep = Row[row]; ep; ep = ep->next) 
	r[row] += ep->value * v[ep->col];
    }
  }
  
  T multVMV( SMatrix<T> &V ) {
    // calculates and return  Vt * M * V,  where 'V' is a vector
    if (V.nCol != 1) { cerr << "Error (SMatrix::multVMV): Invalid vector" << endl; return 0; }
    if (nRow != nCol || nCol != V.nRow) { cerr << "Error (SMatrix::multVMV): Invalid size" << endl; return 0; }
    SMatrixEntry<T> *ep;
    T result = 0, temp;
    for (int row = 0; row < nRow; row++) {
      temp = 0;
      for (ep = Row[row]; ep; ep = ep->next)
	temp += ep->value * V(ep->col);
      result += temp * V(row);
    }
    return result;
  }
  
  T multVMV( T v[] ) {
    // calculates and return  Vt * M * V,  where 'V' is a vector
    if (nRow != nCol) { cerr << "Error (SMatrix::multVMV): Invalid size" << endl; return 0; }
    SMatrixEntry<T> *ep;
    T result = 0, temp;
    for (int row = 0; row < nRow; row++) {
      temp = 0;
      for (ep = Row[row]; ep; ep = ep->next)
	temp += ep->value * v[ep->col];
      result += temp * v[row];
    }
    return result;
  }
  
  T dot( SMatrix<T> &A) {
    if (  nCol != 1) { cerr << "Error (SMatrix::Dot): Invalid vector size" << endl; return 0; }
    if (A.nCol != 1) { cerr << "Error (SMatrix::Dot): Invalid vector size" << endl; return 0; }
    T sum = 0;
    for (int row = 0; row < nRow; row++)  
      if (Row[row] && A.Row[row]) sum += Row[row]->value * A.Row[row]->value;
    return sum;
  }
  
  
  // ----------------------------------------------------------------
  // block copying and adding
  // ----------------------------------------------------------------
public:
  void CopyBlock(T b[], int size, int row, int col) {
    int i, j, pos;
    SMatrixEntry<T> *ep;
    if (nRow == 1 || nCol == 1) {	// vector
      ep = getEntry(row, col);
      for (i = 0; i < size; i++) { ep->value = b[i]; if (i<size-1) ep = getNextEntry(ep); }
    } else {				// matrix
      for (i = pos = 0; i < size; i++) {
	ep = getEntry(row + i, col);
	for (j = 0; j < size; j++) { ep->value = b[pos++]; if (j<size-1) ep = getNextEntry(ep); }
      }
    }
  }
  void addBlock(T b[], int size, int row, int col) {
    int i, j, pos;
    SMatrixEntry<T> *ep;
    if (nRow == 1 || nCol == 1) {	// vector
      ep = getEntry(row, col);
      for (i = 0; i < size; i++) { ep->value += b[i]; if (i<size-1) ep = getNextEntry(ep); }
    } else {				// matrix
      for (i = pos = 0; i < size; i++) {
	ep = getEntry(row + i, col);
	for (j = 0; j < size; j++) { ep->value += b[pos++]; if (j<size-1) ep = getNextEntry(ep); }
      }
    }
  }
private:
  SMatrixEntry<T>* getNextEntry(SMatrixEntry<T>* ep) {
    // find the next entry with same row and next column index.
    // create the entry if necessary.
    if (nCol == 1) {		// column vector
      int row2 = ep->row+1;
      if (!Row[row2]) { Row[row2] = new SMatrixEntry<T>(row2, 0); total++; }
      return Row[row2];
    } else {			// matrix or row vector
      if (!ep->next || ep->next->col != ep->col+1) {
	SMatrixEntry<T> *np = new SMatrixEntry<T>(ep->row, ep->col+1);
	np->next = ep->next;
	ep->next = np;
	total++;
      }
      return ep->next;
    }
  }
  
  
  // ----------------------------------------------------------------
  // other functions
  // ----------------------------------------------------------------
public:  
  void printInfo(bool detail=false) {
    printf("SparseMatrix ( %d x %d :  total = %d )\n", nRow, nCol, total);
    if (nRow > 0 && Row == NULL) printf("  Error: Row[] array shouldn't be NULL\n");
    
    // check validity
    int row, last_col, count = 0;
    SMatrixEntry<T> *ep;
    for (row = 0; row < nRow; row++) {
      last_col = -1;
      if (detail) cout << "  row " << row << " : ";
      for (ep = Row[row]; ep; ep = ep->next) {
	if (ep->row != row || ep->col <= last_col || ep->col >= nCol) {
	  printf("  Error: Invalid row/col (%d %d) order at row %d\n", ep->row, ep->col, row);
	  last_col = ep->col;
	}
	if (detail) cout << "(" << ep->row << "," << ep->col << "; " << ep->value << ")  ";
	count++;
      }
      if (detail) cout << endl;
    }
    if (total != count)
      printf("  Error: Invalid total number (%d), %d entries actually found\n", total, count);
  }
  
  // -----------------------------------------------------------------
  // file I/O
  // -----------------------------------------------------------------
public:
  void writeFile(char *filename) {
    ofstream ofs(filename);
    ofs << "SparseMatrix  " << nRow << " rows  " << nCol << " cols" << "  (" << total << ")" << endl;
    int row, col;
    SMatrixEntry<T> *ep;
    for (row = 0; row < nRow; row++) {
      ep = Row[row];
      for (col = 0; col < nRow; col++) {
	if (ep && ep->col == col) {
	  if (ep->value == 0) ofs << "\t0\t"; 
	  else            ofs << setw(12) << ep->value << "\t";
	  ep = ep->next;
	} else {
	  ofs << "\t-\t";
	}
      }
      ofs << "\n";
    }
  }
  bool checkFile(char *filename) {
    FILE *fp = fopen(filename, "r");
    if (fp == NULL) return false;
    char str[256];  int row, col, total;
    int  n = fscanf(fp, "%s %d rows %d cols  (%d)\n", str, &row, &col, &total);
    fclose(fp);
    return (n == 4 && strcmp(str, "SparseMatrix") == 0);
  }
  bool readFile(char *filename) {
    FILE *fp = fopen(filename, "r");
    if (fp == NULL) return false;
    int row, col, total, count = 0;
    char str[256];  
    int  n = fscanf(fp, "%s %d rows %d cols  (%d)\n", str, &row, &col, &total);
    if (n != 4 || strcmp(str, "SparseMatrix") != 0) return false;
    set( row, col );
    for (int i = 0; i < nRow; i++)
      for (int j = 0; j < nCol; j++) {
	fscanf(fp, "%s", str);
	if (str[0] == '-' && str[1] == '\0' || str[0] == 'N') continue;
	(*this)(i, j) = (T)atof(str);
	count++;
      }
    fclose(fp);
    return true;
  }
};


// ===================================================================
// Output
// ===================================================================

template <class T>
static ostream &operator<< (ostream &os, SMatrix<T> &m)
{
  int row, col, j;
  if (m.nRow == 0 || m.nCol == 0) {
    cout << "[ ] (" << m.nRow << "x" << m.nCol << ")";	// [ ] (0x0)
  } else if (m.nCol == 1) {		// column vector
    os << "[ ";
    for (row = 0; row < m.nRow; row++)
      os << setw(12) << m.Row[row]->value << " ";
    os << "]T ";
  } else if (m.nRow == 1) {		// row vector
    os << "[ ";
    SMatrixEntry<T> *ep = m.Row[0];
    for (col = 0; col < m.nCol; col++) {
      if (ep && col == ep->col) {
	os << setw(12) << ep->value << " ";
	ep = ep->next;
      } else os << setw(12) << '-' << " ";
    }
    os << "]  ";
  } else {				// matrix
    for (row = 0; row < m.nRow; row++) {
      os << (row == 0 ? "[ " : "  ");
      SMatrixEntry<T> *ep = m.Row[row];
      for (col = 0; col < m.nCol; col++) {
	if (ep && col == ep->col) {
	  os << setw(12) << ep->value << " ";
	  ep = ep->next;
	} else
	  os << setw(12) << " " << " ";
      }
      os << (row == m.nRow-1 ? "]  " : "\n");
    }
  }
  return os;
}

template <class T>
void PrintMatrices(SMatrix<T> &m1, SMatrix<T> &m2, char *comment = NULL)
{
  int row, i, j;
  row = (m1.nRow > m2.nRow) ? m1.nRow : m2.nRow;
  if (comment)  cout << "SparseMatrices: " << comment << endl;
  else          cout << "SparseMatrices: " << endl;
  for (i = 0; i < row; i++) {
    
    // m1
    if (i == 0 && m1.nRow == 0 && m1.nCol == 0) cout << "[]";
    else if (i == 0) cout << "[ ";
    else             cout << "  ";
    if (i < m1.nRow) {
      SMatrixEntry<T> *ep = m1.Row[i];
      for (j = 0; j < m1.nCol; j++) {
	if (j > 0) cout << " ";
	if (ep && j == ep->col) {
	  cout << setw(12) << ep->value;
	  ep = ep->next;
	} else cout << setw(12) << " ";
      }
    } else cout << setw(12) << " ";
    if (i == 0)              cout << "    ";
    else if (i == m1.nRow-1) cout << " ]  ";
    else                     cout << "    ";

    // m2
    if (i == 0 && m2.nRow == 0 && m2.nCol == 0) cout << "[]";
    else if (i == 0) cout << "[ ";
    else             cout << "  ";
    if (i < m2.nRow) {
      SMatrixEntry<T> *ep = m2.Row[i];
      for (j = 0; j < m2.nCol; j++) {
	if (j > 0) cout << " ";
	if (ep && j == ep->col) {
	  cout << setw(12) << ep->value;
	  ep = ep->next;
	} else cout << setw(12) << " ";
      }
    } else cout << setw(12) << " ";
    if (i == 0)              cout << "    ";
    else if (i == m2.nRow-1) cout << " ]  ";
    else                     cout << "    ";
    
    cout << endl;
  }
}

template <class T>
void PrintMatrices(SMatrix<T> &m1, SMatrix<T> &m2,
		   SMatrix<T> &m3, char *comment = NULL)
{
  int row, i, j;
  row = ( (m1.nRow > m2.nRow) ?
	  (m1.nRow > m3.nRow ? m1.nRow : m3.nRow) :
	  (m2.nRow > m3.nRow ? m2.nRow : m3.nRow) );
  if (comment)  cout << "SparseMatrices: " << comment << endl;
  else          cout << "SparseMatrices: " << endl;
  for (i = 0; i < row; i++) {
    
    // m1
    if (i == 0 && m1.nRow == 0 && m1.nCol == 0) cout << "[]";
    else if (i == 0) cout << "[ ";
    else             cout << "  ";
    if (i < m1.nRow) {
      SMatrixEntry<T> *ep = m1.Row[i];
      for (j = 0; j < m1.nCol; j++) {
	if (j > 0) cout << " ";
	if (ep && j == ep->col) {
	  cout << setw(12) << ep->value;
	  ep = ep->next;
	} else cout << setw(12) << " ";
      }
    } else cout << setw(12) << " ";
    if (i == 0)              cout << "    ";
    else if (i == m1.nRow-1) cout << " ]  ";
    else                     cout << "    ";

    // m2
    if (i == 0 && m2.nRow == 0 && m2.nCol == 0) cout << "[]";
    else if (i == 0) cout << "[ ";
    else             cout << "  ";
    if (i < m2.nRow) {
      SMatrixEntry<T> *ep = m2.Row[i];
      for (j = 0; j < m2.nCol; j++) {
	if (j > 0) cout << " ";
	if (ep && j == ep->col) {
	  cout << setw(12) << ep->value;
	  ep = ep->next;
	} else cout << setw(12) << " ";
      }
    } else cout << setw(12) << " ";
    if (i == 0)              cout << "    ";
    else if (i == m2.nRow-1) cout << " ]  ";
    else                     cout << "    ";
    
    // m3
    if (i == 0 && m3.nRow == 0 && m3.nCol == 0) cout << "[]";
    else if (i == 0) cout << "[ ";
    else             cout << "  ";
    if (i < m3.nRow) {
      SMatrixEntry<T> *ep = m3.Row[i];
      for (j = 0; j < m3.nCol; j++) {
	if (j > 0) cout << " ";
	if (ep && j == ep->col) {
	  cout << setw(12) << ep->value;
	  ep = ep->next;
	} else cout << setw(12) << " ";
      }
    } else cout << setw(12) << " ";
    if (i == 0)              cout << "    ";
    else if (i == m3.nRow-1) cout << " ]  ";
    else                     cout << "    ";

    cout << endl;
  }
}


}	// end of namespace MTX


#endif  // MTX_SMATRIX_NEW_HPP


