
//
// MLH::TrainingData<X_t>
//
// Jaeil Choi
// Last modified in Sep, 2009
//


#ifndef MLH_TRAINING_DATA_HPP
#define MLH_TRAINING_DATA_HPP

#include <iostream>
#include <cmath>
#include <cstdarg>
#include "mlh_attribute_info.hpp"

namespace MLH {
  
template <class X_t>
class TrainingData {
public:
  AttributeInfo<X_t> *ai;	// attribute info (pointer only)
  int		sx;		// size of each feature vector
  int		nt;		// number of feature vectors
  int		nl;		// number of labels
  // training data
  X_t		**tx;		// list of all data (feature vectors)
  int		*ty;		// list of labels for each data
  double	*w;		// weight for each data
public:
  TrainingData() : ai(NULL), sx(0), nt(0), nl(0), tx(NULL), ty(NULL), w(NULL) {}
  TrainingData(AttributeInfo<X_t> *a) : nt(0), tx(NULL), ty(NULL), w(NULL) { setAttributeInfo(a); }
  ~TrainingData() { clear(); }
  void clear(void) {
    if (tx) { for (int i=0; i<nt; i++) if (tx[i]) free( tx[i] ); }
    if (tx) free(tx);
    if (ty) free(ty);
    if (w)  free(w);
    tx = NULL;  ty = NULL;  w = NULL;
    nt = 0;
  }
  inline bool isReady(void) { return (ai && sx>0 && nt>0 && nl>0 && tx && ty); }
  
  // -----------------------------------------------------------------
  // training data values
  // -----------------------------------------------------------------
public:
  inline X_t& attr(int t, int a) { return tx[t][a]; }
  inline int& label(int t)       { return ty[t]; }
  
  void setAttributeInfo(AttributeInfo<X_t> *a) { 
    clear();
    this->ai = a;
    this->sx = a->getAttrCount();
    this->nl = a->getLabelCCount(); 
  }
  bool initTrainingData(int nt) {
    if (!ai) { fprintf(stderr, "Error (initTrainingData): AttributeInfo NOT set\n"); return false; }
    if (this->tx && this->ty && this->nt == nt) {
      for (int t=0; t < nt; t++) memset( tx[t], 0, sx*sizeof(X_t) );
      memset( ty, 0, nt*sizeof(X_t) );
      return true;
    }
    clear();
    tx = (X_t**)calloc(nt, sizeof(X_t*));
    ty = (int*) calloc(nt, sizeof(int) );
    if (!tx || !ty) { clear(); return false; }
    for (int t=0; t < nt; t++) {
      tx[t] = (X_t*)calloc( sx, sizeof(X_t) );
      if (!tx[t]) { clear(); return false; }
    }
    this->nt = nt;
    return true;
  }
  void copyTrainingDataFromBuffers(X_t **tx, int *ty, int nt) {
    initTrainingData( nt );
    for (int t = 0; t < nt; t++) {
      //for (int a = 0; a < sx; a++) this->tx[t][a] = tx[t][a];
      memcpy( this->tx[t], tx[t], sx*sizeof(X_t) );
      this->ty[t] = ty[t];
    }
    this->nl = nl;
  }
  
  // -----------------------------------------------------------------
  // weights
  // -----------------------------------------------------------------
public:
  void initTrainingWeights(void) {
    if (!isReady()) return;
    if (w) { free(w); w = NULL; }
    w = (double*)malloc(nt * sizeof(double));
    for (int i=0; i<nt; i++) w[i] = 1.0;
  }
  void setTrainingWeightsWithBuffer(double weights[]) {
    if (!isReady()) return;
    if (!w) initTrainingWeights();
    for (int i=0; i<nt; i++) w[i] = weights[i];
  }
  
  // -----------------------------------------------------------------
  // file I/O
  // -----------------------------------------------------------------
public:
  bool writeTrainingDataFile(const char *filename, char *cmmt=NULL) {
    if (!isReady()) return false;
    FILE *fp = fopen( filename, "w+" );
    if (!fp) { printf("Error (writeTrainingDataFile): cannot open file '%s'\n", filename); return false; }
    ai->writeAttributeInfo( fp );
    fprintf(fp, "# TrainingData count %d with feature vector size %d\n", nt, sx);
    fprintf(fp, "# idx label : ");
    int  i, j;
    for (j=0; j < sx; j++) fprintf(fp, "%12s ", ai->getAttrName(j));
    fprintf(fp, "\n");
    for (i=0; i < nt; i++) {
      fprintf(fp, "%5d %5d : ", i, ty[i]);
      for (j=0; j < sx; j++) {
	if (tx[i][j]==UNKNOWN) fprintf(fp, "%12s ", " ? ");
	else fprintf(fp, "%12g ", tx[i][j]);
      }
      fprintf(fp, "\n");
    }
    fclose(fp);
    return true;
  }
  
  bool readTrainingDataFile(const char *filename, bool use_name=false) {
    if (!ai) { fprintf(stderr, "Error (initTrainingData): AttributeInfo NOT set\n"); return false; }
    clear();
    FILE *fp = fopen( filename, "r" );
    if (!fp) { printf("Error (readTrainingDataFile): cannot open file '%s'\n", filename); return false; }
    bool ret = ai->readAttributeInfo(fp);
    setAttributeInfo( ai );
    if (!ret) return false;
    int  i, j, wc, blen, ntt=0;
    char buf[2048], *ws[400];
    X_t  **ttx=NULL;  int *tty=NULL;
    if (use_name) {	// data file with descriptive names ----------
      while (fgets(buf, 2040, fp)) {
	ttx = (X_t**)realloc(ttx, (ntt+1) * sizeof(X_t*));
	tty = (int*) realloc(tty, (ntt+1) * sizeof(int) );
	bool neww=true;
	blen = strlen(buf);
	for (i = wc = 0; i < blen; i++) {
	  if (buf[i]==',' || buf[i]=='\n') { neww = true; buf[i] = '\0'; }
	  else if (neww) { if (buf[i]==' ') i++; ws[wc++] = buf+i; neww=false; }
	}
	if (wc < sx || ntt >= 2048) continue;
	ttx[ntt] = (X_t*)calloc( sx, sizeof(X_t) );		// training case
	for (i = 0; i < sx; i++)
	  ttx[ntt][i] = (ai->getAttr(i)->discrete ? (X_t)ai->convS2I(i, ws[i]) : (X_t)atof(ws[i]));
	tty[ntt] = ai->convS2I(-1, ws[i]);
	ntt++;
      }
    } else {		// data file with numbers --------------------
      fgets(buf, 2040, fp);	// the first line
      if (sscanf(buf, "# TrainingData count %d", &ntt) != 1) {
	fprintf(stderr, "Error (AttributeInfo::readTraniningData): invalid file format\n");
	return false;
      }
      char value[1024];
      fgets(buf, 2040, fp);	// the second line
      ttx = (X_t**)calloc(ntt, sizeof(X_t*));
      tty = (int*) calloc(ntt, sizeof(int) );
      for (i=0; i<ntt; i++) {
	ttx[i] = (X_t*)calloc( sx, sizeof(X_t) );		// training case
	fscanf( fp, " %d %d :", &j, tty+i );
	for (j=0; j<sx; j++) {
	  if (fscanf( fp, " %s", value ) != 1) break;
	  if (value[0]=='?') ttx[i][j] = UNKNOWN;
	  else               ttx[i][j] = (X_t)atof(value);
	}
	if (j<sx) break;
      }
      ntt = i;
    }
    if (ntt <= 0) { fprintf(stderr, "Error (AttributeInfo::readTraniningData): no data at all\n"); return false; }
    copyTrainingDataFromBuffers( ttx, tty, ntt );
    fclose(fp);
    return true;
  }
  
  // -----------------------------------------------------------------
  // file I/O
  // -----------------------------------------------------------------
public:
  void printInfo(char *cmmt=NULL) {
//     if (!isReady())
//       printf("%s is NOT ready\n", (cmmt ? cmmt : "TrainingData"));
//     else 
      printf("%s (%d attributes %d labels) x %d %s weights\n", 
	     (cmmt ? cmmt : "TrainingData"), sx, nl, nt, (w ? "with":"without"));
  }
  
  void printTrainingData(char *cmmt=NULL, int tlst[]=NULL, int nlst=0) {
    // Print the training dataset.
    //   When 'tlst[]' and 'nlst' are given, print the selected subset only.
    if (!isReady()) return;
    printf("%s : %d data %s weights\n", (cmmt ? cmmt : "MLH::Training Data"), nt, (w ? "with":"without"));
    ai->printInfo("  Attrb: ", true);
    for (int i=0; i<nt; i++) {
      if (tlst && i >= nlst) break;
      int idx = (tlst ? tlst[i] : i);
      printf("  %03d:  ", idx);
      for (int a=0; a<sx; a++) {
	if (tx[idx][a]==UNKNOWN) printf("%5s ", " ? ");
	else printf("%5.2f ", tx[idx][a]);  // input attributes
      }
      printf("=>  %d  ", (int)ty[idx]);		// output label
      if (w) printf("(w=%.2f)", w[idx]);	// weight
      printf("\n");
    }
  }
};
  
}

#endif  // MLH_TRAINING_DATA_HPP

