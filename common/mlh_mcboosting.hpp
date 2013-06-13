
//
// MLH::MCBoosting<X_t,H_t> class template for multi-class classification
//
// Jaeil Choi
// last modified in July, 2009
//
// Template types
//   X_t : type of training input data (float or double)
//   H_t : type of the weak learner (for example: 'MLH::DTree<X_t>')
//         The weak learner must implement following functions:
//         - void  clear(void);
//         - void  setAttributeInfo( AttributeInfo<> *ap );
//         - void  setTrainingData ( TrainingData<>  *tp );
//         - bool  train( int label, int wlopt, int dlist[]=NULL, int dsize=0 );
//         - bool  test ( X_t x[], double prob[2]=NULL, double *confidence=NULL );
//         - int   testOnTData(char *cmmt=NULL, int lst[]=NULL, int nlst=0);
//         - void  writeFileSingleLine(FILE *fp);
//         - bool  readFileSingleLine(FILE *fp);
//         - char* getSingleLineInfo(char *buffer);
// Member functions
//   - void  setAttributeInfo( AttributeInfo<> *ap );
//   - void  setTrainingData ( TrainingData<>  *tp );
//   - void  train(int nh, int wlopt, int dlist[]=NLL, int dsize=0);
//   - int   test (X_t x[], double prob[]=NULL);
//   - double testForProb(X_t x[], int lbl);
//   - int   testOnTData(char *cmmt=NULL, int lst[]=NULL, int nlst=0);
//   - bool  writeFile(const char *fname, const char *cmmt=NULL);
//   - bool  readFile(const char *fname);
// 

#ifndef MLH_MCBOOSTING_HPP
#define MLH_MCBOOSTING_HPP

#include <iostream>
#include <fstream>
#include <cmath>
#include "util_word_breaker.hpp"
#include "util_clock.hpp"
#include "mlh_pdf.hpp"
#include "mlh_attribute_info.hpp"
#include "mlh_training_data.hpp"
#include "mlh_boosting.hpp"

  
namespace MLH {
  
template <class X_t, class H_t>
class MCBoosting {

  // prerequisite information ----------------------------------------
public:
  AttributeInfo<X_t>	*ai;		// attribute information
  TrainingData<X_t>	*tt;		// training data
private:
  AttributeInfo<X_t>	ainfo;		// attribute information (local)
  TrainingData<X_t>	tdata;		// training data dummy
  // experts (weak learners) -----------------------------------------
public:
  int			nc;		// number of classes (ai->getLabelCCount())
  Boosting<X_t, H_t>	*mcc;		// multi-class classifiers
  // etc -------------------------------------------------------------
public:
  bool			verbose;
  
public:
  MCBoosting() : ai(&ainfo), tt(&tdata), nc(0), mcc(NULL), verbose(false) {}
  ~MCBoosting() { clear(); }
  void clear(void) {
    if (mcc) {
      for (int i=0; i<nc; i++) mcc[i].clear();
      // delete(mcc);  mcc = NULL; ////
      free(mcc);  mcc = NULL;
    }
    nc = 0;
  }
  inline bool  isReady(void) { return (mcc && mcc[0].isReady()); }
  
  // -----------------------------------------------------------------
  // Setting Data Information
  // -----------------------------------------------------------------
public:	// use the data information given by the user
  void setAttributeInfo(AttributeInfo<X_t> *ainfo=NULL) {
    // Use provided attribute information. (It might be shared with others)
    this->ai = ( ainfo ? ainfo : &this->ainfo );
    if (!ai->isReady()) printf("Warning (MCBoosting::setAttributeInfo): AttributeInfo is not ready yet\n");
  }
  void setTrainingData(TrainingData<X_t> *tp=NULL) {
    // Use provided training data. (It might be shared with others)
    this->tt = ( tp ? tp : &this->tdata );
    if (!tt->isReady()) printf("Warning (MCBoosting::setTrainingData): TrainingData is not ready yet\n");
  }
  
  // -----------------------------------------------------------------
  // Training
  // -----------------------------------------------------------------
  
public:
  bool train(int nh, int wlopt, int dlist[]=NULL, int dsize=0) {
    // Train binary classifier with given list of training data.
    //   nh      : number of iterations (hypotheses)
    //   wlopt   : weak learner option (max. depth for decision tree, example)
    //   dlist[] : list of indices to be used for training
    //   dsize   : number of data  to be used for training
    if (!ai || !ai->isReady())  { 
      if (verbose) fprintf(stderr, "[MCBst] Error: call setAttributeInfo() first.\n"); 
      return false; 
    }
    if (!tt || !tt->isReady())  { 
      if (verbose) fprintf(stderr, "[MCBst] Error: call setTrainingData() first.\n"); 
      return false; 
    }
    clear();
    this->nc = ai->getLabelCCount();
    // this->mcc = new Boosting<X_t, H_t>[ nc ];  ////
    this->mcc = (Boosting<X_t, H_t>*) calloc ( nc, sizeof(Boosting<X_t,H_t>) );
    for (int c = 0; c < nc; c++) {
      if (verbose) printf("[MCBst] Training binary classifier for class %d\n", c);
      mcc[c].setAttributeInfo( ai );
      mcc[c].setTrainingData( tt );
      mcc[c].train( c, nh, wlopt, dlist, dsize );
    }
    return true;
  }
  
  // -----------------------------------------------------------------
  // Testing
  // -----------------------------------------------------------------
public:
  double testForProb(X_t x[], int lbl) {
    // Find the probability of the label 'lbl' for the data 'x[]'.
    double prob[128];  test( x, prob );
    return prob[lbl];
  }
  int test(X_t x[], double prob[]=NULL) {
    // Calculate the probability for each class, and return the best label (class).
    int    c, maxc=0;
    double pr[2]={0,0}, psum=0, maxp=0, pvalues[128];
    if (!prob) prob = pvalues;
    for (c=0; c<nc; c++) {
      mcc[c].test( x, pr );
      prob[c] = pr[1];
      psum   += pr[1];
      if (pr[1] > maxp) { maxp = pr[1]; maxc = c; }
    }
    for (c=0; c<nc; c++)  prob[c] /= psum;
    return maxc;
  }
  int testOnTData(char *cmmt=NULL, int lst[]=NULL, int nlst=0) {
    // Test on a subset of the training data (for debugging)
    if (!tt || !tt->isReady()) return 0;
    if (nlst <= 0 || !lst) { nlst = tt->nt; lst = NULL; }
    printf("%s(%d) : <", (cmmt ? cmmt : "TestOnDataset"), nlst);
    int  i, cnt=0, tidx;
    double wsum=0, prob[128], psum=0;
    for (i = 0; i < nlst; i++) {
      tidx = ( lst ? lst[i] : i );
      int  ret  = test( tt->tx[tidx], prob );
      bool succ = (ret == tt->ty[tidx]);
      if (nlst<=5000) { 
	printf("%c", (succ ? '.':('0'+tt->ty[tidx])));
	if (i%10==9) printf(" ");
      }
      psum += prob[ret];
      if (succ) { cnt++;  if (tt->w) wsum += tt->w[tidx]; }
    }
    if (nlst>5000) printf(" TooManyToShow ");
    if (tt->w) printf(">  %d%% success (%.0f%% w)  avgPr=%.4f\n", cnt*100/nlst, wsum*100, psum/nlst);
    else       printf(">  %d%% success  avgPr=%.4f\n", cnt*100/nlst, psum/nlst);
    return cnt;
  }
  
  // -----------------------------------------------------------------
  // file I/O
  // -----------------------------------------------------------------
public:
  bool writeFile(const char *fname, const char *cmmt=NULL) {
    if (!ai || !ai->isReady()) { 
      if (verbose) fprintf(stderr, "[MCBst] Error: call setAttributeInfo() first.\n"); 
      return false; 
    }
    FILE *fp = fopen( fname, "w+" );      // open the file
    if (fp == NULL) { 
      if (verbose) fprintf(stderr, "Error (MLH::MCBoosting::writeFile): cannot open file '%s'\n", fname);
      return false;
    }
    bool ret = writeFile( fp, cmmt );
    fclose( fp );
    return ret;
  }
  bool writeFile(FILE *fp, const char *cmmt=NULL) {
    if (!fp) return false;
    if (!ai || !ai->isReady()) { 
      if (verbose) fprintf(stderr, "[MCBst] Error: call setAttributeInfo() first.\n"); 
      return false; 
    }
    if (!isReady()) { 
      if (verbose) fprintf(stderr, "[MCBst] Error: classifier not ready\n");
      return false; 
    }
    // write general information
    fprintf( fp, "# %s\n", (cmmt ? cmmt : "") );
    fprintf( fp, "# Multi-Class Boosting Trained Knowledge\n" );
    fprintf( fp, "# \n" );
    ai->writeAttributeInfo( fp );
    int last_nt = (nc>0 && mcc ? mcc[0].last_nt : 0);
    UTIL::Clock clock;
    if (last_nt > 0) 
      fprintf(fp, "Training  : trained on %d data (saved at %s)\n", last_nt, clock.getLocalTime(4) );
    else
      fprintf(fp, "Training  : training data unknown (saved at %s)\n", clock.getLocalTime(4) );
    // write hypotheses
    for (int c = 0; c < nc; c++) {
      fprintf( fp, "Class  %d  : %d experts  (%s==%s)\n", c, mcc[c].nh, ai->getLabelName(), ai->convI2S(-1,c) );
      for (int i=0; i < mcc[c].nh; i++) {
	fprintf( fp, "Ex%02d (%8.6f) : ", i, mcc[c].hw[i] );
	mcc[c].mh[i].writeFileSingleLine( fp );
      }
    }
    return true;
  }

  bool readFile(const char *fname) {
    // open the file
    FILE *fp = fopen( fname, "r" );
    if (fp == NULL) { 
      if (verbose) fprintf(stderr, "Error (MCBoosting::readFile): cannot open file '%s'\n", fname);
      return false; 
    }
    bool ret = readFile( fp );
    fclose(fp);
    return ret;
  }
  bool readFile(FILE *fp) {
    if (!fp) return false;
    int  i, idx, c, n, classifiers=0;
    char line[1024]; 
    clear();
    if (!ai || !ai->readAttributeInfo( fp )) { // read attribute information
      if (verbose) fprintf(stderr, "Error (MCBoosting::readFile): invalid AttributeInfo\n");
      return false;
    }
    setAttributeInfo();		// set  attribute info
    fgets(line, 250, fp);	// read training info
    this->nc = ai->getLabelCCount();
    // this->mcc = new Boosting<X_t, DTree<X_t> >[nc];  ////
    this->mcc = (Boosting<X_t, H_t>*) calloc ( nc, sizeof(Boosting<X_t,H_t>) );
    bool ret[128];
    UTIL::WordBreaker wb;
    while (fgets(line, 250, fp)) {
      wb.parse( line, " :()\t\r\n", false, '#' );
      if (wb.count() < 1) continue;  // skip empty lines
      if (strcmp(wb[0],"Class")==0) {	// -----------
	c = atoi(wb[1]);
	n = atoi(wb[2]);
	if (c < 0 || c >= ai->getLabelCCount() || n<1) break;
	mcc[c].setAttributeInfo( ai );
	mcc[c].initHypotheses( n );
	for (i=0; i < mcc[c].nh; i++) {
	  // read the weight of the hypothesis
	  if (fscanf( fp, "Ex%d (%lf) : ", &idx, mcc[c].hw+i ) != 2) break;
	  // read the hypothesis
	  if (! mcc[c].mh[i].readFileSingleLine( fp )) break;
	}
	for (mcc[c].hwsum=i=0; i < mcc[c].nh; i++) mcc[c].hwsum += mcc[c].hw[i];
	ret[c] = (i == mcc[c].nh);
	if (++classifiers == nc) break;
      } else break;
    }
    for (c = 0; c < nc; c++) if (!ret[c]) break;
    return (c==nc);
  }
  
  // -----------------------------------------------------------------
  // miscellaneous
  // -----------------------------------------------------------------
public:
  void printInfo(char *cmmt=NULL, bool single_line=false) {
    if (single_line) {
      if (!isReady()) printf("%s: NOT READY\n", (cmmt ? cmmt : "MCBoosting"));
      else {
	printf("%s: %2d attr, %d label => ( ", 
	       (cmmt ? cmmt : "MCBoosting"), ai->getAttrCount(), ai->getLabelCCount());
	for (int c=0; c<nc; c++) printf("%d ", mcc[c].nh);
	printf(") experts\n");
      }
    } else {
      char blank[80], str[256];  int indent, i, c;
      for (indent=0; indent<70; indent++) if (!cmmt || cmmt[indent]!=' ') break;
      memset(blank, ' ', indent);  blank[indent]='\0';
      printf("%s  with %d classes and { ", (cmmt ? cmmt : "MCBoosting"), nc);
      for (c = 0; c < nc; c++) printf("%d ", mcc[c].nh);
      printf("} experts\n");
      if (ai && ai->isReady()) ai->printInfo("  Att : ", true);
      if (tt && tt->isReady()) tt->printInfo("  TrD : ");
      for (c = 0; c < nc; c++) {
	printf("%s  Class %d : %d experts  (%s==%s)\n", blank, c, mcc[c].nh, ai->getLabelName(), ai->convI2S(-1,c));
	for (i = 0; i < mcc[c].nh; i++)
	  printf("%s  Ex%02d : a=%.4f  %s \n", blank, i, mcc[c].hw[i], mcc[c].mh[i].getSingleLineInfo(str));
      }
    }
  }
  
};


}	// end of namespace MLH

#endif // MLH_MCBOOSTING_HPP
