
//
// MLH::AdaBoost<X_t,H_t> class template
//
// Jaeil Choi
// last modified in June, 2008
// 
// WARNING! : THIS CLASS IS NOT COMPLETE. DO NOT USE.
//
// ===================================================================
// This class is obsolete. Use MLH::Boosting<X,H> instead.
// ===================================================================
//
// Template types
//   X_t : type of training input data (float or double)
//   H_t : type of the weak learner (for example: 'MLH::DecisionTree<X_t>')
//         The weak learner must implement following functions:
//         - void  setAttributeInfo( AttributeInfo<> *ap );
//         - void  setTrainingData ( TrainingData<>  *tp );
//         - bool  train(int label, int wlopt, int dlist[]=NULL, int dsize=0);
//         - int   test(X_t x[], double *prob);
//         - char* getSingleLineInfo(char str[]);
//         - void  writeFileSingleLine(FILE *fp);
//         - bool  readFileSingleLine(FILE *fp);
// Member functions
//   - void setAttributeInfo( AttributeInfo<> *ap );
//   - void setTrainingData ( TrainingData<>  *tp );
//   - void train  (int label, int nh, int wlopt, int dlist[]=NLL, int dsize=0);
//   - void trainM1(int nh, int wlopt, int dlist[]=NLL, int dsize=0);
//   - int  test   (X_t *x, double *prob);
//   - int  testM1 (X_t *x, double *prob);
// 

#ifndef MLH_ADABOOST_HPP
#define MLH_ADABOOST_HPP

#include <iostream>
#include <fstream>
#include <cmath>
#include "util_word_breaker.hpp"
#include "util_clock.hpp"
#include "mlh_pdf.hpp"
#include "mlh_attribute_info.hpp"

// ===================================================================
// AdaBoost
// ===================================================================
  
namespace MLH {
  
template <class X_t, class H_t>
class AdaBoost {

  // prerequisite information ----------------------------------------
public:
  AttributeInfo<X_t>	*ai;		// attribute information
  TrainingData<X_t>	*tt;		// training data
private:
  AttributeInfo<X_t>	ainfo;		// attribute information (local)
  TrainingData<X_t>	tdata;		// training data dummy
  int			sx;		// size of each feature vector
  // experts (weak learners) -----------------------------------------
public:
  int		nh;		// number of hypotheses
  H_t		*h;		// pointers to experts (WeakLearners)
  double	*hw;		// weights for experts (WeakLearners)
  double	hwsum;
  // test result -----------------------------------------------------
public:
  double	*prv;
  // etc -------------------------------------------------------------
public:
  bool  	 verbose;
  int		last_nt;
  
public:
  AdaBoost() : ai(&ainfo), tt(&tdata), sx(0), nh(0), 
	       h(NULL), hw(NULL), prv(NULL), verbose(false), last_nt(0) {}
  ~AdaBoost() { clear(); }
  void clear(void) {
    for (int i=0; i<nh && h; i++) h[i].clear();
    if (h)  free(h);   h = NULL;     nh = 0;
    if (hw) free(hw);  hw = NULL;
  }
  inline bool  isReady(void) { return nh > 0; }
  inline int   getFeatureVectorSize(void) { return sx; }
  inline int   getAttributeCount(void) { return ai->getAttrCount(); }
  inline char* getLabelName(int lbl) { return ai->convI2S( -1, lbl ); }
  
  // -----------------------------------------------------------------
  // Setting Data Information
  // -----------------------------------------------------------------
public:	// use the data information given by the user
  void setAttributeInfo(AttributeInfo<X_t> *info=NULL) {
    this->ai = ( info ? info : &this->ainfo );
    if (!ai->isReady()) printf("Warning (AdaBootst::setAttributeInfo): AttributeInfo was not set appropriately\n");
    if (this->prv) { free(this->prv); this->prv = NULL; }
    this->prv = (double*)calloc( ai->getLabelCCount(), sizeof(double) );
    sx = ai->getAttrCount();
  }
  void setTrainingData(TrainingData<X_t> *tdp=NULL) {
    // Use provided training data. (It might be shared with others)
    this->tt = (tdata ? tdp : &tdata);
    if (!tt->isReady()) printf("Warning (AdaBootst::setTrainingData): TrainingData is not ready yet\n");
  }
  
  // -----------------------------------------------------------------
  // Training
  // -----------------------------------------------------------------
  
private:
  void initHypotheses(int nh) {
    if (!ai || !ai->isReady()) return;
    clear();
    this->nh = nh;
    //h    = new H_t[ nh ];  // does not work, don't know why.
    h  = (H_t*) calloc( nh, sizeof(H_t) );
    hw = (double*) calloc( nh, sizeof(double) );
    for (int i=0; i<nh; i++) {
      h[i].setAttributeInfo( ai );
      if (tt && tt->isReady()) h[i].setTrainingData( tt );
    }
  }
  
  // -----------------------------------------------------------------
  // Training for binary classification
  // -----------------------------------------------------------------
  
public:
  bool train(int nh, int wlopt, int dlist[]=NULL, int dsize=0, X_t D[]=NULL) {
    // Train binary classifier with given list of training data.
    //   nh      : number of iterations (hypotheses)
    //   wlopt   : weak learner option (max. depth for decision tree, example)
    //   dlist[] : list of indices to be used for training
    //   dsize   : number of data  to be used for training
    //   D[dsize]: initial weights for the training set [NULL]
    // "A decision-theoretic generalization of on-line learning and an application to boosting"
    //   by Yoav Freund and Robert E. Schapire, 1997
    if (!ai || !ai->isReady()) { if (verbose) fprintf(stderr, "[AdBst] Error: call setAttributeInfo() first.\n"); return false; }
    if (!tt || !tt->isReady()) { if (verbose) fprintf(stderr, "[AdBst] Error: call setTrainingData() first.\n"); return false; }
    if (ai->getLabelCCount()!=2) { if (verbose) fprintf(stderr, "[AdBst] Error: data is not binary\n"); return false; }
    int    t, i, *sdl=NULL; 
    // set pointer to the selected data list
    if (dlist) sdl = dlist;
    else {
      dsize = tt->nt;
      sdl   = (int*) malloc( tt->nt * sizeof(int) );
      for (i=0; i<tt->nt; i++) sdl[i] = i;
    }
    if (verbose) printf("[AdBst] training started with %d data  (WeakLearner option: %d)\n", dsize, wlopt);
    int  *y2 = (int*) calloc( dsize, sizeof(int) );
    initHypotheses( nh );
    // initialize weights for the training data
    double e, *w = (double*) calloc( dsize, sizeof(double) ), prob=0;
    if (D) for (i=0; i<dsize; i++) w[i] = D[i];
    else   for (i=0; i<dsize; i++) w[i] = 1.0 / dsize;
    MLH::PDF<double> wpdf(dsize, w, true);
    
    for (t = 0; t < nh; t++) {
      if (verbose) { printf("[AdBst] Expert %d \n", t); wpdf.printInfo("  Weights"); }
      if (!wpdf.normalize()) { this->nh = t;  break; }
      // run weak learner for the hypothesis at current iteration -- h[t]
      //h[t].verbose = this->verbose;
      h[t].train( wlopt, sdl, dsize, w );  // set max depth to 1 (at most two questions)
      h[t].verbose = false;
      // calculate the error of the hypothesis h[t]
      for (e = i = 0; i < dsize; i++) {
	y2[i] = h[t].test( tt->tx[sdl[i]], &prob );
	if (y2[i] != tt->ty[sdl[i]])  e += prob * w[i];
      }					// e must be in the range of [ 0, 0.5 ]
      if (e == 0) { hw[t] = 1.0; for (i=0; i<dsize; i++) w[i] = 0; }
      else {		// weak learner returns -1 or +1 (false or true)
	hw[t] = 0.5 * log( (1-e)/e );	// importance of the hypothesis [ 0, inf ]
	// adjust the weight vector (smaller when ty[i] == y2[i])
	for (i = 0; i < dsize; i++)
	  w[i] = w[i] * exp( -hw[t] * (tt->ty[sdl[i]] == y2[i] ? +1.0 : -1.0) );
      }
      if (verbose) {
	h[t].testOnTData("  TestExpertOnTData", sdl, dsize);
	char str[256];
	printf("  Ex%02d : e=%.4f alpha=%5.2f  %s \n", t, e, hw[t], h[t].getSingleLineInfo(str));
	if (hw[t] < 0) printf("  Error! - something wrong with the weak learner\n");
      }
    }
    for (hwsum = t = 0; t < nh; t++) hwsum += hw[t];
    free(w);   free(y2);
    last_nt = dsize;
    if (!dlist && sdl) free(sdl);
    return true;
  }
  
  // -----------------------------------------------------------------
  // Testing for binary classification
  // -----------------------------------------------------------------
private:  
  bool testPrivate(X_t x[]) {
    return true;
  }
public:  
  bool test(X_t x[], double *prob=NULL) {
    if (!x || nh<=0) { if (prob) *prob = 0.5; return false; }
    if (ai->getLabelCCount()!=2)  { if (verbose) fprintf(stderr, "[AdBst] Error: it's not a binary test\n"); return false; }
    double sum = 0;
    for (int i = 0; i < nh; i++) {
      sum += hw[i] * (h[i].test( x ) ? +1 : -1);
    }
    if (hwsum > 0) sum /= 2*hwsum;  // scale into [ -2.0 ~ +2.0 ]
    if (prob) *prob = exp(sum) / (exp(sum) + exp(-sum));
    return (sum >= 0 ? true : false);
  }
  int testOnTData(char *cmmt=NULL, int lst[]=NULL, int nlst=0, double w[]=NULL ) {
    // Test on a subset of the training data. (for debugging)
    if (!tt || !tt->isReady()) return 0;
    if (nlst <= 0 || !lst) { nlst = tt->nt; lst = NULL; }
    printf("%s(%d) : <", (cmmt ? cmmt : "TestOnDataset"), nlst);
    int  i, cnt=0, tidx;  double wsum=0;
    for (i = 0; i < nlst; i++) {
      tidx = ( lst ? lst[i] : i );
      bool succ = (test( tt->tx[tidx] ) == tt->ty[tidx]);
      if (nlst<=5000) { 
	printf("%c", (succ ? '.':('0'+tt->ty[tidx])));
	if (i%10==9) printf(" ");
      }
      if (succ) { cnt++;  if (w) wsum+=w[tidx]; }
    }
    if (nlst>5000) printf(" TooManyToShow ");
    if (w) printf(">  %d%% success (%.0f%% w)\n", cnt*100/nlst, wsum*100);
    else   printf(">  %d%% success\n", cnt*100/nlst);
    return cnt;
  }
  
  // -----------------------------------------------------------------
  // AdaBoost for a multi-classs classification  Y = {1 ... K}
  // -----------------------------------------------------------------
  
public:
  bool trainM1(int dsize, int dlist[], int nh, X_t D[]=NULL) {
    // Train multi-class classifier with given list of training data.
    //   dsize   : number of data  to be used for training
    //   dlist[] : list of indices to be used for training
    //   nh      : number of iterations (hypotheses)
    //   D[dsize]: initial weights for the training set [NULL]
    // "A decision-theoretic generalization of on-line learning and an application to boosting"
    //   by Yoav Freund and Robert E. Schapire, 1997
    if (!ai || !ai->isReady()) { if (verbose) fprintf(stderr, "[AdBst] Error: call setAttributeInfo() first.\n"); return false; }
    if (!tt || !tt->isReady()) { if (verbose) fprintf(stderr, "[AdBst] Error: call setTrainingData() first.\n"); return false; }
    if (ai->getLabelCCount()==2) { if (verbose) fprintf(stderr, "[AdBst] Error: data is binary. Not for 'trainM1()'\n"); return false; }
    if (verbose) printf("[AdBst] train() started\n");
    int   i, t;
    int   *y2 = (int*) calloc( dsize, sizeof(int) );
    initHypotheses( nh );
    // initialize weights for the training data
    double e, *w = (double*) calloc( dsize, sizeof(double) );
    if (D) for (i=0; i<dsize; i++) w[i] = D[i];
    else   for (i=0; i<dsize; i++) w[i] = 1.0 / dsize;
    MLH::PDF<double> wpdf( dsize, w );
    for (t = 0; t < nh; t++) {
      if (verbose) printf("  iter=%02d\n", t);
      if (!wpdf.normalize()) {
	if (verbose) printf("[AdBst] all weights are zero, train stopped at iter = %d\n", t);
	this->nh = t;  break;
      }
      // run weak learner for the hypothesis at current iteration -- h[t]
      h[t].verbose = this->verbose;
      h[t].train( 0, 0, dsize, dlist, w );   ////  INCOMPLETE
      // calculate the error of the hypothesis h[t]
      for (e = i = 0; i < dsize; i++) {
	y2[i] = h[t].test( tt->tx[dlist[i]] );
	if (y2[i] != tt->ty[dlist[i]])  e += w[i];
      }
      if (verbose) printf("  wlearner tested\n" );
      if (e > 0.5) { 
	if (verbose) printf("[AdBst] wlearner has too big error=%g at iteration %d\n", e, t );
	break;
      }
      // weight of the hypothesis h[t] (less than 1)
      hw[t] = e / ( 1 - e );
      if (verbose) {
	char str[256];
	printf("  Ex%d : e=%.4f beta=%5.2f  %s \n", t, e, hw[t], h[t].getSingleLineInfo(str));
      }
      // set the new weights vector
      for (i = 0; i < dsize; i++)
	if (y2[i] == tt->ty[dlist[i]]) w[i] = w[i] * hw[t];
    }
    this->nh = t;
    if (verbose) printf("[AdBst] train() ended with %d hypotheses\n", this->nh );
    free(w);   free(y2);
    last_nt = dsize;
    return true;
  }
  
  int  testM1(X_t *x) {
    // Evaluate a test data, and return the prediction of its label
    if (ai->getLabelCCount()<3) { if (verbose) fprintf(stderr, "[AdBst] Error: data is binary. Not for 'testM1()'\n"); return false; }
    double sum = 0, max_sum = -FLT_MAX;
    int  i, j, maxy=0;
    for (j = 0; j < ai->getLabelCCount(); j++) {
      for (i = 0; i < nh; i++) {
	if (h[i].test( x ) == j) sum += log( 1.0 / hw[i] );
	if (sum > max_sum) { max_sum = sum; maxy = j; }
      }
    }
    return maxy;
  }
  
  // -----------------------------------------------------------------
  // file I/O
  // -----------------------------------------------------------------
public:
  bool writeFile(const char *fname, const char *cmmt=NULL) {
    if (!ai || !ai->isReady()) { if (verbose) fprintf(stderr, "[AdBst] Error: call setAttributeInfo() first.\n"); return false; }
    FILE *fp = fopen( fname, "w+" );      // open the file
    if (fp == NULL) { 
      std::cerr << "Error (MLH::AdaBoost::writeFile): cannot open file '" << fname << "'" << std::endl;
      return false;
    }
    // write general information
    fprintf( fp, "# %s\n", (cmmt ? cmmt : "") );
    fprintf( fp, "# AdaBoost Trained Knowledge\n" );
    fprintf( fp, "# \n" );
    ai->writeAttributeInfo( fp );
    UTIL::Clock clock;
    if (last_nt > 0)
      fprintf(fp, "Training  : trained on %d data (saved at %s)\n", last_nt, clock.getLocalTime(4) );
    else
      fprintf(fp, "Training  : training data unknown (saved at %s)\n", clock.getLocalTime(4) );
    // write hypotheses
    fprintf( fp, "Experts   : %d\n", nh );
    for (int i=0; i < nh; i++) {
      fprintf( fp, "Ex%02d (%8.6f) : ", i, hw[i] );
      h[i].writeFileSingleLine( fp );
    }
    // close the file
    fclose( fp );
    return true;
  }

  bool readFile(const char *fname) {
    // open the file
    FILE *fp = fopen( fname, "r" );
    if (fp == NULL) { 
      std::cerr << "Error (AdaBoost::readFile): cannot open file '" << fname << "'" << std::endl;
      return false;
    }
    if (!ai || !ai->readAttributeInfo( fp )) { // read attribute information
      fclose( fp );
      std::cerr << "Error (AdaBoost::readFile): invalid AttributeInfo" << std::endl;
      return false;
    }
    setAttributeInfo();
    char line[1024]; 
    int  i, idx, n;
    bool ret = false;
    UTIL::WordBreaker wb;
    while (fgets(line, 250, fp)) {
      wb.parse( line, " :()\t\r\n", false, '#' );
      if (wb.count() < 1) continue;  // skip empty lines
      if        (strcmp(wb[0],"Training")==0) {	// -----------
      } else if (strcmp(wb[0],"Experts")==0) {	// -----------
	if (sscanf(wb[1], "%d", &n) != 1) break;
	initHypotheses( n );
	for (i=0; i < nh; i++) {
	  // read the weight of the hypothesis
	  if (fscanf( fp, "Ex%d (%lf) : ", &idx, hw+i ) != 2) break;
	  // read the hypothesis
	  if (! h[i].readFileSingleLine( fp )) break;
	}
	for (hwsum=i=0; i < nh; i++) hwsum += hw[i];
	ret = (i == nh);
      } else break;
    }
    fclose( fp );	// close the file
    return ret;
  }
  
  // -----------------------------------------------------------------
  // miscellaneous
  // -----------------------------------------------------------------
public:
  void printInfo(char *cmmt=NULL, bool single_line=false) {
    if (single_line) {
      if (nh==0) printf("%s: NOT READY\n", (cmmt ? cmmt : "AdaBoost"));
      else       printf("%s: %d experts (featureSize=%d with %d labels)\n", 
			(cmmt ? cmmt : "AdaBoost"), nh, sx, ai->getLabelCCount());
    } else {
      char blank[80], str[256];  int indent, i;
      for (indent=0; indent<70; indent++) if (!cmmt || cmmt[indent]!=' ') break;
      memset(blank, ' ', indent);  blank[indent]='\0';
      printf("%s (with %d hypotheses) \n", (cmmt ? cmmt : "AdaBoost"), nh);
      if (ai && ai->isReady()) ai->printInfo("  Att : ", true);
      for (i = 0; i < nh; i++)
	printf("%s  Ex%02d : alpha=%.4f  %s \n", blank, i, hw[i], h[i].getSingleLineInfo(str));
    }
  }
  
};


}	// end of namespace MLH

#endif // MLH_ADABOOST_HPP
