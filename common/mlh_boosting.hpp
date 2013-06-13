
//
// MLH::Boosting<X_t,H_t> class template for Binary Classification
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
//         - int   testOnTData( char *cmmt=NULL, int lst[]=NULL, int nlst=0 );
//         - void  writeFileSingleLine(FILE *fp);
//         - bool  readFileSingleLine(FILE *fp);
//         - char* getSingleLineInfo(char *buffer);
// Member functions
//   - void  setAttributeInfo( AttributeInfo<> *ap );
//   - void  setTrainingData ( TrainingData<>  *tp );
//   - void  train(int label, int nh, int wlopt, int dlist[]=NLL, int dsize=0);
//   - bool  test (X_t *x, double prob[2]=NULL);
//   - int   testOnTData(char *cmmt=NULL, int lst[]=NULL, int nlst=0);
//   - bool  writeFile(const char *fname, const char *cmmt=NULL);
//   - bool  readFile(const char *fname);
// 

#ifndef MLH_BOOSTING_HPP
#define MLH_BOOSTING_HPP

#include <iostream>
#include <fstream>
#include <cmath>
#include "util_word_breaker.hpp"
#include "util_clock.hpp"
#include "mlh_pdf.hpp"
#include "mlh_attribute_info.hpp"
#include "mlh_training_data.hpp"


namespace MLH {
  
template <class X_t, class H_t>
class Boosting {

  // prerequisite information ----------------------------------------
public:
  AttributeInfo<X_t>	*ai;		// attribute information
  TrainingData<X_t>	*tt;		// training data
private:
  AttributeInfo<X_t>	ainfo;		// attribute information (local)
  TrainingData<X_t>	tdata;		// training data dummy
  // experts (weak learners) -----------------------------------------
public:
  int			target_label;
  int			nh;		// number of hypotheses
  H_t			*mh;		// hypotheses
  double		*hw;		// weights for hypotheses
  double		hwsum;
  // test result -----------------------------------------------------
public:
  double		prv[2];
  // etc -------------------------------------------------------------
public:
  bool  		verbose;
  int			last_nt;
  
public:
  Boosting() : ai(&ainfo), tt(&tdata), nh(0), mh(NULL), hw(NULL), verbose(false), last_nt(0) {}
  ~Boosting() { clear(); }
  void clear(void) {
    for (int i=0; i<nh && mh; i++) mh[i].clear();
    // if (mh) delete(mh);  mh = NULL;   nh = 0; ////
    if (mh) free(mh);    mh = NULL;   nh = 0;
    if (hw) free(hw);    hw = NULL;
  }
  inline bool  isReady(void) { return ai->isReady() && nh > 0; }
  bool isValid(void) {
    for (int i=0; i<nh; i++) if (!mh[i].isValid()) return false;
    return true;
  }
  
  // -----------------------------------------------------------------
  // Setting Data Information
  // -----------------------------------------------------------------
public:	// use the data information given by the user
  void setAttributeInfo(AttributeInfo<X_t> *ainfo=NULL) {
    // Use provided attribute information. (It might be shared with others)
    this->ai = ( ainfo ? ainfo : &this->ainfo );
    if (!ai->isReady()) printf("Warning (Bootsting::setAttributeInfo): AttributeInfo is not ready yet\n");
  }
  void setTrainingData(TrainingData<X_t> *tp=NULL) {
    // Use provided training data. (It might be shared with others)
    this->tt = ( tp ? tp : &this->tdata );
    if (!tt->isReady()) printf("Warning (Bootsting::setTrainingData): TrainingData is not ready yet\n");
  }
  
  // -----------------------------------------------------------------
  // Training
  // -----------------------------------------------------------------
  
public:
  void initHypotheses(int nh) {
    if (!ai || !ai->isReady()) return;
    clear();
    this->nh = nh;
    // mh = new H_t[nh];			//// 
    mh = (H_t*)calloc(nh, sizeof(H_t));
    hw = (double*) calloc( nh, sizeof(double) );
    for (int i=0; i<nh; i++) {
      mh[i].setAttributeInfo( ai );
      if (tt && tt->isReady()) mh[i].setTrainingData( tt );
      hw[i] = 1.0;
    }
  }
  
public:
  bool train(int label, int nh, int wlopt, int dlist[]=NULL, int dsize=0) {
    // Train binary classifier with given list of training data.
    //   label   : the label to learn about (because this classifier is a binary learner)
    //   nh      : number of iterations (hypotheses)
    //   wlopt   : weak learner option (max. depth for decision tree, example)
    //   dlist[] : list of indices to be used for training
    //   dsize   : number of data  to be used for training
    // "A decision-theoretic generalization of on-line learning and an application to boosting"
    //   by Yoav Freund and Robert E. Schapire, 1997
    if (!ai || !ai->isReady()) { if (verbose) fprintf(stderr, "[Boost] Error: call setAttributeInfo() first.\n"); return false; }
    if (!tt || !tt->isReady()) { if (verbose) fprintf(stderr, "[Boost] Error: call setTrainingData() first.\n"); return false; }
    int    t, i, *sdl=NULL, yy=0, ty, positive=0;
    this->target_label = label;
    // set pointer to the selected data list
    if (dlist) sdl = dlist;
    else {
      dsize = tt->nt;
      sdl   = (int*) malloc( tt->nt * sizeof(int) );
      for (i=0; i<tt->nt; i++) sdl[i] = i;
    }
    for (i=0; i<dsize; i++) if (tt->ty[ sdl[i] ]==target_label) positive++;
    if (verbose) printf("[Boost] training on %d data (class%d: %d) for %d experts with (option=%d)\n", 
			dsize, target_label, positive, nh, wlopt);
    tt->initTrainingWeights();				// initialize training data weights
    for (i=0; i<dsize; i++) tt->w[i] = 1.0 / dsize;
    double *cf = (double*) calloc( dsize, sizeof(double) );  // data confidence
    initHypotheses( nh );				// initialize weak learners
    MLH::PDF<double> wpdf(dsize, tt->w, true);
    double avgcf=0;
    for (t = 0; t < nh; t++) {
      if (verbose) { printf("[Boost] Expert %d / %d \n", t, nh); wpdf.printInfo("  Weights"); }
      if (!wpdf.normalize()) { this->nh = t;  break; }
      // run weak learner for the hypothesis at current iteration -- mh[t]
      mh[t].train( label, wlopt, NULL, 0 );
      mh[t].verbose = false;
      if (verbose) mh[t].testOnTData("  TestExpertOnTData", sdl, dsize );
      double confidence=0, cfsum=0, ecount=0;
      for (i = 0; i < dsize; i++) {
	bool ret = mh[t].test( tt->tx[sdl[i]], NULL, &confidence );
	yy = (ret==false ? -1 : +1);
	cf[i] += yy * confidence;			// in range of [ -10 , +10 ]
	ty = (tt->ty[sdl[i]] != label ? -1 : +1);
	tt->w[i] = 1 / ( 1 + exp( +ty * cf[i] ) );
	cfsum   += 1 / ( 1 + exp( -ty * cf[i] ) );	// for average confidence [ 0, inf ]
	ecount += (ty * cf[i] < 0 ? 1 : 0);		// for training error
      }
      if (verbose) {
	char str[256];
	printf("  MeanConfidence = %.4f   TrainingError = %.4f\n", cfsum/dsize, ecount/dsize);
	printf("  Ex%02d : a=%.4f  %s \n", t, hw[t], mh[t].getSingleLineInfo(str));
      }
      if (fabs(cfsum/dsize - avgcf) < 1e-4) { 
	if (verbose) printf("[Boost] training stopped due to too small confidence increase\n"); 
	this->nh = t+1;  break; 
      }
      avgcf = cfsum / dsize;
    }
    for (hwsum = t = 0; t < nh; t++) hwsum += hw[t];
    free(cf);
    if (!dlist && sdl) free(sdl);
    last_nt = tt->nt;
    return true;
  }
  
  // -----------------------------------------------------------------
  // Testing
  // -----------------------------------------------------------------
public:  
  bool test(X_t x[], double prob[2]=NULL) {
    if (!x)         { fprintf(stderr, "[Boost] Warning: it's NOT trained yet\n"); return false; }
    if (!isReady()) { fprintf(stderr, "[Boost] Error: it's not ready for test\n"); return false; }
    //if (debug) { printf("[Boost] Confidences from %d experts: ", nh); fflush(stdout); }
    double confidence=0, confsum = 0;
    for (int i = 0; i < nh; i++) {
      bool ret = mh[i].test( x, NULL, &confidence );
      confsum += (ret==false ? -1 : +1) * confidence;
      //printf("testResult %d  confidence=%g  consum=%g\n", ret, confidence, confsum);
    }
    //if (debug) { printf("\n"); }
    // probability by Sigmoid conversion from the log-ratio :
    //   Pr [ y==true | x ] = 1 / (1 + e^{-f(x)})
    //   where f(x) is the weighted sum of the log-likelihood ratio of all the experts
    double pr = 1 / (1 + exp(-confsum));
    prv[0] = (1 - pr);
    prv[1] = (pr);
    if (prob) { prob[0] = prv[0]; prob[1] = prv[1]; }
    //printf("prv[2]={%.4f %.2f}\n", prv[0], prv[1]);
    return  (prv[0] > prv[1] ? false : true);
  }
  int testOnTData(char *cmmt=NULL, int lst[]=NULL, int nlst=0) {
    // Test on a subset of the training data. (for debugging)
    if (!tt || !tt->isReady()) return 0;
    if (nlst <= 0 || !lst) { nlst = tt->nt; lst = NULL; }
    printf("%s(%d) : <", (cmmt ? cmmt : "TestOnDataset"), nlst);
    int  i, cnt=0, tidx;
    double wsum=0, prob[2]={0,0}, psum=0;
    for (i = 0; i < nlst; i++) {
      tidx = ( lst ? lst[i] : i );
      bool ret  = test( tt->tx[tidx], prob );
      bool succ = ( ret == true  && tt->ty[tidx] == target_label ||
		    ret == false && tt->ty[tidx] != target_label );
      if (nlst<=5000) { 
	printf("%c", (succ ? '.':('0'+tt->ty[tidx])));
	if (i%10==9) printf(" ");
      }
      psum += prob[ (ret ? 1:0) ];
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
    if (!ai || !ai->isReady()) { if (verbose) fprintf(stderr, "[Boost] Error: call setAttributeInfo() first.\n"); return false; }
    FILE *fp = fopen( fname, "w+" );      // open the file
    if (fp == NULL) { 
      std::cerr << "Error (MLH::Boosting::writeFile): cannot open file '" << fname << "'" << std::endl;
      return false;
    }
    bool ret = writeFile( fp, cmmt );
    fclose( fp );
    return ret;
  }
  bool writeFile(FILE *fp, const char *cmmt=NULL) {
    if (!fp) return false;
    if (!ai || !ai->isReady()) { if (verbose) fprintf(stderr, "[Boost] Error: call setAttributeInfo() first.\n"); return false; }
    // write general information
    fprintf( fp, "# %s\n", (cmmt ? cmmt : "") );
    fprintf( fp, "# Boosting Trained Knowledge\n" );
    fprintf( fp, "# \n" );
    ai->writeAttributeInfo( fp );
    UTIL::Clock clock;
    if (last_nt > 0) 
      fprintf(fp, "Training  : trained on %d data (saved at %s)\n", last_nt, clock.getLocalTime(4) );
    else
      fprintf(fp, "Training  : training data unknown (saved at %s)\n", clock.getLocalTime(4) );
    // write hypotheses
    fprintf( fp, "Class  %d  : %d Experts  (%s==%s)\n", target_label, nh, ai->getLabelName(), ai->convI2S(-1,target_label) );
    for (int i=0; i < nh; i++) {
      fprintf( fp, "Ex%02d (%8.6f) : ", i, hw[i] );
      mh[i].writeFileSingleLine( fp );
    }
    return true;
  }

  bool readFile(const char *fname) {
    // open the file
    FILE *fp = fopen( fname, "r" );
    if (fp == NULL) { std::cerr << "Error (Boosting::readFile): cannot open file '" << fname << "'" << std::endl; return false; }
    bool ret = readFile( fp );
    fclose(fp);
    return ret;
  }
  bool readFile(FILE *fp) {
    if (!fp) return false;
    if (!ai || !ai->readAttributeInfo( fp )) { // read attribute information
      std::cerr << "Error (Boosting::readFile): invalid AttributeInfo" << std::endl;
      return false;
    }
    setAttributeInfo();
    char line[1024]; 
    int  i, idx, c, n;
    bool ret = false;
    UTIL::WordBreaker wb;
    while (fgets(line, 250, fp)) {
      wb.parse( line, " :()\t\r\n", false, '#' );
      if (wb.count() < 1) continue;  // skip empty lines
      if        (strcmp(wb[0],"Training")==0) {	// -----------
      } else if (strcmp(wb[0],"Class")==0) {	// -----------
	c = atoi(wb[1]);
	n = atoi(wb[2]);
	if (c < 0 || c >= ai->getLabelCCount() || n<1) break;
	target_label = c;
	initHypotheses( n );
	for (i=0; i < nh; i++) {
	  // read the weight of the hypothesis
	  if (fscanf( fp, "Ex%d (%lf) : ", &idx, hw+i ) != 2) break;
	  // read the hypothesis
	  mh[i].target_label = target_label;
	  if (! mh[i].readFileSingleLine( fp )) break;
	}
	for (hwsum=i=0; i < nh; i++) hwsum += hw[i];
	ret = (i == nh);
	target_label = c;
	break;
      } else break;
    }
    return ret;
  }
  
  // -----------------------------------------------------------------
  // miscellaneous
  // -----------------------------------------------------------------
public:
  void printExpert(int i, char* cmmt=NULL) { if (i>=0 && i < nh) mh[i].printInfo(cmmt); }
  
  void printInfo(char *cmmt=NULL, bool single_line=false) {
    if (single_line) {
      if (nh==0) printf("%s: NOT READY\n", (cmmt ? cmmt : "Boosting"));
      else       printf("%s: %2d attr, %d label => %2d experts \n", 
			(cmmt ? cmmt : "Boosting"), ai->getAttrCount(), ai->getLabelCCount(), nh);
    } else {
      char blank[80], str[256];  int indent, i;
      for (indent=0; indent<70; indent++) if (!cmmt || cmmt[indent]!=' ') break;
      memset(blank, ' ', indent);  blank[indent]='\0';
      printf("%s (with %d experts) \n", (cmmt ? cmmt : "Boosting"), nh);
      if (ai && ai->isReady()) ai->printInfo("  Att : ", true);
      if (tt && tt->isReady()) tt->printInfo("  TrD : ");
      printf("%s  Class %d : %d experts  (%s==%s)\n", blank, target_label, nh, ai->getLabelName(), ai->convI2S(-1,target_label));
      for (i = 0; i < nh; i++)
	printf("%s  Ex%02d : a=%.4f  %s \n", blank, i, hw[i], mh[i].getSingleLineInfo(str));
    }
  }
  
  
};


}	// end of namespace MLH

#endif // MLH_BOOSTING_HPP
