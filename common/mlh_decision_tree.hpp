
//
// MLH::DecisionTree<X_t> : Decision Tree for Multi-Class Classification
//
// Jaeil Choi
// Last modified in July, 2009
//
// A simple decision tree for multi-class classification problem,
//   (This class was developed as a weak learner for Boosting algorithm.)
//   template data types:
//     X_t        : type of training input data (float or double)
// Member functions
//   - void  clear(void);
//   - void  setAttributeInfo( AttributeInfo<> *ap );
//   - void  setTrainingData ( TrainingData<>  *tp );
//   - bool  train( int maxdepth, int dlist[]=NLL, int dsize=0 );
//   - int   test ( X_t x[], double *prob=NULL, double *confidence=NULL );
//   - int   testOnTData( char *cmmt=NULL, int lst[]=NULL, int nlst=0 );
//   - void  writeFileSingleLine(FILE *fp);
//   - bool  readFileSingleLine(FILE *fp);
//   - char* getSingleLineInfo(char *buffer);
// Member function arguments:
//     ap         : pointer to MLH::AttributeInfo<X_t>, for the attribute information
//     tp         : pointer to MLH::TrainingData<X_t>, for the training data
//     maxdepth   : train(); maximum depth limit of the tree [optional]
//     dlist[]    : train(); list of indices (of tx[] and ty[]) to be used for training [optional]
//     dsize      : train(); size of input index list [optional]
//     w[dsize]   : train(); weights for each of given feature vectors [optional]
//     x[sx]      : test(); a feature vector to be evaluated for test()
//     prob       : test(); calculated probability for the returned label [optional]
//     confidence : test(); log likelihood ratio for the returned label; log(Positive/Negative)/2
//


#ifndef MLH_DECISION_TREE_HPP
#define MLH_DECISION_TREE_HPP

#include <iostream>
#include <cmath>
#include <cstdarg>
#include "mtx_matrix.hpp"
#include "mlh_attribute_info.hpp"
#include "mlh_training_data.hpp"

#define	Log_2		0.69314718055994530942
#define	Log2(x)		((x) <= 0 ? 0.0 : log((double)x) / Log_2)
#define DT2V_COPY(r,a)      do { (r)[0]=(a)[0]; (r)[1]=(a)[1]; } while(0)
#define DT2V_SET(r,x,y)     do { (r)[0]=(x); (r)[1]=(y); } while(0)
#define DT3V_SET(r,x,y,z)   do { (r)[0]=(x); (r)[1]=(y); (r)[2]=(z); } while(0)
#define DT4V_SET(r,x,y,z,w) do { (r)[0]=(x); (r)[1]=(y); (r)[2]=(z); (r)[3]=(w); } while(0)


// ===================================================================
// Training Data
// ===================================================================
  
namespace MLH {
  
template <class X_t>
class DecTreeTraining {
public:
  // training data
  int			sx;
  int			maxdepth;		// depth limit
  // temporary buffers for tree building
  double		*lbldist;	// temporary buffer for label distribution of a leaf node
  double		*tmpInfo;
  double		*tmpGain;
  X_t			*tmpAVal;
public:
  DecTreeTraining() : sx(0), lbldist(NULL), tmpInfo(NULL), tmpGain(NULL), tmpAVal(NULL) {}
  ~DecTreeTraining() { clearTrainingBuffers(); }
  bool isReady(void) { return (sx > 0 && lbldist != NULL); }
  void setupTrainingBuffers(int sx, int nl, int max_depth) {
    // Allocate temporary memory buffers for training
    this->sx = sx;
    lbldist = (double*) malloc( (nl+1) * sizeof(double) );
    tmpInfo = (double*) malloc( sx * sizeof(double) );
    tmpGain = (double*) malloc( sx * sizeof(double) );
    tmpAVal = (X_t*) malloc( sx * sizeof(X_t) );
    maxdepth = (max_depth>0 ? max_depth : 65535);
  }
  void clearTrainingBuffers(void) {
    // Clear temporary memory buffers for training
    if (lbldist) free(lbldist);    lbldist = NULL;
    if (tmpInfo) free(tmpInfo);    tmpInfo = NULL;
    if (tmpGain) free(tmpGain);    tmpGain = NULL;
    if (tmpAVal) free(tmpAVal);    tmpAVal = NULL;
  }
  void printAttrGains(char *cmmt=NULL, bool gain=true, bool info=true) {
    // Print the information gain for each attribute (DURING the training)
    if (!tmpGain || !tmpInfo) return;
    if (gain) {
      printf("%s  AttributeGain[]: ", (cmmt ? cmmt:""));
      for (int a=0; a<sx; a++) printf("%.4f ", tmpGain[a]);
      printf("\n");
    }
    if (info) {
      printf("%s  AttributeInfo[]: ", (cmmt ? cmmt:""));
      for (int a=0; a<sx; a++) printf("%.4f ", tmpGain[a]);
      printf("\n");
    }
  }
};
  
}

// ===================================================================
// Decision Tree Node
// ===================================================================
  
namespace MLH {
  
class PruningData { 
public:
  double cerr, perr; 
public:
  PruningData() {}
  PruningData(double ce, double pe) : cerr(ce), perr(pe) {}
  ~PruningData() {}
};
  
template <class X_t>
class DTNode {
public:
  char			type;		// 'L'eaf, 'D'iscrete, or 'C'ontinous
  bool			visited;	// visited by 'test()' function (OnDecisionPath?)
  union {
    struct {	// Leaf Node
      int		label;		// most frequent class at this node
      double		items;		// num of items  at this node
      double		errors;		// num of errors at this node
      double		*ldist;		// label distribution [nLabels]
    } L;
    struct {	// Decision Node
      int		aidx;		// attribute index
      X_t		value;		// attribute value
      DTNode		*child[2];	// child nodes ( '==' and '!=' or '<' and '>=' )
      float		gain;
      int		flabel;		// most frequent label at this node
    } N;
  };
  PruningData		*pdata;
  // 
public:
  DTNode() { memset(this, 0, sizeof(DTNode<X_t>)); }
  ~DTNode() { clear(); }
  void clear(void) {
    if      (type=='L') { 
      if (L.ldist) free(L.ldist); L.ldist = NULL;
    } else if (type=='D' || type=='C') {
      if (N.child[0]) delete(N.child[0]); N.child[0] = NULL; 
      if (N.child[1]) delete(N.child[1]); N.child[1] = NULL;
    }
    if (pdata) { free(pdata); pdata = NULL; }
  }
  void setLeafNode(int label, double wsum, double wdist[], int nLabels) {
    clear();
    this->type = 'L';
    this->L.label = label;
    this->L.items = wsum;
    if (wdist) {
      this->L.ldist = (double*)calloc( nLabels, sizeof(double) );
      memcpy( this->L.ldist, wdist, nLabels*sizeof(double) );
    } else {
      this->L.ldist = NULL;
    }
  }
  void setDecisionNode(bool discrete, int aidx, X_t value) {
    clear();
    this->type = (discrete ? 'D' : 'C');
    this->N.aidx  = aidx;
    this->N.value = value;
  }
};
  
}

// ===================================================================
// Decision Tree
// ===================================================================
  
namespace MLH {
  
template <class X_t>
class DecisionTree {
  // prerequisite information ----------------------------------------
public:
  AttributeInfo<X_t>	*ai;		// attribute information
  TrainingData<X_t>	*tt;		// training data
private:
  AttributeInfo<X_t>	ainfo;		// attribute information (local)
  TrainingData<X_t>	tdata;		// training data dummy
  // the decision tree -----------------------------------------------
public:
  DTNode<X_t>		*root;		// the root node
  int			nnodes;		// number of decision nodes
  int			nleaf;		// number of leaf node
  int			mdepth;		// maximum depth
  double		totalw;		// total weight
  bool			pruning;
  bool			verbose;
  // etc -------------------------------------------------------------
private:
  DecTreeTraining<X_t>	*ti;
  double		*prv;	// temporary buffer for tree testing
  
public:
  DecisionTree() : ai(&ainfo), tt(&tdata), root(NULL), nnodes(0), nleaf(0), mdepth(0), 
		   pruning(true), verbose(false), ti(NULL), prv(NULL) {}
  ~DecisionTree() { clear(); }
  void clear(void)  { delete(root); root = NULL;  nnodes = nleaf = mdepth = 0; }
  bool isValid(void) { return (root && (root->type=='L' || root->type=='D' || root->type=='C')); }
  int  getTDataCount(void) { return tt->nt; }
  inline char* getLabelName(int lbl) { return ai->convI2S( -1, lbl ); }
  
  // -----------------------------------------------------------------
  // Setting Data Information
  // -----------------------------------------------------------------
public:	// use the data information given by the user
  void setAttributeInfo(AttributeInfo<X_t> *ainfo=NULL) {
    this->ai = ( ainfo ? ainfo : &this->ainfo );
    if (!ai->isReady()) printf("Warning (DecisionTree::setAttributeInfo): AttributeInfo is not ready yet\n");
    if (prv) { free(prv); prv = NULL; }
    prv = (double*)calloc( ai->getLabelCCount(), sizeof(double) );
  }
  void setTrainingData(TrainingData<X_t> *tdp=NULL) {
    this->tt = (tdp ? tdp : &tdata);
    if (!tt->isReady()) printf("Warning (DecisionTreesetTrainingData): TrainingData is not ready yet\n");
  }
  
  // -----------------------------------------------------------------
  // Training
  // -----------------------------------------------------------------
public:
  bool train(int maxdepth, int dlist[]=NULL, int dsize=0, double w[]=NULL) {
    // Create a decision tree using given training datasets, labels, and weights.
    // This function assumes that the training data was already set.
    //   maxdepth    : maximum depth (0: not restricted)
    //   dlist[dsize]: list of indices to be used for training
    //   dsize       : number of data  to be used for training
    //   w[total]    : list of weights for the entire feature vectors (total>=dsize)
    // Note that this function assumes 'setTrainingData()' was already called.
    if (!ai || !ai->isReady()) { printf("Error (DecisionTree::train): DataInfo not set\n"); return false; }
    if (!tt || !tt->isReady()) { printf("Error (DecisionTree::train): Training data not set\n"); return false; }
    if (dsize<0||dsize>tt->nt) { printf("Error (DecisionTree::train): invalid argument\n"); return false; }
    clear();	// clear any existing decision nodes
    ti = new DecTreeTraining<X_t>;
    ti->setupTrainingBuffers( tt->sx, tt->nl, maxdepth );
    int  i, szinfo[4], *tset=NULL, tsize=0, *vset=NULL, vsize=0;
    if (dlist) {             // use everything in the list, without pruning
      tsize = dsize;
      tset  = dlist;
    } else if (!pruning) {   // use everything in the training set, without pruning
      tsize = tt->nt;
      tset  = (int*) malloc( tsize * sizeof(int) );  // training set
      vset  = NULL;  vsize = 0;                      // validation set
      for (i=0; i<tt->nt; i++) tset[i] = i;
    } else {                 // use just half of the training set, with post-pruning
      tsize = vsize = tt->nt / 2;
      tset  = (int*) malloc( (tsize+1) * sizeof(int) );  // training set
      vset  = (int*) malloc( (vsize+1) * sizeof(int) );  // validation set
      for (i=0; i<tt->nt; i++) 
	if (i%2==0) tset[i>>1] = i; else vset[i>>1] = i;
    }
    if (verbose) printf("[DTree] training started on %d data, %s\n", tsize, (w ? "with given weight":"without weights"));
    DT4V_SET( szinfo, 0,0,0,0 );	// szinfo[3]: depth of current node (input)
    // create the decision tree
    root = createDecisionNode( tset, tsize, szinfo );
    totalw = sumupTree( root, szinfo, NULL );
    nnodes = szinfo[0];			// number of decision nodes in the tree (output)
    nleaf  = szinfo[1];			// number of leaf nodes in the tree (output)
    mdepth = szinfo[2];			// maximum depth of the tree (output)
    // prune the decision tree
    if (vset && vsize>0) {
      if (verbose) printf("[DTree] pruning the tree of size (%d+%d %d) using %d test data\n", nnodes, nleaf, mdepth, vsize);
      prune( vset, vsize, tset, tsize, szinfo );
      //printf("[DTree] pruned  the tree of size (%d+%d %d) using %d test data\n", nnodes, nleaf, mdepth, vsize);
    }
    if (!dlist) {
      if (tset) { free(tset); tset = NULL; tsize = 0; }
      if (vset) { free(vset); vset = NULL; vsize = 0; }
    }
    if (ti) { delete(ti); ti = NULL; }
    return true;
  }
  
private:
  DTNode<X_t>* createDecisionNode(int *lst, int nlst, int szinfo[4]) {
    //    lst     : index list of the training data
    //    nlst    : size of the index list
    //    szinfo[0] : number of decision nodes in this subtree (botton-up)
    //    szinfo[1] : number of leaf nodes in the subtree (bottom-up)
    //    szinfo[2] : maximum depth of this subtree (bottom-up)
    //    szinfo[3] : depth of current node (top-down)
    if (lst == NULL || nlst < 1) return NULL;
    char blank[80];  makeBlank(blank, szinfo[3]);
    // create a new decision node
    DTNode<X_t> *np = new DTNode<X_t>;
    // setup the best decision criteria (and its subsets)
    int    i, nLabels = ai->getLabelCCount(), aidx=0;
    int    lidx = findMostFrequentLabel( lst, nlst );
    double wsum = ti->lbldist[nLabels], wchild[2];
    X_t    criteria=0;
    bool   leaf = (ti->lbldist[lidx] == wsum || nlst<=1 || szinfo[3] >= ti->maxdepth);
    bool   unable = (!leaf && (aidx = findBestAttribute( lst, nlst, wsum, &criteria, wchild )) < 0);
    if (wsum <= 0) {			// ERROR =======================
      delete(np); return NULL;
    } else if (leaf || unable) { 	// LEAF ========================
      if (ti->lbldist[lidx]==wsum) np->setLeafNode( lidx, wsum, NULL,        nLabels );
      else                         np->setLeafNode( lidx, wsum, ti->lbldist, nLabels );
      if (verbose) {
	printf("%slabel = %d  with items=%.2f [ ", blank, np->L.label, np->L.items);
	for (i=0; i<nLabels && unable; i++) printf("%.2f ", np->L.ldist[i]);
	printf("]\n");
	if (unable) {
	  printf("Error (DecisionTree::createDecisionNode): unable to find the best attribute\n");
	  //tt->printTrainingData( "Unable", lst, nlst );
	}
      }
      DT3V_SET( szinfo, 0, 1, 0 );
    } else {				// Decision node =============
      int *sublst[2]={NULL,NULL}, lpos, rpos, nsize[2][4];
      sublst[0] = (int*) calloc( nlst, sizeof(int) );
      sublst[1] = (int*) calloc( nlst, sizeof(int) );
      DT4V_SET( nsize[0], 0, 0, 0, szinfo[3]+1 );
      DT4V_SET( nsize[1], 0, 0, 0, szinfo[3]+1 );
      bool discrete = ai->getAttr(aidx)->discrete;
      if (discrete) {	// DISCRETE ====================
	np->setDecisionNode( true, aidx, criteria );
	// separate the list into two sets
	for (i=lpos=rpos=0; i < nlst; i++) {
	  if ( tt->tx[ lst[i] ][ np->N.aidx ] == np->N.value ) sublst[0][ lpos++ ] = lst[i];
	  else                                                 sublst[1][ rpos++ ] = lst[i];
	}
      } else {		// CONTINUOUS ==================
	np->setDecisionNode( false, aidx, criteria );
	// separate the list into two sets
	for (i=lpos=rpos=0; i < nlst; i++) {
	  if ( tt->tx[ lst[i] ][ np->N.aidx ] < np->N.value ) sublst[0][ lpos++ ] = lst[i];
	  else                                                sublst[1][ rpos++ ] = lst[i];
	}
      }
      if (verbose) {
	printf("%snode: ( att[%d] %c %g ) => cnt(%d %d)  ", blank, np->N.aidx, (discrete ? '=':'<'), np->N.value, lpos, rpos);
	ti->printAttrGains("", true, false);
      }
      np->N.flabel   = lidx;			// save the most frequent label of this subtree
      np->N.gain     = ti->tmpGain[aidx];	// save the gain of this decision node
      // process its left/right subsets
      np->N.child[0] = createDecisionNode( sublst[0], lpos, nsize[0] );
      np->N.child[1] = createDecisionNode( sublst[1], rpos, nsize[1] );
      // complete the size and depth of the current subtree
      szinfo[0] =  nsize[0][0] + nsize[1][0] + 1;	// number of decision nodes in the subtree
      szinfo[1] =  nsize[0][1] + nsize[1][1];	// number of leaf nodes in the subtree
      szinfo[2] = (nsize[0][2] > nsize[1][2] ? nsize[0][2]+1 : nsize[1][2]+1);  // max depth
      if (sublst[0]) free(sublst[0]);
      if (sublst[1]) free(sublst[1]);
      ai->getAttr( aidx )->importance += wsum;
    }
    return np;
  }
  
  void makeBlank(char buf[], int i, int s=2) { memset(buf,' ',(i+1)*s);  buf[(i+1)*s]='\0'; }
  
  // -----------------------------------------------------------------
  // finding the optimal decision (question & criteria)
  // -----------------------------------------------------------------
private:
  int findMostFrequentLabel(int *lst, int nlst) {
    int  i, l, cidx, lbest;
    double *lbldist = ti->lbldist;
    int  *ty = tt->ty, nLabels = ai->getLabelCCount();
    for (l = 0; l <= nLabels; l++) lbldist[l] = 0;
    for (i = 0; i < nlst; i++) {
      cidx = lst[i];
      if (cidx < 0 || cidx > tt->nt) printf("error 1\n");
      double wv = ( tt->w ? tt->w[ cidx ] : 1.0 );
      if (!ty || ty[cidx] < 0 || ty[cidx] >= nLabels) printf("error 2\n");
      lbldist[ ty[ cidx ] ] += wv;
      lbldist[ nLabels ]    += wv;
    }
    for (l = 1, lbest=0; l < nLabels; l++)
      if (lbldist[l] > lbldist[lbest]) lbest = l;
    return lbest;
  }
  
  int findBestAttribute(int *lst, int nlst, double wsum, X_t *criteria, double wchild[2]) {
    int i, tidx, aidx, nLabels = ai->getLabelCCount(), sx = ai->getAttrCount();
    bool debug = false; // (nlst==3 && lst && lst[0]==8);
    DT2V_SET( wchild, 0, 0 );
    if (debug) printf("[DTree] findBestAttribute for %d data among %d attributes: \n", nlst, sx);
    int  *tlst = (int*) malloc( nlst * sizeof(int) );	// 'lst[]' should not be shuffled
    memcpy( tlst, lst, nlst * sizeof(int) );		// because of the sync with 'w[]'
    for (aidx = 0; aidx < sx; aidx++) {
      Attribute *att = ai->getAttr( aidx );
      int nvalues = att->d.ncases;
      double wusum = 0, wvsum = 0;
      if (att->discrete) {
	MTX::Matrix<double> mFreq (nvalues, nLabels, true), mFreqC, mPosi, mNega;
	// calculate frequency of each label for each value of this attribute
	for (i=0; i<nlst; i++) {			// mFreq is (ncases x nLabels)
	  tidx = tlst[i];
	  X_t value = tt->tx[tidx][aidx];
	  if (value==UNKNOWN) { wusum += (tt->w ? tt->w[tidx]:1); continue; }
	  else mFreq((int)value,(int)tt->ty[tidx]) += (tt->w ? tt->w[tidx]:1); 
	}
	mFreqC.sumOverRow( mFreq );  // ( 1 x nLabels ) frequency for each class label
	wvsum = wsum - wusum;		
	// find the best attribute value to split with
	double thisInfo, min_info=0;
	int    vidx, min_vidx=-1;
	for (vidx = 0; vidx < nvalues; vidx++) {
	  mPosi.set( 1, nLabels, &mFreq(vidx,0) );	// (1 x nLabels) for positive case
	  mNega.sub( mFreqC, mPosi );			// (1 x nLabels) for negative cases
	  thisInfo = ( calcInfo( mPosi.data, nLabels ) + calcInfo( mNega.data, nLabels ) );
	  if (vidx==0 || thisInfo < min_info) {
	    min_vidx = vidx;  min_info = thisInfo; 
	    DT2V_SET( wchild, mPosi.getRowSum(0), mNega.getRowSum(0) );
	  }
	}
	// remember the gain of the attribute
	ti->tmpGain[aidx] = calcInfo( mFreqC.data, nLabels ) / wvsum  -  min_info / wvsum;
	ti->tmpInfo[aidx] = calcInfo( wchild, 2 ) / wvsum;
	ti->tmpAVal[aidx] = (X_t)min_vidx;
	if (debug) printf("  att[%d] : %.2f/%.2f \n", aidx, ti->tmpGain[aidx], ti->tmpInfo[aidx]);
      } else {
	MTX::Matrix<double> mFreq(2, nLabels, true), mFreqC(1,nLabels, true);
	for (i=0; i<nlst; i++) {
	  tidx = tlst[i];
	  if (tt->tx[tidx][aidx] == UNKNOWN) { wusum += (tt->w ? tt->w[tidx]:1); continue; }
	  mFreqC((int)tt->ty[tidx]) += (tt->w ? tt->w[tidx]:1); 
	}
	wvsum = (wsum - wusum);
	if (debug) printf("  att[%02d] : ", aidx);
	// sort the list 'tlst[]' wrt the attribute values (ascending order)
	sortDataIndices( tlst, aidx, 0, nlst-1 );
	// find the best attribute cut
	mFreq.copySubMatrixFrom( 1, 0, mFreqC );
	double thisInfo, min_info=0;
	int    vidx, min_vidx=-1, trial=0;
	//if (debug) mFreq.printMatrix("%2.0f");
	for (vidx = 0; vidx < nlst-1; vidx++) {
	  tidx = tlst[vidx];
	  if (tt->tx[tidx][aidx] == UNKNOWN) continue;
	  mFreq(0, (int)tt->ty[tidx]) += (tt->w ? tt->w[tidx]:1);  // mFreq is (2 x nLabels)
	  mFreq(1, (int)tt->ty[tidx]) -= (tt->w ? tt->w[tidx]:1);
	  X_t currv = tt->tx[ tlst[vidx+0] ][aidx];
	  X_t nextv = tt->tx[ tlst[vidx+1] ][aidx];
	  if (currv < nextv) {
	    thisInfo = ( calcInfo( &mFreq(0,0), nLabels ) + calcInfo( &mFreq(1,0), nLabels ) );
	    if (vidx==0 || thisInfo < min_info) {
	      min_vidx = vidx;  min_info = thisInfo; 
	      DT2V_SET( wchild, mFreq.getRowSum(0), mFreq.getRowSum(1) );
	    }
	    trial++;
	  }
	  //if (debug) printf("(%.0f) ", tt->tx[ tlst[vidx] ][aidx]); // check the sorting
	}
	if (debug) printf(" trial=%5d   ", trial);
	//if (debug) mFreq.printMatrix("%2.0f");
	if (min_vidx >= 0) {
	  ti->tmpGain[aidx] = calcInfo( mFreqC.data, nLabels ) / wvsum  -  min_info / wvsum;
	  //ti->tmpInfo[aidx] = calcInfo( mFreqV.data, 2 ) / wvsum;
	  ti->tmpInfo[aidx] = calcInfo( wchild, 2 ) / wvsum;
	  ti->tmpAVal[aidx] = (tt->tx[tlst[min_vidx]][aidx] + tt->tx[tlst[min_vidx+1]][aidx]) / 2;
	  double threshCost = Log2( trial ) / nlst;
	  if (debug) printf("%.2f/%.2f/%.2f  ", ti->tmpGain[aidx], ti->tmpInfo[aidx], threshCost);
	  if (ti->tmpGain[aidx] <= threshCost) ti->tmpGain[aidx] = ti->tmpInfo[aidx] = 0;
	} else ti->tmpGain[aidx] = ti->tmpInfo[aidx] = ti->tmpAVal[aidx] = 0;
	if (debug) printf("%.2f/%.2f(min_vidx=%d) \n", ti->tmpGain[aidx], ti->tmpInfo[aidx], min_vidx);
      }
    }
    // if (debug && nlst == tt->nt) tt->printAttrGains();
    // find the best attribute to split with
    int    best_aidx = -1;
    double val=0, best_val=0;
    for (aidx = 0; aidx < sx; aidx++) {
      if (ti->tmpGain[aidx] <= 0 || ti->tmpInfo[aidx] <= 1e-6) continue;
      // For the nodes with multiple children, 'tmpGain[aidx] / tmpInfo[aidx]' is desired.
      val = ti->tmpGain[aidx];
      if (val > best_val) { best_val = val; best_aidx = aidx; }
    }
    if (best_aidx < 0) {
      if (debug) printf(" =>  best attribute NOT found\n");
    } else {
      if (debug) printf(" =>  best(%d)=%g\n", best_aidx, ti->tmpAVal[best_aidx]);
      if (criteria) *criteria = ti->tmpAVal[ best_aidx ];
    }
    free( tlst );
    return best_aidx;
  }
  
  double calcInfo(double freq[], int n) {
    double total=0, sum=0, f=0;
    for (int i=0; i<n; i++) { f = freq[i]; sum += f * Log2(f);  total += f; }
    return (total * Log2(total) - sum);
  }
  
  void sortDataIndices(int lst[], int aidx, int sidx, int eidx) {
    // Sort the array 'lst[]' in the index range [sidx,eidx] (QuickSort)
    if (sidx >= eidx) return;
    int midx = sortDataIndicesPartition( lst, aidx, sidx, eidx );
    sortDataIndices( lst, aidx, sidx, midx );
    sortDataIndices( lst, aidx, midx+1, eidx );
  }
  int  sortDataIndicesPartition(int lst[], int aidx, int sidx, int eidx) {
    X_t  x = tt->tx[ lst[sidx] ][aidx];
    int  i = sidx-1, j = eidx+1;
    while (1) {
      do { j = j - 1; } while ( tt->tx[ lst[j] ][aidx] > x );
      do { i = i + 1; } while ( tt->tx[ lst[i] ][aidx] < x );
      if (i < j) { int tmp = lst[i]; lst[i] = lst[j]; lst[j] = tmp; }
      else return j;
    }
  }
  
  // -----------------------------------------------------------------
  // Pruning
  // -----------------------------------------------------------------
private:
  void prune(int vset[], int vsize, int tset[], int tsize, int szinfo[]) {
    // Reduced Error Pruning with the test data set to prevent overfitting.
    //   vset[vsize] : list of indices to be used for testing (validation set)
    //   vsize       : number of data  to be used for testing (validation set)
    //   tset[tsize] : list of indices to be used for training (training set)
    //   tsize       : number of data  to be used for training (training set)
    //   w[total]    : list of weights for the entire feature vectors (total>vsize)
    double cerror=0, perror=0, max_reduction=0;
    do {
      cerror = evaluateTree( root, vset, vsize, &max_reduction );
      perror = pruneTree( root, max_reduction );
      totalw = sumupTree( root, szinfo, NULL );
      nnodes = szinfo[0];	// number of decision nodes in the tree (output)
      nleaf  = szinfo[1];	// number of leaf nodes in the tree (output)
      mdepth = szinfo[2];	// maximum depth of the tree (output)
      if (verbose && perror < cerror)
	printf("[DTree]   pruned to size (%d %d) (error: %.4f => %.4f with max_reduction=%g)\n", nnodes, nleaf, cerror, perror, max_reduction);
    } while (perror < cerror);
  }
  
  double evaluateTree(DTNode<X_t> *np, int vset[], int vsize, double *max_reduction) {
    // Evaluate the subtree using the test(validation) set, and
    //   return the current error before the pruning.
    PruningData *pdata = np->pdata;
    if (!pdata) { pdata = np->pdata = new PruningData; }
    if (np->type == 'L') {
      pdata->perr = 0;
      for (int i=0; i < vsize; i++) {
	if (tt->ty[ vset[i] ] != np->L.label) pdata->perr += (tt->w ? tt->w[ vset[i] ] : 1.0);
      }
      pdata->cerr = pdata->perr;
    } else {
      int *sublst[2]={NULL,NULL}, i, lpos=0, rpos=0;
      sublst[0] = (int*) calloc( vsize, sizeof(int) );
      sublst[1] = (int*) calloc( vsize, sizeof(int) );
      pdata->perr = 0;
      for (i=lpos=rpos=0; i < vsize; i++) {
	if (ai->getAttr(np->N.aidx)->discrete) {
	  if ( tt->tx[ vset[i] ][ np->N.aidx ] == np->N.value ) sublst[0][ lpos++ ] = vset[i];
	  else                                                  sublst[1][ rpos++ ] = vset[i];
	} else {
	  if ( tt->tx[ vset[i] ][ np->N.aidx ] < np->N.value ) sublst[0][ lpos++ ] = vset[i];
	  else                                                 sublst[1][ rpos++ ] = vset[i];
	}
	if (tt->ty[ vset[i] ] != np->N.flabel) pdata->perr += (tt->w ? tt->w[ vset[i] ] : 1.0);
      }
      pdata->cerr = ( evaluateTree( np->N.child[0], sublst[0], lpos, max_reduction ) + 
		      evaluateTree( np->N.child[1], sublst[1], rpos, max_reduction ) );
      double reduction = pdata->cerr - pdata->perr;
      if (reduction > *max_reduction) *max_reduction = reduction;
    }
    return pdata->cerr;
  }
  
  double pruneTree(DTNode<X_t> *np, double max_reduction) {
    // Prune the subtree where the max_reduction value is found, and
    //   return the new error for the subtree.
    PruningData *pdata = np->pdata;
    if (!np || !pdata || max_reduction<=0) return 1e6;
    double perror=0, reduction = pdata->cerr - pdata->perr;
    if (np->type == 'L') {
      perror = pdata->cerr;
    } else if (reduction == max_reduction) {
      perror = pdata->perr;
      int    nLabels = ai->getLabelCCount(), flabel = np->N.flabel, szinfo[4]={0,0,0,0};
      double *wdist = (double*)calloc( nLabels, sizeof(double) );
      double wsum = sumupTree( np, szinfo, wdist );
      np->clear();
      np->setLeafNode( flabel, wsum, wdist, nLabels );
      np->pdata = new PruningData( perror, perror );
      free( wdist );
      if (verbose)
	printf("[DTree]   prunning subtree of size (%02d + %02d) to a leaf node\n", szinfo[0], szinfo[1]);
    } else {
      perror = ( pruneTree( np->N.child[0], max_reduction ) +
		 pruneTree( np->N.child[1], max_reduction ) );
    }
    return perror;
  }
  
  double sumupTree(DTNode<X_t> *np, int szinfo[4], double *wdist) {
    // 
    if (!np) return 0;
    double wsum = 0;
    if (np->type == 'L') {
      if (np->L.ldist) {
	int  i, nLabels = ai->getLabelCCount();
	for (i=0; i < nLabels; i++) { 
	  if (wdist) wdist[i] += np->L.ldist[i];
	  wsum += np->L.ldist[i];
	}
      } else {
	if (wdist) wdist[ np->L.label ] += np->L.items;
	wsum = np->L.items;
      }
      if (szinfo) DT3V_SET( szinfo, 0, 1, 0 );
    } else {
      int    szL[4]={0,0,0,0}, szR[4]={0,0,0,0};
      double wL = sumupTree( np->N.child[0], szL, wdist );
      double wR = sumupTree( np->N.child[1], szR, wdist );
      wsum = wL + wR;
      if (szinfo) {
	szinfo[0] = szL[0] + szR[0] + 1;	// # of decision nodes in the subtree
	szinfo[1] = szL[1] + szR[1];		// # of leaf nodes in the subtree
	szinfo[2] = (szL[2] > szR[2] ? szL[2]+1 : szR[2]+1);  // max depth
      }
    }
    return wsum;
  }
  
  // -----------------------------------------------------------------
  // Testing
  // -----------------------------------------------------------------
public:
  int   test(X_t x[], double *prob=NULL, double *conf=NULL) {
    // Test the given example, and return the label for it.
    //   prob	: [output] the probability for the output label
    //   conf   : [output] the confidence of the output label (log-likelihood ratio)
    if (!ai->isReady()) { printf("Error (DecisionTree::test): AttributeInfo NOT set\n"); return -1; }
    if (!root) { printf("Error (DecisionTree::test): Tree is NOT ready\n"); return -1; }
    // test the example by calling testRecursive() function which saves the result in prv[]
    //   (This is only for future expansion to support 'unknown' values)
    int nLabels = ai->getLabelCCount();
    memset( prv, 0, nLabels*sizeof(double) );
    testRecursive( x, root, 1.0 );
    // find the best label
    int  i, lbest;  double lblwbest=0, lblwsum=0;
    for (i = lbest = 0; i < nLabels; i++) {
      if (i==0 || prv[i] > prv[lbest]) { lbest = i; lblwbest = prv[i]; }
      lblwsum += prv[i];
    }
    if (lblwsum<=0) return -1;
    if (lblwbest > lblwsum*0.9999) lblwbest = lblwsum * 0.9999;
    if (prob) *prob = lblwbest / lblwsum;                            // probability
    if (conf) *conf = 0.5 * log( lblwbest / (lblwsum - lblwbest) );  // log likelihood ratio
//     if (verbose) {
//       printf("[DTree] TestResult: %d (%s)  { ", lbest, ai->convI2S(-1, lbest) );
//       for (i = 0; i < nLabels; i++) printf("%.4f ", prv[i]);
//       printf("} \n");
//     }
    return lbest;
  }
  
private:
  void testRecursive(X_t x[], DTNode<X_t> *np, double weight) {
    // Test a data (feature vector) and save the probability in 'prv[]'.
    int nLabels = ai->getLabelCCount();
    if        (np->type == 'L') {	// Leaf
      np->visited = true;
      if (np->L.ldist) {
	for (int i = 0; i < nLabels; i++)
	  prv[i] += weight * np->L.ldist[i] / np->L.items;
      } else {
	prv[np->L.label] += weight;
      }
    } else if (np->type == 'D') {	// Discrete
      if (x[np->N.aidx] == np->N.value) testRecursive( x, np->N.child[0], weight ); 
      else                              testRecursive( x, np->N.child[1], weight );
    } else if (np->type == 'C') {	// Continuous
      if (x[np->N.aidx] <  np->N.value) testRecursive( x, np->N.child[0], weight );
      else                              testRecursive( x, np->N.child[1], weight );
    }
  }
  
public:
  int testOnTData(char *cmmt=NULL, int lst[]=NULL, int nlst=0) {
    // Test an expert on the training dataset itself. (for debugging)
    if (!ai || !ai->isReady()) return 0;
    if (!tt || !tt->isReady()) return 0;
    if (nlst <= 0 || !lst) { nlst = tt->nt; lst = NULL; }
    printf("%s(%d) : <", (cmmt ? cmmt : "TestOnTData"), nlst);  fflush(stdout);
    int  i, cnt=0, tidx;  double wsum=0;
    for (i = 0; i < nlst; i++) {
      tidx = (lst ? lst[i] : i);
      bool succ = (this->test( tt->tx[tidx] ) == tt->ty[tidx]);
      if (nlst<=5000) { 
	printf("%c", (succ ? '.':('0'+tt->ty[tidx])));
	if (i%10==9) printf(" ");
      }
      if (succ) { cnt++;  if (tt->w) wsum += tt->w[tidx]; }
    }
    if (nlst>5000) printf(" TooManyToShow ");
    if (tt->w) printf(">  %d%% success (%.0f%% w)\n", cnt*100/nlst, wsum*100);
    else       printf(">  %d%% success\n", cnt*100/nlst);
    return cnt;
  }
  
  // -----------------------------------------------------------------
  // File I/O of the tree
  // -----------------------------------------------------------------
public:
  bool writeFile(const char *fname, const char *cmmt=NULL) {
    FILE *fp = fopen(fname, "w+");
    if (!fp) return false;
    ai->writeAttributeInfo( fp );
    writeFileSingleLine( fp );
    fclose(fp);
    return true;
  }
  bool readFile(const char *fname) {
    FILE *fp = fopen(fname, "r");
    if (!fp) return false;
    if (!ai->readAttributeInfo( fp )) {
      fclose( fp );
      std::cerr << "Error (DecisionTree::readFile): invalid AttributeInfo" << std::endl;
      return false;
    }
    setAttributeInfo();
    bool ret = readFileSingleLine( fp );
    fclose(fp);
    return ret;
  }
  void writeFileSingleLine(FILE *fp, DTNode<X_t> *np=NULL) {
    if (root == NULL) return;
    else if (np == NULL) {
      //if (verbose) printf("DTree: writeFileSingleLine\n");
      fprintf( fp, "DTree(%d+%d,%d): ", nnodes, nleaf, mdepth );
      writeFileSingleLine( fp, root );
      fprintf( fp, "\n" );
    } else {
      if        (np->type == 'L') { 
	fprintf(fp, "(%d!%g:[", (int)(np->L.label), np->L.items); 
	if (np->L.ldist) {
	  for (int i=0; i < ai->getLabelCCount(); i++) 
	    fprintf(fp, "%s%.g", (i==0 ? "":","), (np->L.ldist ? np->L.ldist[i] : -1));
	}
	fprintf(fp, "]) ");
      } else if (np->type == 'D' || np->type == 'C') { 
	if (np->type=='D') fprintf(fp, "(%d=%d) ", np->N.aidx, (int)np->N.value); 
	else               fprintf(fp, "(%d<%g) ", np->N.aidx, np->N.value); 
	if (np->N.child[0]) writeFileSingleLine( fp, np->N.child[0] );
	if (np->N.child[1]) writeFileSingleLine( fp, np->N.child[1] );
      }
    }
  }
  bool readFileSingleLine(FILE *fp) {
    clear();
    char buf[80];
    fscanf( fp, "%s", buf);
    if (sscanf( buf, "DTree(%d+%d,%d):", &nnodes, &nleaf, &mdepth ) != 3) {
      fprintf(stderr, "Error (DecisionTree:readFileSingleLine): invalid contents '%s'\n", buf);
      return false;
    }
    root = readFileSingleLineRecursive( fp );
    if (!root) return false;
    fscanf( fp, "\n" );
    return true;
  }
private:
  DTNode<X_t>* readFileSingleLineRecursive( FILE *fp ) {
    DTNode<X_t> *np = new DTNode<X_t>;
    char buf[160], first='\0';
    fscanf( fp, "%s ", buf );
    int  i, k, blen=strlen(buf), nLabels = ai->getLabelCCount();
    for (i=0; i<blen; i++) {
      char bb[]={ '(', ')', ',', '[', ']', '!', '=', '<' };
      for (k=0; k<8; k++) if (buf[i]==bb[k]) break;
      if (!first && k<8 && buf[i]!='(') first = buf[i];
      if  (k<8)  buf[i] = ' ';
    }
    //if  (debug) printf("read '%s', first='%c'\n", buf, first);
    if  (first=='!') {	// leaf node
      int     idx, wc;  double wsum;
      sscanf( buf, " %d %lf ", &idx, &wsum );
      int     label = idx;
      for (i=1, wc=0; i<blen; i++)
	if (buf[i]!=' ' && buf[i-1]==' ') {
	  if (wc >= 2) prv[wc-2] = atof( buf + i );
	  wc++;
	}
      np->setLeafNode( label, wsum, (wc==nLabels+2 ? prv:NULL), nLabels );
      //if (debug) printf("leaf %d wsum=%.4f  wc=%d/%d\n", label, wsum, wc, nLabels);
    } else {		// decision node 
      int aidx;  double criteria;
      sscanf( buf, " %d %lf ", &aidx, &criteria );
      if (first=='=') {
	np->setDecisionNode( true, aidx, (X_t)criteria );
	np->N.child[0] = readFileSingleLineRecursive( fp );
	np->N.child[1] = readFileSingleLineRecursive( fp );
      } else if (first=='<') {
	np->setDecisionNode( false, aidx, (X_t)criteria );
	np->N.child[0] = readFileSingleLineRecursive( fp );
	np->N.child[1] = readFileSingleLineRecursive( fp );
      } else {
	printf("Warning (DecisionTree::readFileSingleLineRecursive): invalid node type\n");
	delete(np); np = NULL;
      }
      //if (debug) printf("Node [%d] %c %.g ?\n", aidx, first, criteria);
    }
    return np;
  }
  
  // -----------------------------------------------------------------
  // Etc
  // -----------------------------------------------------------------
public:
  char* getSingleLineInfo(char buffer[], DTNode<X_t> *np=NULL) {
    // Get a single line information of the tree.
    if (!np) { sprintf(buffer, "DTree[%2d+%2d,%d] ", nnodes, nleaf, mdepth); np = root; }
    if        (np->type == 'L') {
      sprintf(buffer+strlen(buffer), "(%s) ", ai->convI2S(-1,np->L.label));
    } else if (np->type == 'D' || np->type == 'C') {
      if (np->type == 'D') 
	sprintf(buffer+strlen(buffer), "(%s=%s) ", ai->getAttrName(np->N.aidx), ai->convI2S(np->N.aidx, (int)np->N.value));
      else 
	sprintf(buffer+strlen(buffer), "(%s<%g) ", ai->getAttrName(np->N.aidx), np->N.value);
      if (np==root) getSingleLineInfo( buffer, np->N.child[0] );
      if (np==root) getSingleLineInfo( buffer, np->N.child[1] );
    } else {
      sprintf(buffer, "(?) ");
    }
    return buffer;
  }
  void printTrainingData(char *cmmt=NULL, int tlst[]=NULL, int nlst=0) {
    if (tt && tt->isReady()) tt->printTrainingData(cmmt, tlst, nlst);
  }
public:
  void printInfo(char *cmmt=NULL, bool detail=true) {
    printf("%s with (%d+%d) nodes, maxd=%d, totalw=%g\n", 
	   (cmmt ? cmmt : "DecisionTree"), nnodes, nleaf, mdepth, totalw);
    if (ai && ai->isReady()) ai->printInfo("  Attrb: ", true);
    if (tt && tt->nt>0 && tt->tx)  printf("  TData: %d training data set \n", tt->nt);
    if (root) {
      printf("  DTree: %d decision and %d leaf nodes with maximum depth %d\n", nnodes, nleaf, mdepth);
      if (detail) printInfoRecursive( root, 4 );
    } else printf("  Tree is NOT ready\n");
  }
private:
  void printInfoRecursive(DTNode<X_t> *np, int indent=0) {
    char blank[80];  memset(blank,' ',(indent<79?indent:79));  blank[(indent<79?indent:79)]='\0';
    if (np == NULL) {
      printf("%sError (MLH::DTREE::printInfo): NULL tree node\n", blank);
    }
    if        (np->type == 'L') {
      double pr = (np->L.ldist ? (np->L.ldist[np->L.label]/np->L.items) : 1.0);
      printf(" =>  %s (w=%.2f, p=%.2f)\n", ai->convI2S(-1, np->L.label), np->L.items, pr);
    } else if (np->type == 'D') {
      int aidx = np->N.aidx;
      printf("%s(%s == %s) ", blank, ai->getAttr(aidx)->name, ai->convI2S(aidx, (int)np->N.value));
      if (np->N.child[0]->type=='L') printInfoRecursive( np->N.child[0] );
      else { printf(" ?\n"); printInfoRecursive(np->N.child[0], indent+2); }
      if (ai->getAttr(aidx)->d.ncases>2)
	printf("%s(%s != %s) ", blank, ai->getAttr(aidx)->name, ai->convI2S(aidx, (int)np->N.value));
      else
	printf("%s(%s == %s) ", blank, ai->getAttr(aidx)->name, ai->convI2S(aidx, ((int)np->N.value + 1)%2));
      if (np->N.child[1]->type=='L') printInfoRecursive( np->N.child[1] );
      else { printf(" ?\n"); printInfoRecursive(np->N.child[1], indent+2); }
    } else if (np->type == 'C') {
      printf("%s(%s <  %g) ", blank, ai->getAttr(np->N.aidx)->name, (float)np->N.value);
      if (np->N.child[0]->type=='L') printInfoRecursive( np->N.child[0] ); 
      else { printf(" ?\n"); printInfoRecursive(np->N.child[0], indent+2); }
      printf("%s(%s >= %g) ", blank, ai->getAttr(np->N.aidx)->name, (float)np->N.value);
      if (np->N.child[1]->type=='L') printInfoRecursive( np->N.child[1] );
      else { printf(" ?\n"); printInfoRecursive(np->N.child[1], indent+2); }
    } else printf("%sWarning (DecisionTree::printInfo): invalid node type\n", blank);
  }
  
};  
  
}	// end of namespace MLH

#endif  // MLH_DECISION_TREE_HPP


// ===================================================================
#if 0	// beginning of the example
// ===================================================================

#include <iostream>
#include "../common/mlh_decision_tree.hpp"
using namespace std;
int main(int argc, char **argv)
{
  if (argc<2) return EXIT_FAILURE;
  MLH::DecisionTree<double,bool>	dtree;
  dtree.readTrainingDataFile( argv[1], true );
  dtree.train();
  dtree.printInfo();
  double x1[]={ 0, 80, 60, 0 }, x2[]={ 2, 75, 77, 1}, prob=0;
  printf("Test 1: %s  p=%.4f\n", dtree.getLabelName(dtree.test(x1, &prob)), prob);
  printf("Test 2: %s  p=%.4f\n", dtree.getLabelName(dtree.test(x2, &prob)), prob);
  dtree.writeFile("output1.txt");
  return EXIT_SUCCESS;
}
// Contents of the input text file  --------------
// Input     : 4
//         0 : outlook              (D:3:sunny:overcast:rain)  i=0.00
//         1 : temperature          (C:0:100)  i=0.00
//         2 : humidity             (C:0:100)  i=0.00
//         3 : windy                (D:2:false:true)  i=0.00
// Output    : Play?                (D:2:DoNotPlay:Play)
// sunny, 85, 85, false, DoNotPlay
// sunny, 80, 90, true, DoNotPlay
// overcast, 83, 78, false, Play
// rain, 70, 96, false, Play
// rain, 68, 80, false, Play
// rain, 65, 70, true, DoNotPlay
// overcast, 64, 65, true, Play
// sunny, 72, 95, false, DoNotPlay
// sunny, 69, 70, false, Play
// rain, 75, 80, false, Play
// sunny, 75, 70, true, Play
// overcast, 72, 90, true, Play
// overcast, 81, 75, false, Play
// rain, 71, 80, true, DoNotPlay
// -----------------------------------------------

// ===================================================================
#endif	// end of the example
// ===================================================================
