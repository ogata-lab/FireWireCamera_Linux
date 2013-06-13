
//
// MLH::DTree<X_t> : Decision Tree for Binary Classification
//
// Jaeil Choi
// Last modified in Sep, 2009
//
// A simple decision tree for binary classification problem,
//   (This class was developed as a weak learner for Boosting algorithm.)
//   template data types:
//     X_t        : type of training input data (float or double)
// Member functions
//   - void  clear(void);
//   - void  setAttributeInfo( AttributeInfo<> *ap );
//   - void  setTrainingData ( TrainingData<>  *tp );
//   - bool  train( int label, int maxdepth, int dlist[]=NLL, int dsize=0 );
//   - bool  test ( X_t x[], double prob[2]=NULL, double *confidence=NULL );
//   - int   testOnTData( char *cmmt=NULL, int lst[]=NULL, int nlst=0 );
//   - void  writeFileSingleLine(FILE *fp);
//   - bool  readFileSingleLine(FILE *fp);
//   - char* getSingleLineInfo(char *buffer);
// Member function arguments:
//     ap         : pointer to MLH::AttributeInfo<X_t>, for the attribute information
//     tp         : pointer to MLH::TrainingData<X_t>, for the training data
//     label      : target label to learn for. (Note that training data may have multiple labels)
//     maxdepth   : train(); maximum depth limit of the tree [optional]
//     dlist[]    : train(); list of indices (of tx[] and ty[]) to be used for training [optional]
//     dsize      : train(); size of input index list [optional]
//     w[dsize]   : train(); weights for each of given feature vectors [optional]
//     x[sx]      : test(); a feature vector to be evaluated for test()
//     prob       : test(); calculated probability for the returned label [optional]
//     confidence : test(); log likelihood ratio for the returned label; log(Positive/Negative)/2
//


#ifndef MLH_DTREE_HPP
#define MLH_DTREE_HPP

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
class DTreeTrainInfo {
public:
  int		sx;
  int		maxdepth;	// depth limit
  // temporary buffers for tree building
  double	retdist[3];	// buffer for return value {false,true} distribution of a leaf node
  double	*tmpInfo;
  double	*tmpGain;
  X_t		*tmpAVal;
public:
  DTreeTrainInfo() : tmpInfo(NULL), tmpGain(NULL), tmpAVal(NULL) {}
  DTreeTrainInfo(int max_depth, int sx) : tmpInfo(NULL), tmpGain(NULL), tmpAVal(NULL) { setTrainInfo(max_depth, sx); }
  ~DTreeTrainInfo() { clearTrainInfo(); }
  void setTrainInfo(int max_depth, int sx) {
    // Allocate temporary memory buffers for training
    tmpInfo = (double*) malloc( sx * sizeof(double) );
    tmpGain = (double*) malloc( sx * sizeof(double) );
    tmpAVal = (X_t*) malloc( sx * sizeof(X_t) );
    this->sx = sx;
    this->maxdepth = (max_depth>0 ? max_depth : 65535);
  }
  void clearTrainInfo(void) {
    // Clear temporary memory buffers for training
    if (tmpInfo) free(tmpInfo);  tmpInfo = NULL;
    if (tmpGain) free(tmpGain);  tmpGain = NULL;
    if (tmpAVal) free(tmpAVal);  tmpAVal = NULL;
    this->sx = this->maxdepth = 0;
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
  
class DTPruningData { 
public:
  double cerr, perr; 
public:
  DTPruningData() {}
  DTPruningData(double ce, double pe) : cerr(ce), perr(pe) {}
  ~DTPruningData() {}
};
  
template <class X_t>
class DTreeNode {
public:
  char			type;		// 'L'eaf, 'D'iscrete, or 'C'ontinous
  bool			visited;	// visited by 'test()' function (OnDecisionPath?)
  double		items;		// num of items (or total weights) at this node
  union {
    struct {	// Leaf Node
      bool		ret;		// most frequent class at this node
      double		errors;		// num of errors at this node
      double		rdist[2];	// return value {false,true} distribution
    } L;
    struct {	// Decision Node
      int		aidx;		// attribute index
      X_t		value;		// attribute value for splitting
      DTreeNode<X_t>	*child[2];	// child nodes ( '==' and '!=' or '<' and '>=' )
      float		gain;
      bool		fret;		// most frequent return value {false,true}
    } N;
  };
  DTPruningData		*pdata;
  // 
public:
  DTreeNode() { memset(this, 0, sizeof(DTreeNode<X_t>)); }
  ~DTreeNode() { clear(); }
  void clear(void) {
    if (type=='D' || type=='C') {
      if (N.child[0]) delete(N.child[0]); N.child[0] = NULL; 
      if (N.child[1]) delete(N.child[1]); N.child[1] = NULL;
    }
    if (pdata) { free(pdata); pdata = NULL; }
  }
  void setLeafNode(bool ret, double wsum, double wdist[]) {
    clear();
    this->type = 'L';
    this->items = wsum;
    this->L.ret   = ret;
    if    (wdist) { this->L.rdist[0] = wdist[0];  this->L.rdist[1] = wdist[1]; }
    else if (ret) { this->L.rdist[0] = 0;         this->L.rdist[1] = wsum;     }
    else          { this->L.rdist[0] = wsum;      this->L.rdist[1] = 0;        }
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
class DTree {
  // prerequisite information ----------------------------------------
public:
  AttributeInfo<X_t>	*ai;		// attribute information
  TrainingData<X_t>	*tt;		// training data
private:
  AttributeInfo<X_t>	ainfo;		// attribute information dummy
  TrainingData<X_t>	tdata;		// training data dummy
  // the decision tree -----------------------------------------------
public:
  int			target_label;
  DTreeNode<X_t>	*root;		// the root node
  int			nnodes;		// number of decision nodes
  int			nleaf;		// number of leaf node
  int			mdepth;		// maximum depth
  double		totalw;		// total weight
  bool			pruning;
  bool			verbose;
  // etc -------------------------------------------------------------
private:
  DTreeTrainInfo<X_t>	*ti;		// pointer to training data and buffers
  double		prv[2];		// temporary buffer for tree testing
  
public:
  DTree() : ai(&ainfo), tt(&tdata), root(NULL), nnodes(0), nleaf(0), mdepth(0), 
	    pruning(false), verbose(false), ti(NULL) {}
  ~DTree() { clear(); }
  void clear(void)  {
    if (root) delete(root); root = NULL;  
    if (ti)   delete(ti);   ti = NULL;
    nnodes = nleaf = mdepth = 0; 
  }
  bool isValid(void) { return (root && (root->type=='L' || root->type=='D' || root->type=='C')); }
  inline char* getLabelName(int lbl) { return ai->convI2S( -1, lbl ); }
  
  // -----------------------------------------------------------------
  // Setting Prerequisite Information
  // -----------------------------------------------------------------
public:
  void setAttributeInfo(AttributeInfo<X_t> *ap=NULL) {
    // Use provided attribute information. (It might be shared with others)
    this->ai = ( ap ? ap : &this->ainfo );
    if (!ai->isReady()) printf("Warning (DTree::setAttributeInfo): AttributeInfo is not ready yet\n");
  }
  void setTrainingData(TrainingData<X_t> *tp=NULL) {
    // Use provided training data. (It might be shared with others)
    this->tt = ( tp ? tp : &this->tdata );
    if (!tt->isReady()) printf("Warning (DTree::setTrainingData): TrainingData is not ready yet\n");
  }
  
  // -----------------------------------------------------------------
  // Training
  // -----------------------------------------------------------------
public:
  bool train(int label, int maxdepth, int dlist[]=NULL, int dsize=0) {
    // Create a decision tree using given training datasets, labels, and weights.
    // This function assumes that the training data was already set.
    //   label       : target label to learn for. (Note that there might be multiple labels)
    //   maxdepth    : maximum depth (0: not restricted)
    //   dlist[dsize]: list of indices to be used for training
    //   dsize       : number of data  to be used for training
    // Note that this function assumes 'setTrainingData()' was already called.
    if (!ai || !ai->isReady()) { printf("Error (DTree::train): DataInfo not set\n"); return false; }
    if (!tt || !tt->isReady()) { printf("Error (DTree::train): Training data not set\n"); return false; }
    if (dsize<0||dsize>tt->nt) { printf("Error (DTree::train): invalid argument\n"); return false; }
    clear();	// clear any existing decision nodes
    ti = new DTreeTrainInfo<X_t>( maxdepth, ai->getAttrCount() );
    int  i, szinfo[4], *tset=NULL, tsize=0, *vset=NULL, vsize=0, positive=0;
    this->target_label = label;
    if (dlist) {
      if (!pruning) {   // use the list for the training set, without pruning
	tsize = dsize;
	tset  = dlist;
      } else {          // use the list of the list, with post-pruning
	tsize = vsize = dsize / 2;
	tset  = (int*) malloc( (tsize+1) * sizeof(int) );  // training set
	vset  = (int*) malloc( (vsize+1) * sizeof(int) );  // validation set
	for (i=0; i<dsize; i++)
	  if (i%2==0) tset[i>>1] = dlist[i]; else vset[i>>1] = dlist[i];
      }
    } else {
      if (!pruning) {   // use everything in the training set, without pruning
	tsize = tt->nt;
	tset  = (int*) malloc( tsize * sizeof(int) );  // training set
	for (i=0; i<tt->nt; i++) tset[i] = i;
      } else {          // use just half of the training set, with post-pruning
	tsize = vsize = tt->nt / 2;
	tset  = (int*) malloc( (tsize+1) * sizeof(int) );  // training set
	vset  = (int*) malloc( (vsize+1) * sizeof(int) );  // validation set
	for (i=0; i<tt->nt; i++) 
	  if (i%2==0) tset[i>>1] = i; else vset[i>>1] = i;
      }
    }
    for (i=0; i<tsize; i++) if (tt->ty[ tset[i] ]==target_label) positive++;
    if (verbose) printf("[DTree] training on %d data (class%d: %d), weight=%s pruning=%s\n", 
			tsize, target_label, positive, (tt->w ? "Y":"N"), (pruning ? "Y":"N"));
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
    if ((!dlist || pruning) && tset) { free(tset); tset = NULL; tsize = 0; }
    if ((          pruning) && vset) { free(vset); vset = NULL; vsize = 0; }
    if (ti) { delete(ti); ti = NULL; }
    updateNodeTotalItems( root );
    return true;
  }
  
private:
  DTreeNode<X_t>* createDecisionNode(int *lst, int nlst, int szinfo[4]) {
    //    lst     : index list of the training data
    //    nlst    : size of the index list
    //    szinfo[0] : number of decision nodes in this subtree (botton-up)
    //    szinfo[1] : number of leaf nodes in the subtree (bottom-up)
    //    szinfo[2] : maximum depth of this subtree (bottom-up)
    //    szinfo[3] : depth of current node (top-down)
    if (lst == NULL || nlst < 1) return NULL;
    char blank[80];  makeBlank(blank, szinfo[3]);
    // create a new decision node
    DTreeNode<X_t> *np = new DTreeNode<X_t>;
    // setup the best decision criteria (and its subsets)
    int    i, aidx=0;
    bool   ret = getMostLikelyResult( lst, nlst );
    double wsum = ti->retdist[2], wchild[2];
    X_t    criteria=0;
    bool   leaf = (ti->retdist[(ret ? 1:0)] == wsum || nlst<=1 || szinfo[3] >= ti->maxdepth);
    bool   unable = (!leaf && (aidx = findBestAttribute( lst, nlst, wsum, &criteria, wchild )) < 0);
    if (wsum <= 0) {			// ERROR ---------------------
      delete(np); return NULL;
    } else if (leaf || unable) { 	// LEAF ----------------------
      if (ti->retdist[(ret ? 1:0)]==wsum) np->setLeafNode( ret, wsum, NULL );
      else                                np->setLeafNode( ret, wsum, ti->retdist );
      if (verbose) {
	printf("%sret = %d  with items=%.2f [%.2f %.2f]\n", blank, np->L.ret, np->items, np->L.rdist[0], np->L.rdist[1]);
	if (unable) {
	  printf("Error (DTree::createDecisionNode): unable to find the best attribute\n");
	  //tt->printTrainingData( "Unable", lst, nlst );
	}
      }
      DT3V_SET( szinfo, 0, 1, 0 );
    } else {				// Decision node -------------
      int *sublst[2]={NULL,NULL}, lpos, rpos, nsize[2][4];
      sublst[0] = (int*) calloc( nlst, sizeof(int) );
      sublst[1] = (int*) calloc( nlst, sizeof(int) );
      DT4V_SET( nsize[0], 0, 0, 0, szinfo[3]+1 );
      DT4V_SET( nsize[1], 0, 0, 0, szinfo[3]+1 );
      bool discrete = ai->getAttr(aidx)->discrete;
      if (discrete) {			// DISCRETE ------------------
	np->setDecisionNode( true, aidx, criteria );
	// separate the list into two sets
	for (i=lpos=rpos=0; i < nlst; i++) {
	  X_t value = tt->tx[ lst[i] ][ np->N.aidx ];
	  if (value == UNKNOWN) { 
	    sublst[0][ lpos++ ] = lst[i]; 
	    sublst[1][ rpos++ ] = lst[i];
	    if (tt->w) tt->w[ lst[i] ] /= 2;
	  } else if (value == np->N.value) sublst[0][ lpos++ ] = lst[i];
	  else                             sublst[1][ rpos++ ] = lst[i];
	}
      } else {				// CONTINUOUS ----------------
	np->setDecisionNode( false, aidx, criteria );
	// separate the list into two sets
	for (i=lpos=rpos=0; i < nlst; i++) {
	  X_t value = tt->tx[ lst[i] ][ np->N.aidx ];
	  if (value == UNKNOWN) {
	    sublst[0][ lpos++ ] = lst[i]; 
	    sublst[1][ rpos++ ] = lst[i];
	    if (tt->w) tt->w[ lst[i] ] /= 2;
	  } else if ( value < np->N.value ) sublst[0][ lpos++ ] = lst[i];
	  else                              sublst[1][ rpos++ ] = lst[i];
	}
      }
      if (verbose) {
	printf("%snode: ( att[%d] %c %g ) => cnt(%d %d)  ", blank, np->N.aidx, (discrete ? '=':'<'), np->N.value, lpos, rpos);
	ti->printAttrGains("", true, false);
      }
      np->items      = wsum;
      np->N.fret     = ret;			// save the most frequent label of this subtree
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
      for (i=0; i < nlst && tt->w; i++)	// restore the original weight values for 'UNKNOWN'
	if (tt->tx[ lst[i] ][ np->N.aidx ] == UNKNOWN) tt->w[ lst[i] ] *= 2;
      ai->getAttr( aidx )->importance += wsum;
    }
    return np;
  }
  
  void makeBlank(char buf[], int i, int s=2) { memset(buf,' ',(i+1)*s);  buf[(i+1)*s]='\0'; }
  
  // -----------------------------------------------------------------
  // finding the optimal decision (question & criteria)
  // -----------------------------------------------------------------
private:
  bool getMostLikelyResult(int *lst, int nlst) {
    double *retdist = ti->retdist;
    int    i, cidx, *ty = tt->ty;
    retdist[0] = retdist[1] = retdist[2] = 0;
    for (i = 0; i < nlst; i++) {
      cidx = lst[i];
      if (cidx < 0 || cidx > tt->nt) printf("error 1\n");
      double wv = ( tt->w ? tt->w[ cidx ] : 1.0 );
      if (ty[ cidx ] == target_label) retdist[ 1 ] += wv;  // true
      else                            retdist[ 0 ] += wv;  // false
      retdist[ 2 ] += wv;
    }
    return (retdist[0] > retdist[1] ? false : true);
  }
  
  int findBestAttribute(int *lst, int nlst, double wsum, X_t *criteria, double wchild[2]) {
    int i, tidx, aidx, sx = ai->getAttrCount();
    bool debug = false; // (nlst==3 && lst && lst[0]==8);
    DT2V_SET( wchild, 0, 0 );
    if (debug) printf("[DTree] findBestAttribute for %d data among %d attributes: \n", nlst, sx);
    int  *tlst = (int*) malloc( nlst * sizeof(int) );	// 'lst[]' should not be shuffled
    memcpy( tlst, lst, nlst * sizeof(int) );		// because of the sync with 'w[]'
    for (aidx = 0; aidx < sx; aidx++) {
      Attribute *att = ai->getAttr( aidx );
      if (att->discrete) {
	int nvalues = att->d.ncases;
	MTX::Matrix<double> mFreq(nvalues, 2, true), mFreqC, mPosi, mNega;
	// calculate frequency for {false,true} for each value of this attribute
	for (i=0; i<nlst; i++) {		// mFreq is (ncases x 2)
	  tidx = tlst[i];
	  int ridx = (tt->ty[tidx] == target_label ? 1 : 0);
	  mFreq((int)tt->tx[tidx][aidx], ridx) += (tt->w ? tt->w[tidx]:1); 
	}
	mFreqC.sumOverRow( mFreq );		// (1 x 2) frequency for {false,true}
	// find the best attribute value to split with
	double thisInfo, min_info=0;
	int    vidx, min_vidx=-1;
	for (vidx = 0; vidx < nvalues; vidx++) {
	  mPosi.set( 1, 2, &mFreq(vidx,0) );	// (1 x 2) for positive(<=) case
	  mNega.sub( mFreqC, mPosi );		// (1 x 2) for negative(> ) cases
	  thisInfo = ( calcInfo( mPosi.data, 2 ) + calcInfo( mNega.data, 2 ) );
	  if (vidx==0 || thisInfo < min_info) {
	    min_vidx = vidx;  min_info = thisInfo; 
	    DT2V_SET( wchild, mPosi.getRowSum(0), mNega.getRowSum(0) );
	  }
	}
	// remember the gain of the attribute
	ti->tmpGain[aidx] = calcInfo( mFreqC.data, 2 ) / wsum  -  min_info / wsum;
	ti->tmpInfo[aidx] = calcInfo( wchild, 2 ) / wsum;
	ti->tmpAVal[aidx] = (X_t)min_vidx;
	if (debug) printf("  att[%2d] : %.2f/%.2f \n", aidx, ti->tmpGain[aidx], ti->tmpInfo[aidx]);
      } else {
	double freq[2][2]={{0,0},{0,0}}, freqc[2]={0,0};
	MTX::Matrix<double> mFreq(2, 2, true), mFreqC(1, 2, true);
	if (debug) { printf("  att[%2d] : ", aidx); fflush(stdout); }
	for (i=0; i<nlst; i++) {
	  int  ridx = (tt->ty[tlst[i]] == target_label ? 1 : 0);
	  freqc[ridx] += (tt->w ? tt->w[tlst[i]]:1); 
	}
	// sort the list 'tlst[]' wrt the attribute values (ascending order)
	sortDataIndices( tlst, aidx, 0, nlst-1 );
	// find the best attribute cut
	DT2V_COPY( freq[1], freqc );
	double thisInfo, min_info=0;
	int    vidx, min_vidx=-1, trial=0;
	for (vidx = 0; vidx < nlst-1; vidx++) {
	  tidx = tlst[vidx];
	  int ridx = (tt->ty[tidx] == target_label ? 1 : 0);
	  freq[0][ridx] += (tt->w ? tt->w[tidx]:1);
	  freq[1][ridx] -= (tt->w ? tt->w[tidx]:1);
	  X_t currv = tt->tx[ tlst[vidx+0] ][aidx];
	  X_t nextv = tt->tx[ tlst[vidx+1] ][aidx];
	  if (currv < nextv) {
	    thisInfo = ( calcInfo( freq[0], 2 ) + calcInfo( freq[1], 2 ) );
	    if (vidx==0 || thisInfo < min_info) {
	      min_vidx = vidx;  min_info = thisInfo; 
	      DT2V_SET( wchild, (freq[0][0]+freq[0][1]), (freq[1][0]+freq[1][1]) );
	    }
	    trial++;
	  }
	  //if (debug) printf("(%.0f) ", tt->tx[ tlst[vidx] ][aidx]); // check the sorting
	}
	if (debug) { printf(" trial=%5d   ", trial); fflush(stdout); }
	if (min_vidx >= 0) {
	  ti->tmpGain[aidx] = calcInfo( freqc, 2 ) / wsum  -  min_info / wsum;
	  ti->tmpInfo[aidx] = calcInfo( wchild, 2 ) / wsum;
	  ti->tmpAVal[aidx] = (tt->tx[tlst[min_vidx]][aidx] + tt->tx[tlst[min_vidx+1]][aidx]) / 2;
	  double threshCost = Log2( trial ) / nlst;
	  if (debug) { printf("%.2f/%.2f/%.2f  ", ti->tmpGain[aidx], ti->tmpInfo[aidx], threshCost); fflush(stdout); }
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
  
  double evaluateTree(DTreeNode<X_t> *np, int vset[], int vsize, double *max_reduction) {
    // Evaluate the subtree using the test(validation) set, and
    //   return the current error before the pruning.
    DTPruningData *pdata = np->pdata;
    if (!pdata) { pdata = np->pdata = new DTPruningData; }
    if (np->type == 'L') {
      pdata->perr = 0;
      for (int i=0; i < vsize; i++) {
	if ( np->L.ret == true  && tt->ty[ vset[i] ] == target_label ||
	     np->L.ret == false && tt->ty[ vset[i] ] != target_label )
	  pdata->perr += (tt->w ? tt->w[ vset[i] ] : 1.0);
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
	if ( np->N.fret == true  && tt->ty[ vset[i] ] == target_label ||
	     np->N.fret == false && tt->ty[ vset[i] ] != target_label )
	  pdata->perr += (tt->w ? tt->w[ vset[i] ] : 1.0);
      }
      pdata->cerr = ( evaluateTree( np->N.child[0], sublst[0], lpos, max_reduction ) + 
		      evaluateTree( np->N.child[1], sublst[1], rpos, max_reduction ) );
      double reduction = pdata->cerr - pdata->perr;
      if (reduction > *max_reduction) *max_reduction = reduction;
    }
    return pdata->cerr;
  }
  
  double pruneTree(DTreeNode<X_t> *np, double max_reduction) {
    // Prune the subtree where the max_reduction value is found, and
    //   return the new error for the subtree.
    DTPruningData *pdata = np->pdata;
    if (!np || !pdata || max_reduction<=0) return 1e6;
    double perror=0, reduction = pdata->cerr - pdata->perr;
    if (np->type == 'L') {
      perror = pdata->cerr;
    } else if (reduction == max_reduction) {
      perror = pdata->perr;
      int    flabel = np->N.fret, szinfo[4]={0,0,0,0};
      double wdist[2];
      double wsum = sumupTree( np, szinfo, wdist );
      np->clear();
      np->setLeafNode( flabel, wsum, wdist );
      np->pdata = new DTPruningData( perror, perror );
      if (verbose)
	printf("[DTree]   prunning subtree of size (%02d + %02d) to a leaf node\n", szinfo[0], szinfo[1]);
    } else {
      perror = ( pruneTree( np->N.child[0], max_reduction ) +
		 pruneTree( np->N.child[1], max_reduction ) );
    }
    return perror;
  }
  
  double sumupTree(DTreeNode<X_t> *np, int szinfo[4], double wdist[2]) {
    // 
    if (!np) return 0;
    double wsum = 0;
    if (np->type == 'L') {
      if (szinfo) DT3V_SET( szinfo, 0, 1, 0 );
      if (wdist) { wdist[0] += np->L.rdist[0]; wdist[1] += np->L.rdist[1]; }
      wsum = np->L.rdist[0] + np->L.rdist[1];
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
  bool  test(X_t x[], double prob[2]=NULL, double *conf=NULL) {
    // Test the given example, and return the label for it.
    //   prob	: [output] the probability for {false,true}
    //   conf   : [output] the confidence of the output label (log-likelihood ratio)
    if (!ai->isReady()) { printf("Error (DTree::test): AttributeInfo NOT set\n"); return -1; }
    if (!root) { printf("Error (DTree::test): Tree is NOT ready\n"); return -1; }
    // test the example by calling testRecursive() function which saves the result in prv[]
    //   (This is only for future expansion to support 'unknown' values)
    prv[0] = prv[1] = 0;
    testRecursive( x, root, 1.0 );
    double prsum = prv[0] + prv[1];
    if (prob) { prob[0] = prv[0] / prsum;  prob[1] = prv[1] / prsum; }
    if (prv[0] > prv[1]) {
      if (prv[0] > prsum*0.9999) prv[0] = prsum*0.9999;
      if (conf) *conf = 0.5 * log( prv[0]/(prsum-prv[0]) );	// log likelihood ratio
      //if (verbose)  printf("[DTree] TestResult: 'N' { %.4f %.4f }\n", prv[0], prv[1] );
      return false;
    } else {
      if (prv[1] > prsum*0.9999) prv[1] = prsum*0.9999;
      if (conf) *conf = 0.5 * log( prv[1]/(prsum-prv[1]) );	// log likelihood ratio
      //if (verbose)  printf("[DTree] TestResult: 'Y' { %.4f %.4f }\n", prv[0], prv[1] );
      return true;
    }
  }
  
private:
  void testRecursive(X_t x[], DTreeNode<X_t> *np, double weight) {
    // Test a data (feature vector) and save the probability in 'prv[]'.
    if (!np) return;
    else if   (np->type == 'L') {		// Leaf node
      np->visited = true;
      prv[0] += weight * np->L.rdist[0] / np->items;
      prv[1] += weight * np->L.rdist[1] / np->items;
    } else if (x[np->N.aidx] == UNKNOWN) {	// Unknown value
      double wch0 = (np->N.child[0] ? np->N.child[0]->items : 0);
      double wch1 = (np->N.child[1] ? np->N.child[1]->items : 0);
      testRecursive( x, np->N.child[0], weight * wch0 / np->items ); 
      testRecursive( x, np->N.child[1], weight * wch1 / np->items ); 
    } else if (np->type == 'D') {		// Discrete value
      if (x[np->N.aidx] == np->N.value) testRecursive( x, np->N.child[0], weight ); 
      else                              testRecursive( x, np->N.child[1], weight );
    } else if (np->type == 'C') {		// Continuous value
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
      bool ret  = this->test( tt->tx[tidx] );
      bool succ = ( ret == true  && tt->ty[tidx]==target_label ||
		    ret == false && tt->ty[tidx]!=target_label );
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
      std::cerr << "Error (DTree::readFile): invalid AttributeInfo" << std::endl;
      return false;
    }
    setAttributeInfo();
    bool ret = readFileSingleLine( fp );
    fclose(fp);
    return ret;
  }
  void writeFileSingleLine(FILE *fp, DTreeNode<X_t> *np=NULL) {
    if (root == NULL) return;
    else if (np == NULL) {
      //if (verbose) printf("DTree: writeFileSingleLine\n");
      fprintf( fp, "DTree(%02d+%02d,%d): ", nnodes, nleaf, mdepth );
      writeFileSingleLine( fp, root );
      fprintf( fp, "\n" );
    } else {
      if        (np->type == 'L') { 
	fprintf(fp, "(%s!%g:[%g,%g]) ", (np->L.ret ? "Y":"N"), np->items, np->L.rdist[0], np->L.rdist[1]); 
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
      fprintf(stderr, "Error (DTree:readFileSingleLine): invalid contents '%s'\n", buf);
      return false;
    }
    root = readFileSingleLineRecursive( fp );
    if (!root) return false;
    fscanf( fp, "\n" );
    updateNodeTotalItems( root );
    return true;
  }
private:
  DTreeNode<X_t>* readFileSingleLineRecursive( FILE *fp ) {
    DTreeNode<X_t> *np = new DTreeNode<X_t>;
    char buf[160], first='\0';
    fscanf( fp, "%s ", buf );
    int  i, k, blen=strlen(buf);
    for (i=0; i<blen; i++) {
      char bb[]={ '(', ')', ',', '[', ']', '!', '=', '<' };
      for (k=0; k<8; k++) if (buf[i]==bb[k]) break;
      if (!first && k<8 && buf[i]!='(') first = buf[i];
      if  (k<8)  buf[i] = ' ';
    }
    //if  (debug) printf("read '%s', first='%c'\n", buf, first);
    if  (first=='!') {	// leaf node
      int     wc;  char rch;  double wsum;
      sscanf( buf, " %c %lf ", &rch, &wsum );
      for (i=1, wc=0; i<blen; i++)
	if (buf[i]!=' ' && buf[i-1]==' ') {
	  if (wc >= 2) prv[wc-2] = atof( buf + i );
	  wc++;
	}
      bool dec = (rch=='0'||rch=='1' ? (rch-'0' == target_label) :
		  rch=='Y' ? true : false);
      np->setLeafNode( dec, wsum, prv );
      //if (debug) printf("leaf %s wsum=%.4f  wc=%d\n", (ret ? "Y":"N"), wsum, wc);
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
	printf("Warning (DTree::readFileSingleLineRecursive): invalid node type\n");
	delete(np); np = NULL;
      }
      //if (debug) printf("Node [%d] %c %.g ?\n", aidx, first, criteria);
    }
    return np;
  }
  
  double updateNodeTotalItems(DTreeNode<X_t> *np) {
    if (!np) return 0;
    else if (np->type == 'L') return np->items;
    else if (np->type == 'N') {
      np->items = ( updateNodeTotalItems( np->N.child[0] ) +
		    updateNodeTotalItems( np->N.child[1] ) );
      return np->items;
    } else return 0;
  }
  
  // -----------------------------------------------------------------
  // Etc
  // -----------------------------------------------------------------
public:
  char* getSingleLineInfo(char buffer[], DTreeNode<X_t> *np=NULL) {
    // Get a single line information of the tree.
    if (!np) { sprintf(buffer, "DTree[%2d+%2d,%d] ", nnodes, nleaf, mdepth); np = root; }
    if        (np->type == 'L') {
      sprintf(buffer+strlen(buffer), "(%s) ", (np->L.ret ? "Y":"N"));
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
	   (cmmt ? cmmt : "BinaryDTree"), nnodes, nleaf, mdepth, totalw);
    if (ai && ai->isReady()) ai->printInfo("  Attrb: ", true);
    if (tt && tt->nt>0 && tt->tx)  printf("  TData: %d training data set \n", tt->nt);
    if (root) {
      printf("  TAttr: trained for label %d (%s)\n", target_label, ai->convI2S(-1,target_label));
      printf("  DTree: %d decision and %d leaf nodes with maximum depth %d\n", nnodes, nleaf, mdepth);
      if (detail) printInfoRecursive( root, 4 );
    } else printf("  Tree is NOT ready\n");
  }
private:
  void printInfoRecursive(DTreeNode<X_t> *np, int indent=0) {
    char blank[80];  memset(blank,' ',(indent<79?indent:79));  blank[(indent<79?indent:79)]='\0';
    if (np == NULL) {
      printf("%sError (MLH::DTREE::printInfo): NULL tree node\n", blank);
    }
    if        (np->type == 'L') {
      double pr = (np->L.rdist[np->L.ret]/np->items);
      printf(" =>  %s (w=%.2f, p=%.2f)\n", (np->L.ret ? "Y":"N"), np->items, pr);
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
    } else printf("%sWarning (DTree::printInfo): invalid node type\n", blank);
  }
  
};  
  
}	// end of namespace MLH

#endif  // MLH_DTREE_HPP


// ===================================================================
#if 0	// beginning of the example
// ===================================================================

#include <iostream>
#include "../common/mlh_dtree.hpp"
using namespace std;
int main(int argc, char **argv)
{
  if (argc<2) return EXIT_FAILURE;
  MLH::DTree<double>	dtree;
  dtree.readTrainingDataFile( argv[1], true );
  dtree.train();
  dtree.printInfo();
  double x1[]={ 0, 80, 60, 0 }, x2[]={ 2, 75, 77, 1}, prob=0;
  printf("Test 1: %s  p=%.4f\n", dtree.getLabelName(dtree.test(x1, &prob)), prob);
  printf("Test 2: %s  p=%.4f\n", dtree.getLabelName(dtree.test(x2, &prob)), prob);
  dtree.writeFile("output1.txt");
  return EXIT_SUCCESS;
}
// Contents of the input text file  'golf.data' --------------
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
