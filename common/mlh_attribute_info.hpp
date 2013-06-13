
//
// MLH::AttributeInfo<X_t>
//
// Jaeil Choi
// Last modified in Sep, 2009
//


#ifndef MLH_ATTRIBUTE_INFO_HPP
#define MLH_ATTRIBUTE_INFO_HPP

#include <iostream>
#include <cmath>
#include <cstdarg>

namespace MLH {
  
typedef struct {
  char name[31];
  bool discrete;
  union {
    struct { int ncases; char **case_names; } d;
    struct {  float minv;   float maxv;   } c;
  };
  double importance;	// importantan measure after the training
} Attribute;
  
// const int   UNKNOWN=(-59999);
const float UNKNOWN=(-3.39183647E+38F);
  
template <class X_t>
class AttributeInfo {
private:
  int		sx;
  Attribute	*dta;	// description of each data entry
  Attribute	dtl;	// description of the label
  int		aidx;	// current attribute index 
public:
  AttributeInfo() : sx(0), dta(NULL), aidx(0) { memset(&dtl,0,sizeof(Attribute)); dtl.discrete = true; }
  ~AttributeInfo() { clear(); }
  void  clear(void) { 
    for (int i=0; i<sx && dta; i++) {		// clear attribute information
      if (!dta[i].discrete) continue;
      if (dta[i].d.case_names) {
	for (int j=0; j<dta[i].d.ncases; j++) 
	  if (dta[i].d.case_names[j]) free(dta[i].d.case_names[j]);
	free(dta[i].d.case_names);
	dta[i].d.case_names = NULL;
      }
      dta[i].discrete = false;
    }
    if (dtl.d.case_names) {			// clear label information
      for (int j=0; j<dtl.d.ncases; j++) 
	if (dtl.d.case_names[j]) free(dtl.d.case_names[j]);
      free(dtl.d.case_names);
      dtl.d.case_names = NULL;
    }
    if (dta) free(dta);  dta = NULL; sx = 0;
  }
  inline bool       isReady(void) { return (sx>0 && dtl.d.ncases>=2 && sx==aidx); }
  inline int        getAttrCount(void) { return sx; }
  inline Attribute* getAttr(int idx) { return dta+idx; }
  inline char*      getAttrName(int idx) { return dta[idx].name; }
  inline int        getAttrCCount(int idx) { return dta[idx].d.ncases; }
  inline Attribute* getLabel(void)  { return &dtl; }
  inline char*      getLabelName(void)  { return dtl.name; }
  inline int        getLabelCCount(void) { return dtl.d.ncases; }
  void clearImportance(void) { for (int i=0; i<sx && dta; i++) dta[i].importance = 0; }
  
public:
  bool initAttributes(int sx) {
    clear();
    this->sx = sx;
    this->aidx = 0;
    dta = (Attribute*)calloc(sx,sizeof(Attribute));
    return (dta != NULL);
  }
  bool setAttributeDisc(const char *name, int ncases, char **case_names ) { 
    strncpy(dta[aidx].name, name, 19);
    dta[aidx].discrete = true;
    dta[aidx].d.ncases = ncases;
    dta[aidx].d.case_names = case_names;
    aidx++;
    return true;
  }
  bool setAttributeDisc(const char *name, int ncases, const char *case_name, ... ) { 
    if (!dta || aidx < 0 || aidx >= sx || ncases<2 || ncases > 100) return false;
    char **case_names = (char**)calloc( ncases, sizeof(char*) );
    if (!case_names) return false;
    case_names[0] = (char*)calloc( (strlen(case_name)+1), sizeof(char) );
    strcpy( case_names[0], case_name );
    va_list valist;
    va_start( valist, case_name );
    for (int i=1; i < ncases; i++) {
      char *cp = (char*)va_arg( valist, char* );
      case_names[i] = (char*)calloc( (strlen(cp)+1), sizeof(char) );
      strcpy( case_names[i], cp );
    }
    va_end( valist );
    return setAttributeDisc( name, ncases, case_names );
  }
  bool setAttributeCont(const char *name, float minv, float maxv) { 
    if (!dta || aidx < 0 || aidx >= sx) return false;
    strncpy(dta[aidx].name, name, 19);
    dta[aidx].discrete = false;
    dta[aidx].c.minv = (float)minv;
    dta[aidx].c.maxv = (float)maxv;
    aidx++;
    return true;
  }
  bool setOutputLabel(const char *name, int ncases, char **case_names ) { 
    strncpy(dtl.name, name, 19);
    dtl.discrete = true;
    dtl.d.ncases = ncases;
    dtl.d.case_names = case_names;
    return true;
  }
  bool setOutputLabel(const char *name, int ncases, const char *case_name, ... ) { 
    if (ncases<2 || ncases > 100) return false;
    char **case_names = (char**)calloc( ncases, sizeof(char*) );
    if (!case_names) return false;
    case_names[0] = (char*)calloc( (strlen(case_name)+1), sizeof(char) );
    strcpy( case_names[0], case_name );
    va_list valist;
    va_start( valist, case_name );
    for (int i=1; i < ncases; i++) {
      char *cp = (char*)va_arg( valist, char* );
      case_names[i] = (char*)calloc( (strlen(cp)+1), sizeof(char) );
      strcpy( case_names[i], cp );
    }
    va_end( valist );
    return setOutputLabel( name, ncases, case_names );
  }
  
  char* convI2S(int aidx, int cidx=-1) {
    if (aidx < 0) {	// label
      if (cidx < 0) return (char*)(dtl.name);
      if (!dtl.discrete || cidx >= dtl.d.ncases) return "?";
      return (char*)(dtl.d.case_names[cidx]);
    } else {		// attribute
      if (aidx >= sx) return "?";
      if (cidx < 0) return (char*)(dta[aidx].name);
      if (!dta[aidx].discrete || cidx >= dta[aidx].d.ncases) return "?";
      return (char*)(dta[aidx].d.case_names[cidx]);
    }
  }
  int   convS2I(int aidx, const char *cname) {
    if (aidx < 0) {	// label
      if (!dtl.discrete) return -1;
      for (int c=0; c<dtl.d.ncases; c++) 
	if (strcmp(cname, dtl.d.case_names[c])==0) return c;
      return -1;
    } else {		// attribute
      if (aidx >= sx || !dta[aidx].discrete) return -1;
      for (int c=0; c<dta[aidx].d.ncases; c++) 
	if (strcmp(cname, dta[aidx].d.case_names[c])==0) return c;
      return -1;
    }
  }
  
public:
  void printInfo(char *cmmt=NULL, bool single_line=false) {
    if (!isReady()) return;
    else if (single_line) {
      printf("%s(%d: ", (cmmt? cmmt : "AttributeInfo: "), sx); fflush(stdout);
      for (int j=0; j<sx; j++) printf("%s ", getAttr(j)->name);
      printf(") => (%s)\n", getLabel()->name);
    } else {
      printf("%s (sx=%d  aidx=%d)\n", (cmmt ? cmmt : "AttributeInfo"), sx, aidx);
      writeAttributeInfo( stdout, "  " );
    }
  }
  void printAttributeValues(X_t *attr_values, const char *format="%8.4f", double label=-1, int beginning_from=0) {
    // Print the attribute values with given format and label,
    //   beginning from the specified attribute.
    if (!isReady() || !attr_values) return;
    int  j, size=8, prec=4;
    char anfmt[10], avfmt[10];
    if (!format || sscanf(format, "%%%d.%d", &size, &prec)!=2) { size=8; prec=4; }
    sprintf(anfmt, "%%%ds ", size);
    sprintf(avfmt,  "%%%d.%df ", size, prec);
    if (label>=0) printf("  %-12s : ", getLabelName());
    else          printf("  AttrValues : ");
    for (j=beginning_from; j < sx; j++) printf(anfmt, getAttrName(j));	// attribute name
    printf("\n");
    if (label>=0) printf("  %.2f %7s : ", label, dtl.d.case_names[ (int)(label+0.5) ]);
    else          printf("             : ");
    for (j=beginning_from; j < sx; j++) {
      if (attr_values[j]==UNKNOWN) printf(anfmt, " ? ");
      else printf(avfmt, attr_values[j]);	// attribute value
    }
    printf("\n");
  }
  
  // -----------------------------------------------------------------
  // File I/O of the data format
  // -----------------------------------------------------------------
public:
  void writeAttributeInfo(FILE *fp, char *head=NULL) {
    if (!fp) return;
    int  i, c, natt = getAttrCount();
    fprintf( fp, "%sInput     : %d\n", (head ? head:""), natt);
    for (i=0; i<natt; i++) {
      Attribute *att = getAttr(i);
      fprintf( fp, "%s      %3d : %-20s ", (head ? head:""), i, att->name);
      if (att->discrete) {
	fprintf( fp, "(D:%d", att->d.ncases);
	for (c=0; c<att->d.ncases; c++) fprintf( fp, ":%s", convI2S( i, c ) );
	fprintf( fp, ")  i=%.4f\n", att->importance);
      } else {
	fprintf( fp, "(C:%g:%g)  i=%.4f\n", att->c.minv, att->c.maxv, att->importance);
      }
    }
    fprintf( fp, "%sOutput    : %-20s (D:%d", (head ? head:""), dtl.name, dtl.d.ncases );
    for (c=0; c<dtl.d.ncases; c++) fprintf( fp, ":%s", convI2S( -1, c ) );
    fprintf( fp, ")\n");
  }
  bool readAttributeInfo(FILE *fp) {
    int i, natt, k;  char str[512];
    do {
      if (!fgets( str, 512-1, fp )) { return false; }
      if (sscanf( str, "Input     : %d", &natt ) == 1) break;
      else if (str[0] == '#') { if(str[strlen(str)-1]!='\n') fgets(str, 512-1, fp); }
      else { return false; }
    } while (true);
    // if (fscanf( fp, "Input(%d): ", &natt ) != 1) return false;
    initAttributes( natt );
    for (i=0; i<natt; i++) {
      if (fscanf( fp, " %d : ", &k ) != 1) break;
      if (!readFileAttrInfo( fp, i )) break;
    }
    if (i < natt) { clear(); return false; }
    fscanf(fp, "\nOutput    : ");
    if (!readFileAttrInfo( fp, -1 )) { clear(); return false; }
    fscanf(fp, "\n");
    return true;
  }
private:
  bool readFileAttrInfo(FILE *fp, int aidx) {
    char str[512], name[80], sd[4], sn1[40], sn2[40], cname[80];
    if (!fgets( str, 500, fp )) { clear(); return false; }
    int  i, slen=strlen(str), pos;
    for (i=0; i<slen; i++) if (str[i]=='(' || str[i]==')' || str[i] == ':') str[i] = ' ';
    if (sscanf( str, "%s %s %s %s", name, sd, sn1, sn2) != 4) { clear(); return false; }
    if (sd[0] == 'D') {
      int  ncases = atoi(sn1);
      char **cnames = (char**) calloc(ncases, sizeof(char*));
      for (pos=i=0; pos<slen && i<3; pos++) { if (str[pos]==' '&&str[pos+1]!=' ') i++; }
      for (i=0; i<ncases && pos<slen; i++) {
	if (sscanf( str+pos, "%s", cname) != 1) { clear(); return false; }
	cnames[i] = (char*) calloc( (strlen(cname)+1), sizeof(char) );
	strcpy( cnames[i], cname );  pos += strlen(cname)+1;
      }
      if (i < ncases) { clear(); return false; }
      else if (aidx>=0) setAttributeDisc( name, ncases, cnames );
      else              setOutputLabel( name, ncases, cnames );
    } else {
      if (aidx>=0) setAttributeCont( name, (float)atof(sn1), (float)atof(sn2) );
    }
    return true;
  }
  
};
  
}

#endif  // MLH_ATTRIBUTE_INFO_HPP


// ===================================================================
#if 0	// beginning of the example
// ===================================================================

#include <iostream>
#include "../common/mlh_attribute_info.hpp"
#include "../common/mlh_decision_tree.hpp"
using namespace std;
int main(int argc, char **argv)
{
  if (argc<2) return EXIT_FAILURE;
  MLH::AttributeInfo<double> info;
  info.initAttributes( 4 );
  info.setAttributeDisc( "outlook",     3, "sunny", "overcast", "rain" );
  info.setAttributeCont( "temperature", 0, 100 );
  info.setAttributeCont( "humidity",    0, 100 );
  info.setAttributeDisc( "windy",       2, "false", "true" );
  info.setOutputLabel  ( "Play?",       2, "DoNotPlay", "Play" );
  MLH::DecisionTree<double,bool>	dtree;
  dtree.setAttributeInfo( &info );
  dtree.printInfo();
  dtree.readTrainingDataFile( argv[1], true );
  dtree.train( 1, 5 );
  dtree.printInfo();
  return EXIT_SUCCESS;
}
// Contents of the input text file  --------------
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
