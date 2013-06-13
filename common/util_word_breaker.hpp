
//
// UTIL::WordBreaker class for seperating words in a string
//
// Jaeil Choi
// last modified in June, 2006
//
// For usage, read the example code at the end of this file.
//

#ifndef UTIL_WORD_BREAKER_HPP
#define UTIL_WORD_BREAKER_HPP

#include <iostream>

namespace UTIL {
  
class WordBreaker {
public:
  char **words, *copy;
  int  wcount;		// the number of words
public:
  WordBreaker() : words(NULL), copy(NULL), wcount(0) {}
  WordBreaker(const char *str, char *br=" \t\r\n", char comment='\0') 
    : words(NULL), copy(NULL), wcount(0) { parse( (char*)str, br, true, comment ); }
  ~WordBreaker() { clear(); }
  void clear(void) { if (words) free(words); if (copy) free(copy); words=NULL; copy=NULL; wcount=0; }
  
  int  count(void) { return wcount; }
  char* operator[](int idx) { return (words && idx>=0 && idx<wcount ? words[idx] : (char*)""); }
  //char* operator[](int idx) { return (words && idx<wcount ? words[idx] : NULL); }
  
public:
  int parse(char *str, char *br=" \t\r\n", bool keep=true, char comment='\0') {
    // Break the string 'str' into words, with break character in 'br'.
    // If 'keep' is false, the input string will be overwritten. (It will crash if 'str' is 'const char*'.)
    // Characters 'comment' and following characters will be cut out and ignored.
    clear();
    if (!str) return 0;
    int  i, n, len = (int)strlen(str);
    bool in;
    if (keep) { copy = (char*)malloc((len+1)*sizeof(char));  strcpy(copy, str);  str = copy; }
    for (i = n = 0, in = false; i < len; i++) {			// count the words
      if (str[i] == comment) break;
      else if (isBreaker(str[i], br)) in = false;
      else { if (!in) n++;   in = true; }
    }
    wcount = n;
    words = (char**) calloc( wcount, sizeof(char*) );	// create word pointers
    for (i = n = 0, in = false; i < len; i++) {			// break into words
      //std::cout << "processing '" << str[i] << "'   at i = " << i << std::endl;
      if (str[i] == comment) { str[i] = '\0'; break; }
      else if (isBreaker(str[i], br)) { str[i] = '\0'; in = false; }
      else { if (!in) words[n++] = str+i;  in = true; }
    }
    return wcount;
  }
  
  char* parseFirst(char *str, char *first, char *br=" \t\r\n", char comment='\0') {
    // Copy the first word of the string 'str' in 'first', and return the rest.
    // The 'comment' character and following characters will be cut out and ignored.
    int  i, n = 0, len = (int)strlen(str);
    bool in = false, copied = false;
    char *second = NULL;
    for (i = 0; i < len; i++) {
      if (str[i] == comment) { str[i] = '\0'; break; }	// cut out the comment
      else if (isBreaker(str[i], br)) {
	if (in && !copied) { first[n] = '\0';  copied = true;}
	in = false;
      } else {
	if (!copied) first[n++] = str[i];		// copy the first word
	else if (!in && !second) second = str+i;	// beginning of the second word
	in = true;
      }
    }
    return second;
  }
  
  bool isBreaker(char c, char *breakers=NULL) {
    if (breakers == NULL) breakers = " \t\r\n";
    int i, n = (int)strlen(breakers);
    for (i = 0; i < n; i++) if (c == breakers[i]) break;
    return (i < n ? true : false);
  }
  
  bool isEmpty(char *str, char comment='\0') {
    // Check if the input string 'str' is empty (consists of white characters),
    // Note that the part commented out is not considered.
    for (int i = 0; str[i] && str[i] != comment; i++) 
      if (str[i]!=' ' && str[i]!='\t' && str[i]!='\n') return false;
    return true;
  }
  void trim(char *str) { // remove white characters on both ends
    int len, i, idx_start, idx_end;
    len = (int)strlen(str);
    for (i = 0; str[i]; i++) if (str[i]!=' ' && str[i]!='\t' && str[i]!='\n') break;
    if (str[i]) idx_start = i;
    else { str[0] = '\0'; return; }
    for (i = len-1; i >= 0; i--) if (str[i]!=' ' && str[i]!='\t' && str[i]!='\n') break;
    idx_end = i;
    for (i = 0; i < (idx_end - idx_start + 1); i++) str[i] = str[idx_start + i];
    str[i] = '\0';
  }
  
  int getWordPosition(int widx) {
    return (words && widx>=0 && widx<wcount ? words[widx] - copy : -1);
  }
  
  void printInfo(void) {
    printf("words (%d) : ", wcount);
    for (int i = 0; i < wcount; i++) printf("'%s' ", words[i]);
    printf("\n");
  }
  
};

}	// end of namespace UTIL

#endif  /* UTIL_WORD_BREAKER_HPP */



#if 0 // =============================================================

#include <iostream>
#include "src/common/util_word_breaker.hpp"
using namespace std;

int main(int argc, char *argv[])
{
  UTIL::WordBreaker wb;
  char str1[120];  sprintf(str1, "I am Tom.");
  char str2[120];  sprintf(str2, "I am Tom. \n I'm a student. # comment");
  char str3[120];  sprintf(str3, " Hello, \tWorld. My name is Jerry. After this# comes the comment.");
  
  wb.parse( str1 );
  wb.printInfo();
  wb.parse( str2, " .,'\t\r\n", true, '#' );
  wb.printInfo();
  
  char first[80], *rest;
  rest = wb.parseFirst( str3, first, NULL, '#' );
  cout << "first : " << first << endl;
  cout << "rest  : " << rest << endl;
  
  return EXIT_SUCCESS;
}

#endif // ============================================================
