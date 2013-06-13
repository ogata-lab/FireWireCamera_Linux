
//
// UTIL::FileName class
//
// Jaeil Choi
// last modified in June, 2003
//

#ifndef UTIL_FILE_NAME_HPP
#define UTIL_FILE_NAME_HPP

#include <iostream>

namespace UTIL {
  
class FileName {
public:
  FileName()  {}
  ~FileName() {}
  
  // -----------------------------------------------------------------
public:
  bool IsEmpty(char *s) {	// 
    for (int i = 0; s[i]; i++) 
      if (s[i]!=' ' && s[i]!='\t' && s[i]!='\n') return false;
    return true;
  }
  void Trim(char *s) {		// remove white characters on both ends
    int len, i, idx_start, idx_end;
    len = (int)strlen(s);
    for (i = 0; s[i]; i++) if (s[i]!=' ' && s[i]!='\t' && s[i]!='\n') break;
    if (s[i]) idx_start = i;
    else { s[0] = '\0'; return; }
    for (i = len-1; i >= 0; i--) if (s[i]!=' ' && s[i]!='\t' && s[i]!='\n') break;
    idx_end = i;
    for (i = 0; i < (idx_end - idx_start + 1); i++) s[i] = s[idx_start + i];
    s[i] = '\0';
  }
  
  // -----------------------------------------------------------------
public:  
  
  char *Ext(char *filename) {
    // find file extension
    if (filename == NULL) return NULL;
    int len, i;
    len = (int)strlen(filename);
    for (i = len-1; i >= 0 && filename[i] != '.'; i--);
    if (i >= 0)  return  filename + i + 1;
    else return filename + len;
  }
  
  char* Base(char *filename, char *buffer=NULL) {
    if (filename == NULL) return NULL;
    if (buffer) strcpy(buffer, filename); else buffer = filename;
    char *ext = Ext(buffer);
    if (ext != NULL) *(ext-1) = '\0';
    return buffer;
  }
  
  char *ReplaceExt(char *filename, char *new_ext, char *buffer=NULL) {
    if (filename == NULL) return NULL;
    if (buffer) strcpy(buffer, filename); else buffer = filename;
    char *ext = Ext(buffer);
    if (ext != NULL) sprintf( ext, "%s", new_ext );
    return buffer;
  }
  
  bool CompareExt(char *filename, char *ext, bool cs=false) {
    if (cs) {
      char *file_ext = Ext(filename);
      LowerCase(ext); if (strcmp(file_ext, ext) == 0) return true;
      UpperCase(ext); if (strcmp(file_ext, ext) == 0) return true;
      return false;
    } else {
      return ( strcmp( Ext(filename), ext ) == 0 );
    }
  }
  
  char *LowerCase(char *filename, char *buffer=NULL) {
    if (filename == NULL) return NULL;
    if (buffer) strcpy(buffer, filename); else buffer = filename;
    int i, len = (int)strlen(buffer);
    for (i = 0; i < len; i++) 
      if (buffer[i] >= 'A' && buffer[i] <= 'Z') 
	buffer[i] = 'a' + (buffer[i] - 'A');
    return buffer;
  }
  
  char *UpperCase(char *filename, char *buffer=NULL) {
    if (filename == NULL) return NULL;
    if (buffer) strcpy(buffer, filename); else buffer = filename;
    int i, len = (int)strlen(buffer);
    for (i = 0; i < len; i++) 
      if (buffer[i] >= 'a' && buffer[i] <= 'z') 
	buffer[i] = 'A' + (buffer[i] - 'a');
    return buffer;
  }
  
  void ReplaceNumber(char *s_old, char *s_new, int number, 
		     bool all=false, char number_chr='#') {
    char format[10], temp[70], *p, *op = s_old, *np = s_new;
    int  format_size, actual_size, i;
    while (true) {
      // find '###' in the string 's_old'
      p = strchr(op, number_chr);
      if (p == NULL) break;
      // copy up to the location 'p'
      while (*op && op < p) *np++ = *op++;
      // replace "###" with specific number
      for (format_size = 0; p[format_size] == number_chr; format_size++);
      sprintf(format, "%%0%dd", format_size);
      sprintf(temp, format, number);
      actual_size = (int)strlen(temp);
      for (i = 0; i < actual_size; i++) *np++ = temp[i];
      op += format_size;
      if (all != true) break;
    }
    // copy remaining characters
    while (*op) *np++ = *op++;
    *np = '\0';
  }
};

}	// end of namespace UTIL

#endif  // UTIL_FILE_NAME_HPP

