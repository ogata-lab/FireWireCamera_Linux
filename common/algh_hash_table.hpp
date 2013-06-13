
//
// ALGH::HashTable<K_t,D_t> class template
// 
// Jaeil Choi
// last modified in July, 2004
//
//
// Usage : 
//   HashTable<KEY_T,DATA_T> htable (size, hash_ftn);
//   htable.insert(key1, data1);
//   htable.lookup(key1, &data);
//
// You can use default hash funtions;
//   "int hash_integer(int)"    for 'int'   key types
//   "int hash_pointer(void*)"  for 'void*' key types
// For other user-defined key types, you need to 
//   a. define the key class, and make everything 'public'
//   b. define constructors and '!=' operator in the class
//   c. define the hash function of the key
//   d. If the key class cannot be copied, define '=' operator for it.
//
// Implementation details :
//   Every hash node has (key, data) pair.
//   Collisions are resolved by chaining, i.e. different key may
//     have same key value, and they are put in the linked list.
//   Default hash function is from
//     http://www.concentric.net/~Ttwang/tech/inthash.htm
//
// For sample codes, look at examples at the end of this file.
// 

#ifndef ALGH_HASH_TABLE_H
#define ALGH_HASH_TABLE_H

#include <iostream>
#include <iomanip>
#include <cstddef>

namespace ALGH {
  
using namespace std;

static int hash_integer(int key) 
{
  // default hash function for 'int' key
  // Usage : HashTable<int,CData> htable( 128, hash_integer );
  key += (key << 12);
  key ^= (key >> 22);
  key += (key << 4);
  key ^= (key >> 9);
  key += (key << 10);
  key ^= (key >> 2);
  key += (key << 7);
  key ^= (key >> 12);
  return key;
}

static int hash_pointer(void* key) 
{
  // default hash function for pointer key
  // Usage : HashTable<void*,CData> htable( 128, hash_pointer );
  int value = ((int)key) / 4;
  value += (value << 12);
  value ^= (value >> 22);
  value += (value << 4);
  value ^= (value >> 9);
  value += (value << 10);
  value ^= (value >> 2);
  value += (value << 7);
  value ^= (value >> 12);
  return value;
}

static int hash_string(char *key)
{
  int value = 0, shift = 4; // multiplied by 16
  while(*key) { value = (value << shift) + *key;  key++; }
  return (int) value;
}


template <class K_t, class D_t>
class HashNode {
public:
  K_t	key;
  D_t	data;
  HashNode<K_t,D_t> *next;
public:
  HashNode(K_t key, D_t data) { this->key = key; this->data = data; }
  ~HashNode() {}
};


template <class K_t, class D_t>
class HashTable {
private:
  HashNode<K_t,D_t>	**nodes;
  int			table_size;	// always  2^n
  int			entry_count;
  
  int  (*hash_ftn)(K_t key);
  
public:
  HashTable() : nodes(NULL), hash_ftn(NULL) {}
  HashTable(int size, int (*hash_ftn) (K_t key)) : nodes(NULL), hash_ftn(NULL) { create(size, hash_ftn); }
  ~HashTable() { clear(); }
  void clear(void) {
    HashNode<K_t,D_t> *np, *next;
    if (nodes) {	// clear data nodes
      for (int i = 0; i < table_size; i++)
	for (np = nodes[i]; np; np = next) { next = np->next; delete(np); }
      delete nodes;
    }
    nodes = NULL;
    this->hash_ftn = NULL;
  }
  
  int  size (void) { return (nodes ? table_size  : 0); }
  int  count(void) { return (nodes ? entry_count : 0); }
  
  void create(int table_size, int  (*hash_ftn) (K_t key)) {
    clear();
    // decide the table_size as power of 2
    for (this->table_size = 8; this->table_size < table_size; this->table_size *= 2);
    // create hash table
    nodes = new HashNode<K_t,D_t>*[this->table_size];
    for (int i = 0; i < this->table_size; i++) nodes[i] = NULL;
    // set functions
    this->hash_ftn = hash_ftn;
    if (hash_ftn == NULL) cerr << "Error : HashTable<> cannot have NULL hash function." << endl;
    entry_count = 0;
  }
    
  bool lookup(K_t key, D_t* datap) {
    HashNode<K_t,D_t> *np;
    np = nodes[ hash_ftn(key) & (table_size-1) ];
    while (np && np->key != key) np = np->next;
    if (np) {
      *datap = np->data;
      return true;
    } else {
      return false;
    }
  }
  
  void insert(K_t key, D_t data, bool AllowDuplicate = true) {
    HashNode<K_t,D_t> *np;
    int pos = hash_ftn(key) & (table_size-1);
    
    if (!AllowDuplicate) {  // check if it's already inserted
      np = nodes[ pos];
      while (np && np->key != key) np = np->next;
      if (np) np->data = data;
    }
      
    np = new HashNode<K_t,D_t>(key, data);
    if (np == NULL) cout << "Error : HashTable<> cannot allocate a hash node ! " << endl;
    np->next = nodes[pos];
    nodes[pos] = np;
    entry_count++;
  }
  
  void remove(K_t key) {
    int pos;
    pos = hash_ftn(key) & (table_size-1);
    HashNode<K_t,D_t> *np, *prev = NULL;
//     prev = (HashNode<K_t,D_t>*) ((char*)(nodes + pos) - (size_t)&((HashNode<K_t,D_t>*)0)->next);
//     np = nodes[pos];
//     while (np && np->key != key) { prev = np; np = np->next; }
//     // remove the node
//     if (np) { prev->next = np->next; delete(np); }
    np = nodes[pos];
    while (np && np->key != key) { prev = np; np = np->next; }
    if (np) {	// remove the node 
      if (np == nodes[pos]) nodes[pos] = np->next;
      else                  prev->next = np->next;
      delete(np); 
    }
    entry_count--;
  }
  
  void initialize(void) {
    // clear all the nodes in the table
    HashNode<K_t,D_t> *np, *next;
    if (nodes == NULL) return;
    for (int i = 0; i < table_size; i++) {
      for (np = nodes[i]; np; np = next) { next = np->next; delete(np); }
      nodes[i] = NULL;
    }
    entry_count = 0;
  }
  
  void printInfo(bool detail = false) {
    int    i, ecount = 0, max = 0, empty = 0;
    double sum = 0.0, sum2 = 0.0, avg;
    
    cout << "HashTable (size=" << table_size << ")" << endl;
    if (nodes == NULL) return;
    for (i = 0; i < table_size; i++) {
      HashNode<K_t,D_t> *np = nodes[i];
      ecount = 0;
      if (detail) cout << "  " << setfill('0') << setw(3) << i << " : ";
      while (np) { 
	if (detail) cout << np->key << "(" << np->data << ") ";  
	ecount++;  np = np->next; 
      }
      if (detail) cout << endl;

      sum   += ecount;
      sum2  += ecount * ecount;
      if (ecount > max) max = ecount;
      if (ecount == 0) empty++;
    }

    avg = sum / table_size;
    printf("  Load factor : %.1f(%.0f/%d)  (var=%.1f max=%d  %d empty slots)\n",
	   avg, sum, table_size, sum2/table_size - avg*avg, max, empty);
  }
  
  // -----------------------------------------------------------------
  // iterating through all the entries
  // -----------------------------------------------------------------
private:
  int			go_index;
  HashNode<K_t,D_t>	*go_node;
public:
  void goFirst(void) {		// go to the first entry
    for (int i = 0; i < table_size; i++)  
      if (nodes[i]) { go_index = i;  go_node = nodes[i];  return; }
    go_index = 0;  go_node = NULL;
  }
  void goNext(void) {		// move on to the next entry
    if (go_node == NULL) return;
    else if (go_node->next) { 
      go_node = go_node->next; 
      return; 
    } else {
      for (int i = go_index+1; i < table_size; i++)  
	if (nodes[i]) { go_index = i;  go_node = nodes[i];  return; }
    }
    go_index = 0;  go_node = NULL;
  }
  D_t  getCurr(K_t *key=NULL) {	// return the data of current entry
    if (go_node == NULL) return NULL;
    if (key) *key = go_node->key; 
    return go_node->data;
  }
  void remCurr(void) {		// remove current entry, and move on
    if (go_node == NULL) return;
    HashNode<K_t,D_t> *curr = go_node;
    goNext();  // the next entry become the 'current'
    remove( curr->key );
  }
  void dummy(void) {
    if (false && hash_ftn == NULL) {
      this->hash_ftn = hash_integer;
      this->hash_ftn = hash_pointer;
      this->hash_ftn = hash_string;
    }
  }
  
};


}	// end of namespace ALGH

#endif // ALGH_HASH_TABLE_H

//   unsigned int Hash32(unsigned int key) {  // by Robert Jenkins
//     key += (key << 12);
//     key ^= (key >> 22);
//     key += (key << 4);
//     key ^= (key >> 9);
//     key += (key << 10);
//     key ^= (key >> 2);
//     key += (key << 7);
//     key ^= (key >> 12);
//     return key;
//   }
//   int Hash32(int key) {		// by Thomas Wang
//     key += ~(key << 15);
//     key ^=  (key >>> 10);
//     key +=  (key << 3);
//     key ^=  (key >>> 6);
//     key += ~(key << 11);
//     key ^=  (key >>> 16);
//     return key;
//   }
//   long Hash64(long key) {	// by Thomas Wang
//     key += ~(key << 32);
//     key ^= (key >>> 22);
//     key += ~(key << 13);
//     key ^= (key >>> 8);
//     key += (key << 3);
//     key ^= (key >>> 15);
//     key += ~(key << 27);
//     key ^= (key >>> 31);
//     return key;
//   }
//   long Hash64_2(long key) {
//     long c1 = 0x6e5ea73858134343L;
//     long c2 = 0xb34e8f99a2ec9ef5L;
//     key ^= ((c1 ^ key) >>> 32);
//     key *= c1;
//     key ^= ((c2 ^ key) >>> 31);
//     key *= c2;
//     key ^= ((c1 ^ key) >>> 32);
//     return key;
//   }
  

// -------------------------------------------------------------------

#if 0
#include <iostream>
#include "algh_hash_table.hpp"
using namespace std;

// sample code for integer key ---------------------------------------

int main(void)
{
  ALGH::HashTable<int,int> ht(5, hash_integer);
  
  ht.insert(101, 1);
  ht.insert(102, 2);
  ht.insert(103, 3);
  ht.insert(104, 4);
  ht.insert(105, 42);
//   ht.remove(105);

  ht.printInfo();
  
  int v;
  if (ht.lookup(101, &v)) printf("lookup '101' : '%d'\n", v);
  if (ht.lookup(102, &v)) printf("lookup '102' : '%d'\n", v);
  if (ht.lookup(104, &v)) printf("lookup '104' : '%d'\n", v);
  if (ht.lookup(105, &v)) printf("lookup '105' : '%d'\n", v);

  return EXIT_SUCCESS;
}

// sample code for multiple integer keys -----------------------------

class HashKey { 
public: 
  int k1, k2; 
  HashKey() {}
  HashKey(int k1, int k2)     { this->k1 = k1; this->k2 = k2; }
  bool operator!= (HashKey k) { return (k1 != k.k1 || k2 != k.k2); }
};
static int hash_ftn(HashKey key) 
{ 
  return hash_integer(key.k1) + hash_integer(key.k2); 
}
  
int main(void)
{
  ALGH::HashTable<HashKey,char> ht2(5, hash_ftn);
  
  ht2.insert( HashKey(101, 11), 'a' );
  ht2.insert( HashKey(102, 12), 'b' );
  ht2.insert( HashKey(103, 13), 'c' );
  ht2.insert( HashKey(104, 14), 'd' );
  ht2.insert( HashKey(105, 15), 'e' );
  ht2.remove( HashKey(15, 105) );

  ht2.printInfo();
  
  char c;
  if (ht2.lookup(HashKey(101, 11), &c)) printf("lookup '101,11' : '%c'\n", c);
  if (ht2.lookup(HashKey(12, 102), &c)) printf("lookup '102,12' : '%c'\n", c);
  if (ht2.lookup(HashKey(104, 14), &c)) printf("lookup '104,14' : '%c'\n", c);
  if (ht2.lookup(HashKey(105, 15), &c)) printf("lookup '105,15' : '%c'\n", c);

  return EXIT_SUCCESS;
}

#endif

// -------------------------------------------------------------------
