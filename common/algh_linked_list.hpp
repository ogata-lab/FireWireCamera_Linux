
//
// ALGH::LinkedList<T> class template
//
// Jaeil Choi
// last modified in May, 2003
//
// - add, find, and remove a node at any position in the list,
//   using add(), insert(), append(), getNth(), or remove()
//
// Example
//   LinkedList<char> list;
//   list.add('C');
//   list.add('A');
//   list.append('D');
//   list.insert(1, 'B');
//   list.remove(0);
//   cout << "list: ";
//   for (list.goFirst(); (c = list.getCurr()); list.goNext())  cout << c << " ";
//   // OR
//   for (int i = 0; i < list.count; i++) {
//     char c = list.getNth(i);
//     cout << c << " " << flush;
//   }
//   cout << endl;	// result : "list: B C D"
//

#ifndef ALGH_LINKED_LIST_HPP
#define ALGH_LINKED_LIST_HPP

#include <iostream>

namespace ALGH {
  
using namespace std;

template <class T>
class LinkedListNode {
public:
  T		obj;
  LinkedListNode<T>	*next;
};

template <class T>
class LinkedList {
public:
  int		count;
  LinkedListNode<T>	*head, *tail, *current, *previous;
  
public:
  LinkedList() : count(0), head(NULL), tail(NULL), current(NULL), previous(NULL) {}
  ~LinkedList() { clear();  }
  void	clear(void) {
    LinkedListNode<T> *p = head, *next;
    while (p) { next = p->next;  free(p);  p = next; }
    count = 0;
    head = current = previous = NULL;
  }
  void	clearAndFree(void) {
    LinkedListNode<T> *p = head, *next;
    while (p) { 
      if ((void*)(p->obj)) free((void*)(p->obj));
      next = p->next;  free(p);  p = next; 
    }
    count = 0;
    head = current = previous = NULL;
  }
  
  // -----------------------------------------------------------------
  // iterating the list
  // -----------------------------------------------------------------
  void goFirst(void)   { previous = NULL;  current = head; }
  void goLast(void)    { previous = NULL;  current = tail; }
  void goNext(void)    { if (current) { previous = current; current = current->next; } }
  void goNth(int idx)  { current = NthNode(idx); }
  int  goFind(T v) {
    int i;   LinkedListNode<T> *p = head;
    for (i = 0; p && (p->obj != v); i++) p = p->next;
    if (p && p->obj == v) { current = p;    return i;  } // set 'current' and return the index
    else                  { current = NULL; return -1; }
  }    
  T getCurr(void) { return (current ? current->obj : (T)0); }
  T getPrev(void) { return (previous ? previous->obj : (T)0); }
  T getNext(int n=1) { 
    LinkedListNode<T> *np=current;
    if (np==NULL) { std::cerr << "Warning (getNext()) : current==NULL" << std::endl; return (T)0; }
    for (int i=0; i<n && np; i++, np=np->next);
    return ( np ? np->obj : (T)0 );
  }
  T    remCurr(void) {		// remove the current node
    // This function works in "for(goFirst();(p=getCurr());goNext())" loop only
    LinkedListNode<T> *np = current;
    current = np->next;
    if (previous) previous->next = current;
    if      (np == NULL)    return (T)0;
    else if (np == head) {  head = current; if (np == tail) tail = NULL;  }
    else if (np == tail) {  tail = previous;  } 
    T data = np->obj;
    free(np);
    count--;
    return data;
  }
  void insCurr (T obj) { 
    if (current == NULL) { append(obj); return; }
    LinkedListNode<T> *newp = (LinkedListNode<T>*) calloc(1, sizeof(LinkedListNode<T>));
    newp->obj  = obj;
    newp->next = current->next;
    current->next = newp;
  }
  
  // -----------------------------------------------------------------
  // -----------------------------------------------------------------
  
  T    getHead(void) { return (head ? head->obj : (T)0); }
  T    getNth(int index)  {	// starts from 0
    LinkedListNode<T> *p = NthNode(index);
    if (!p) return (T)0;
    return p->obj;
  }
  
  int	add(T obj) {			// add new node at the start
    LinkedListNode<T> *newp = (LinkedListNode<T>*) calloc(1, sizeof(LinkedListNode<T>));
    newp->obj  = obj;
    newp->next = head;
    if (head == NULL) { count = 0; tail = newp; }
    head = newp;
    count++;
    return 0;
  }
  
  int	append(T obj) {			// add new node at the end
    LinkedListNode<T> *newp = (LinkedListNode<T>*) calloc(1, sizeof(LinkedListNode<T>));
    newp->obj  = obj;
    newp->next = NULL;
    if (head == NULL) {
      head = tail = newp;
    } else {
      tail->next = newp;
      tail = newp;
    }
    return count++;
  }
  
  void	insert(int index, T obj) {	// insert node at the given index position
    LinkedListNode<T> *prev, *curr = NthNodeWithPrev(index, &prev);
    if (curr == head) add(obj);
    else {
      LinkedListNode<T> *newp = (LinkedListNode<T>*) calloc(1, sizeof(LinkedListNode<T>));
      newp->obj  = obj;
      newp->next = curr;
      prev->next = newp;
      if (curr == NULL) tail = newp;
      count++;
    }
  }
  
  int  find(T obj) {
    LinkedListNode<T> *p = head;  int i = 0;
    for (p = head; p; p = p->next, i++) if (p->obj == obj) return i;
    return -1;
  }
  
  inline T  pop(void) { return remove(0); }
  
  T  remove(int index) {		// remove a node in the middle
    LinkedListNode<T> *prev, *curr;
    if (index == 0) {
      if ((curr = head) == NULL) return (T)0;
      head = curr->next;
      if (curr == tail) tail = NULL;
    } else {
      if ((curr = NthNodeWithPrev(index, &prev)) == NULL) return (T)0;
      prev->next = curr->next;
      if (curr == tail) tail = prev;
    }
    T data = curr->obj;
    free(curr);
    count--;
    return data;
  }
  
  void findAndRemove(T obj) {
    LinkedListNode<T> *prev = NULL, *curr = head;
    while (curr && curr->obj != obj) { prev = curr; curr = curr->next; }
     if (curr == NULL) return;
    else if (curr == head) head = curr->next;
    else prev->next = curr->next;
    if (curr == tail) tail = prev;
     free(curr);
    count--;
  }    
  
  void  merge(LinkedList<T>& llist, int (*compare)(T a, T b)=NULL) {
    // Merge two list, assuming both list are in appropriate order,
    LinkedListNode<T> *mp = this->head;  // list of existing nodes
    LinkedListNode<T> *lp = llist.head;  // list of given new nodes
    this->head = this->tail = NULL;  this->count = 0;
    llist.head = llist.tail = NULL;  llist.count = 0;
    LinkedListNode<T> *curr = NULL, *next=NULL, *tmp;
    int compare_result;
    while (mp || lp) {
      // decide the next node
      if      (mp == NULL) compare_result = +1; // { next = lp; lp = lp->next; }
      else if (lp == NULL) compare_result = -1; // { next = mp; mp = mp->next; }
      else if (compare == NULL) { 
	if      (mp->obj < lp->obj) compare_result = -1;
	else if (mp->obj > lp->obj) compare_result = +1;
	else                        compare_result =  0;
      } else compare_result = compare(mp->obj, lp->obj);
      // set the next node
      switch (compare_result) {
      case -1: next = mp; mp = mp->next;  break;
      case +1: next = lp; lp = lp->next;  break;
      case  0: next = mp; mp = mp->next;  tmp=lp->next; delete(lp); lp=tmp;  break;
      }
      // add the next node to the list
      next->next = NULL;
      if (curr == NULL) curr = head = tail = next;
      else { curr->next = tail = next; curr = next; }
      count++; 
    }
  }
  
  void	print(void (*ftn)(T), char* comment = NULL) {
    std::cout << "LIST(" << count << "): " << std::flush;
    LinkedListNode<T> *p = head;
    while (p) {
      if (ftn) (ftn)(p->obj);
      else std::cout << p->obj << std::flush;
      std::cout << " ";
      p = p->next;
    }
    if (comment) std::cout << " " << comment;
    std::cout << std::endl;
  }
  
  LinkedListNode<T>* NthNode(int index) {
    LinkedListNode<T> *p = head;
    if (index < 0)  index = count + index;
    for (int i = 0; i < index && p; i++) p = p->next;
    return p;
  }
  
  LinkedListNode<T>* NthNodeWithPrev(int index, LinkedListNode<T>** pp) {
    *pp = NULL;
    LinkedListNode<T> *p = head;
    if (index < 0)  index = count + index;
    for (int i = 0; i < index && p; i++) { *pp = p; p = p->next; }
    return p;
  }
  
  T* getListInArray(T array[]=NULL) {
    if (array == NULL) array = (T*) malloc( count * sizeof(T) );
    int  i = 0;  T obj; 
    for (goFirst(); (obj = getCurr()); goNext()) array[i++] = obj;
    return array;
  }
  
  void swapList(LinkedList<T> *list) {
    LinkedList<T> tmp;
    memcpy( &tmp, this, sizeof(LinkedList<T>) );
    memcpy( this, list, sizeof(LinkedList<T>) );
    memcpy( list, &tmp, sizeof(LinkedList<T>) );
    memset( &tmp,    0, sizeof(LinkedList<T>) );
  }
};


template <class T>
static ostream &operator<< (ostream &os, const LinkedList<T> &lst) {
  os << "[ ";
  LinkedListNode<T> *p = lst.head;
  while (p) { cout << (T)(p->obj) << " ";  p = p->next; }
  os << "]";
  return os;
}

}	// end of namespace ALGH

#endif  /* ALGH_LINKED_LIST_HPP */


