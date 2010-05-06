#include <stdio.h>    /* for printf */
#include <stdlib.h>   /* for malloc */
 
typedef struct node {
	double data;
	struct node *next; /* pointer to next element in list */
} LLIST;
 
LLIST *list_add(LLIST **p, double i);
void list_remove(LLIST **p);
LLIST **list_search(LLIST **n, double i);
void list_print(LLIST *n);
void list_free(LLIST *n);

