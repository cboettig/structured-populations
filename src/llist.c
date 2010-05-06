#include "llist.h"
 
LLIST *list_add(LLIST **p, double i)
{
	if (p == NULL)
		return NULL;
 
	LLIST *n = malloc(sizeof(LLIST));
	if (n == NULL)
		return NULL;
 
	n->next = *p; /* the previous element (*p) now becomes the "next" element */
	*p = n;       /* add new empty element to the front (head) of the list */
	n->data = i;
 
	return *p;
}
 
void list_remove(LLIST **p) /* remove head */
{
	if (p != NULL && *p != NULL)
	{
		LLIST *n = *p;
		*p = (*p)->next;
		free(n);
	}
}
 
LLIST **list_search(LLIST **n, double i)
{
	if (n == NULL)
		return NULL;
 
	while (*n != NULL)
	{
		if ((*n)->data == i)
		{
			return n;
		}
		n = &(*n)->next;
	}
	return NULL;
}
 
void list_print(LLIST *n)
{
	if (n == NULL)
	{
		printf("list is empty\n");
	}
	while (n != NULL)
	{
		printf("print %p %p %g\n", n, n->next, n->data);
		n = n->next;
	}
}

void list_free(LLIST *n)
{
	if (n == NULL)
	{
		printf("list is empty\n");
	}
	while (n != NULL)
	{
		list_remove(&n);
	}
}
  
