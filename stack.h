/*
 * stack.c
 * stack module summary:
 *
 * module that contains a struct of a simple stack implemented by linked list.
 *
 *
 * initialize				- creates a new stack
 * push						- pushes a node into the stack
 * pop						- pops a node from the stack
 * isEmpty					- checks if the stack is empty
 */

#ifndef STACK_H_
#define STACK_H_

#include "errorChecking.h"
typedef struct node node;

struct node
{
	void* data;
	node* next;
};

/*returns 1 if empty stack, else 0 */
int isEmpty(node *s);

/*pushes data into the stack, returns new head */
node* push(node* head,int* data);

/*writes top of stack into element, pops its out of stack, and returns new head */
node* pop(node *head,int** element);

/* initializes empty stack */
void initialize(node** s);

#endif
