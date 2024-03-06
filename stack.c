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

#include "stack.h"
#include <stdio.h>
#include <stdlib.h>

/*initialize the empty stack */
void initialize(node** head)
{
    *head = NULL;
}

/* push node to stack */
node* push(node* head,int* data)
{
    node* temp;
    temp = malloc(sizeof(node));
    checkNull(temp);
    temp->data = data;
    temp->next = head;
    head = temp;
    return head;
}

/* pop node from stack, returns next element */
node* pop(node *head,int **element)
{
    node* temp;
    temp = head;
    *element = head->data;
    head = head->next;
    free(temp);
    return head;
}

/*returns 1 if the stack is empty, 0 else */
int isEmpty(node* head)
{
    return (head == NULL ? 1 : 0);
}

