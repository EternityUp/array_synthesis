#include "print_2d_ptr.h"
#include <stdio.h>

void print_1d_ptr(int *ptr, int len)
{
  for (int i = 0; i< len; i++)
    printf("%d	", ptr[i]);
  return;
}



void print_2d_ptr(int **ptr, int rows, int cols)
{
  for (int i = 0; i < rows; i++)
  {
    printf("print %dth row\n", i);
    print_1d_ptr(ptr[i], cols);
    printf("\n");
  }
  return;
}