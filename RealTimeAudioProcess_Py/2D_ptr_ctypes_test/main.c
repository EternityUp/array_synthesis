#include <stdio.h>
#include <stdlib.h>
#include "print_2d_ptr.h"

int main(int argc, char *argv[])
{
  int **ptr;
  int rows = 2;
  int cols = 3;
  ptr = (int **)malloc(sizeof(int*) * rows);
  for (int i = 0; i < rows; i++)
  {
    ptr[i] = (int*)malloc(sizeof(int) * cols);
    for (int j = 0; j < cols; j++)
      ptr[i][j] = i + j;
  }
  
  print_2d_ptr(ptr, rows, cols);
  
  
  for (int i = 0; i < rows; i++)
    free(ptr[i]);
  
  free(ptr);
  
  return 0;
}