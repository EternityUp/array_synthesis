/* File: example.i */
%module example

%{
#define SWIG_FILE_WITH_INIT
#include "example.h"
%}

%inline %{
extern int My_variable;
extern double density;
%}

%{
extern char *path;
%}
%immutable;
extern char *path;
%mutable;

#define PI 3.14159
#define VERSION "1.0"
enum Beverage { ALE, LAGER, STOUT, PILSNER };

%constant int FOO = 42;

struct Vector {
  double x, y, z;
};

struct Foo {
  %immutable;
  int xx;
  char *path;
  %mutable;
  float yy;
};

struct Bar {
    int  x[16];
};











int fact(int n);
