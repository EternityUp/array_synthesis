/* File: example.h */

#define PI 3.14159
#define VERSION "1.0"

enum Beverage { ALE, LAGER, STOUT, PILSNER };


struct Vector {
  double x, y, z;
};

struct Foo {
  int xx = 10;
  char *path = "/home/xzc";
  float yy = 1.0f;
};

struct Bar {
    int  x[16];
};






int fact(int n);

