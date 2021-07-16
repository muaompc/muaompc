#include <stdio.h>
#include <stdlib.h>

struct hh {int i; double d;};

int main()
{
  struct hh **tt;
  struct hh *uu[2];
  struct hh ee[2];
  struct hh *oo;
  struct hh h1 = {1, 1.};
  struct hh h2 = {2, 2.};
  struct hh h3 = {3, 3.};
  struct hh h4 = {4, 4.};

  oo = (struct hh*)calloc(2, sizeof(struct hh));
  uu[0] = &h1;
  uu[1] = &h4;

  ee[0] = h1;
  ee[1] = h4;
  oo[0] = h1;
  oo[1] = h4;

  printf("%f \n", uu[0]->d);
  printf("%f \n", uu[1]->d);
  printf("%f \n", ee[0].d);
  printf("%f \n", ee[1].d);
  printf("%f \n", oo[0].d);
  printf("%f \n", oo[1].d);
}
