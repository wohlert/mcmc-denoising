#include <stdlib.h>
#include <stdio.h>

int main(int argc, const char *argv[])
{
  for(int i = 0; i < 10000; i++) {
    printf("%f\n", drand48());
  }

  return 0;
}
