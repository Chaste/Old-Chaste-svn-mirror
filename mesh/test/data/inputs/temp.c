#include <stdio.h>
int main(int argc, char *argv[])
{
  int i, j,n;
  n=0;
  for (i=0;i<=20;i++)
    {
    for (j=0;j<=20;j++)
      {
      printf("%i\t%f\t%f\n", n++, i/200.0, j/200.0);
    }
  }
}
