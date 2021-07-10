#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX 100
#define abs(n) ((n) > 0 ? (n) : -(n))

#define erro(vExato, vAprox) ((abs(vExato) > 1) ? (abs(vExato - vAprox) / abs(vExato)) : (abs(vExato - vAprox)))

void printArray(int n, double b[])
{
   for (int i = 0; i < n; i++)
      printf("%10.4lf\n", b[i]);
   printf("\n");
}

void print2d_array(int n, double a[][MAX])
{
   for (int i = 0; i < n; i++)
   {
      for (int j = 0; j < n; j++)
         printf("%10.4lf ", a[i][j]);
      printf("\n");
   }
   printf("\n");
}

void readArray(int n, double v[])
{
   for (int i = 0; i < n; i++)
      scanf("%lf", &v[i]);
}

void read2d_array(int n, double a[][MAX])
{
   for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
         scanf("%lf", &a[i][j]);
}