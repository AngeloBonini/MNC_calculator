#include <stdio.h>
#include<stdbool.h>
#include<math.h>
#define abs_normaInf(n) ((n) > 0 ? (n) : -(n))

double normaInf(int n, double v[])
{
   double max = 0;
   for (int i = 0; i < n; i++)
      if (abs_normaInf(v[i]) > max)
         max = abs_normaInf(v[i]);
   return max;
}