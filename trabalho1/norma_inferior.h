#include <stdio.h>
#include<stdbool.h>
#include<math.h>
#define abs_(n) ((n) > 0 ? (n) : -(n))

double normaInf(int n, double v[])
{
   double max = 0;
   for (int i = 0; i < n; i++)
      if (abs_(v[i]) > max)
         max = abs(v[i]);
   return max;
}