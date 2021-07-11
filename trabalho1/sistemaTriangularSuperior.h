#define MAX 100

#include "matriz.h"

bool sistemaTriangularSuperior(int n, double a[][MAX], double b[], double x[])
{
   if (determinante(n, a) == 0)
      return false;

   for (int i = n - 1; i >= 0; i--)
   {
      double s = 0;
      for (int j = i + 1; j < n; j++)
         s += a[i][j] * x[j];

      x[i] = (b[i] - s) / a[i][i];
   }
   return true;
}