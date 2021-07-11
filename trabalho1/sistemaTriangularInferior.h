#include <math.h>
#include <stdio.h>
#include <stdbool.h>
#define MAX 100

void subMatriz_TriangularInferior(int ordem, double a[][MAX], int coluna, double sub[][MAX])
{
   for (int i = 1; i < ordem; i++)
      for (int j = 0, k = 0; j < ordem; j++, k++)
      {
         if (j != coluna)
            sub[i - 1][k] = a[i][j];
         else
            k--;
      }
}
double determinante_TriangularInferior(int ordem, double a[][MAX])
{
   double s = 0;
   if (ordem == 1)
      return a[0][0];
   else
   {
      double sub[MAX][MAX];
      for (int j = 0; j < ordem; j++)
         if (a[0][j] != 0)
         {
            subMatriz_TriangularInferior(ordem, a, j, sub);
            s += a[0][j] * pow(-1, j) * determinante_TriangularInferior(ordem - 1, sub);
         }
   }
   return s;
}
bool sistemaTriangularInferior(int n, double a[][MAX], double b[], double x[])
{
   if (determinante_TriangularInferior(n, a) == 0)
      return false;

   for (int i = 0; i < n; i++)
   {
      double s = 0;
      for (int j = 0; j < i; j++)
         s += a[i][j] * x[j];

      x[i] = (b[i] - s) / a[i][i];
   }
   return true;
}