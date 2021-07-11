#include <math.h>
#include <stdio.h>
#define MAX 100


void subMatriz_TriangularSuperior(int ordem, double a[][MAX], int coluna, double sub[][MAX])
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
double determinante_TriangularSuperior(int ordem, double a[][MAX])
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
            subMatriz_TriangularSuperior(ordem, a, j, sub);
            s += a[0][j] * pow(-1, j) * determinante_TriangularSuperior(ordem - 1, sub);
         }
   }
   return s;
}

bool sistemaTriangularSuperior(int n, double a[][MAX], double b[], double x[])
{
   if (determinante_TriangularSuperior(n, a) == 0)
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