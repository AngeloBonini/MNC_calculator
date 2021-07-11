#include <stdio.h>
#include<stdbool.h>
#include<math.h>
#define MAX 100

void subMatriz_LU(int ordem, double a[][MAX], int coluna, double sub[][MAX])
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

double determinante_LU(int ordem, double a[][MAX])
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
            subMatriz_LU(ordem, a, j, sub);
            s += a[0][j] * pow(-1, j) * determinante_LU(ordem - 1, sub);
         }
   }
   return s;
}
bool temSubMatrizesNaoSingulares_LU(int ordem, double matriz[][MAX])
{
   for (int i = 1; i <= ordem; i++)
      if (determinante_LU(i, matriz) == 0)
         return false;
   return true;
}
bool sistemaTriangularSuperior_LU(int n, double a[][MAX], double b[], double x[])
{
   if (determinante_LU(n, a) == 0)
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
bool sistemaTriangularInferior_LU(int n, double a[][MAX], double b[], double x[])
{
   if (determinante_LU(n, a) == 0)
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

bool decomposicaoLU(int n, double a[][MAX], double b[], double x[])
{
   double u[MAX][MAX], l[MAX][MAX], y[MAX], s;

   if (!temSubMatrizesNaoSingulares_LU(n, a))
      return false;

   for (int p = 0; p < n; p++)
   {

      for (int j = p; j < n; j++)
      {
         s = 0;
         for (int k = 0; k < p; k++)
            s += l[p][k] * u[k][j];

         u[p][j] = a[p][j] - s;
      }

      for (int i = p; i < n; i++)
      {
         s = 0;
         for (int k = 0; k < p; k++)
            s += l[i][k] * u[k][p];

         l[i][p] = (a[i][p] - s) / u[p][p];
      }
   }

   sistemaTriangularInferior_LU(n, l, b, y);
   sistemaTriangularSuperior_LU(n, u, y, x);
   return true;
}