
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#define MAX 100

void subMatriz_cholesky(int ordem, double a[][MAX], int coluna, double sub[][MAX])
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
double determinante_cholesky(int ordem, double a[][MAX])
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
            subMatriz_cholesky(ordem, a, j, sub);
            s += a[0][j] * pow(-1, j) * determinante_cholesky(ordem - 1, sub);
         }
   }
   return s;
}

void transposta_cholesky(int ordem, double a[][MAX])
{
   double aux;
   for (int i = 0; i < ordem; i++)
      for (int j = i + 1; j < ordem; j++)
      {
         aux = a[i][j];
         a[i][j] = a[j][i];
         a[j][i] = aux;
      }
}

bool checaSimetria_cholesky(int ordem, double matriz[][MAX])
{
   for (int i = 0; i < ordem; i++)
      for (int j = i + 1; j < ordem; j++)
         if (matriz[i][j] != matriz[j][i])
            return false;
   return true;
}

bool checaPositivaDefinida_cholesky(int ordem, double matriz[][MAX])
{
   for (int i = 1; i <= ordem; i++)
      if (determinante_cholesky(i, matriz) <= 0)
         return false;
   return true;
}
bool sistemaTriangularInferior_cholesky(int n, double a[][MAX], double b[], double x[])
{
   if (determinante_cholesky(n, a) == 0)
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


bool sistemaTriangularSuperior_cholesky(int n, double a[][MAX], double b[], double x[])
{
   if (determinante_cholesky(n, a) == 0)
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

bool cholesky(int n, double a[][MAX], double b[], double x[]){
   double l[MAX][MAX], y[MAX], s;

   if (!checaSimetria_cholesky(n, a) || !checaPositivaDefinida_cholesky(n, a))
      return false;

   for (int j = 0; j < n; j++)
   {

      s = 0;
      for (int k = 0; k < j; k++)
         s += pow(l[j][k], 2);

      l[j][j] = sqrt(a[j][j] - s);

      for (int i = j + 1; i < n; i++)
      {
         s = 0;
         for (int k = 0; k < j; k++)
            s += l[i][k] * l[j][k];

         l[i][j] = (a[i][j] - s) / l[j][j];
      }
   }

   sistemaTriangularInferior_cholesky(n, l, b, y);

   transposta_cholesky(n, l);
   sistemaTriangularSuperior_cholesky(n, l, y, x);
   return true;
};