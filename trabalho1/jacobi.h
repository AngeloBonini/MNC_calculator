#include <math.h>
#include <stdio.h>
#include<stdbool.h>
#define MAX 100
#define abs_(n) ((n) > 0 ? (n) : -(n))

void subMatriz_jacobi(int ordem, double a[][MAX], int coluna, double sub[][MAX])
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

double determinante_jacobi(int ordem, double a[][MAX])
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
            subMatriz_jacobi(ordem, a, j, sub);
            s += a[0][j] * pow(-1, j) * determinante_jacobi(ordem - 1, sub);
         }
   }
   return s;
}
void copiaVetor_jacobi(int ordem, double a[], double copia[])
{
   for (int i = 0; i < ordem; i++)
      copia[i] = a[i];
}

double *diferencaVet(int n, double v1[], double v2[]){
   double *v = (double *)malloc(sizeof(double) * n);
   for (int i = 0; i < n; i++)
      v[i] = v1[i] - v2[i];
   return v;
};

bool checaDiagonalPrincipal_jacobi(int ordem, double matriz[][MAX])
{
   for (int i = 0; i < ordem; i++)
      if (matriz[i][i] == 0)
         return false;
   return true;
}
bool criterioLinhas_jacobi(int n, double a[][MAX]){
   for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
         if (j != i && abs(a[i][j] / a[i][i]) >= 1)
            return false;
   return true;
};
bool criterioColunas_jacobi(int n, double a[][MAX]){
   for (int j = 0; j < n; j++)
      for (int i = 0; i < n; i++)
         if (i != j && abs(a[i][j] / a[j][j]) >= 1)
            return false;
   return true;
};

double normaInf(int n, double v[])
{
   double max = 0;
   for (int i = 0; i < n; i++)
      if (abs_(v[i]) > max)
         max = abs_(v[i]);
   return max;
}

bool jacobi(int n, double a[][MAX], double b[], double e, double x_ant[], int maxIte, double x[], int *ite)
{
   if (!checaDiagonalPrincipal_jacobi(n, a) || determinante_jacobi(n, a) == 0 || (!criterioLinhas_jacobi(n, a) && !criterioColunas_jacobi(n, a)))
      return false;

   double *v, s;
   for (*ite = 1; *ite <= maxIte; (*ite)++)
   {
      for (int i = 0; i < n; i++)
      {
         s = 0;
         for (int j = 0; j < n; j++)
            if (j != i)
               s += a[i][j] * x_ant[j];

         x[i] = (b[i] - s) / a[i][i];
      }

      v = diferencaVet(n, x, x_ant);
      if (normaInf(n, v) / normaInf(n, x) < e)
         return true;

      copiaVetor_jacobi(n, x, x_ant);
   }
   (*ite)--;
   return true;
}