
#include <math.h>
#include <stdbool.h>
#define MAX 100
#define abs_(n) ((n) > 0 ? (n) : -(n));

bool criterioSassenfeld(int n, double a[][MAX]){
   double beta[MAX], max = 0;

   for (int i = 0; i < n; i++)
   {
      beta[i] = 0;
      for (int j = 0; j < i - 1; j++)
         beta[i] += abs(a[i][j] / a[i][i]) * beta[j];

      for (int j = i + 1; j < n; j++)
         beta[i] += abs(a[i][j] / a[i][i]);

      if (beta[i] >= 1)
         return false;
   }
   return true;
};

bool criterioLinhas(int n, double a[][MAX]){
   for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
         if (j != i && abs(a[i][j] / a[i][i]) >= 1)
            return false;
   return true;
};
double *diferencaVet(int n, double v1[], double v2[]){
   double *v = (double *)malloc(sizeof(double) * n);
   for (int i = 0; i < n; i++)
      v[i] = v1[i] - v2[i];
   return v;
};

void subMatriz_gaussSeidel(int ordem, double a[][MAX], int coluna, double sub[][MAX])
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

double determinante_gaussSeidel(int ordem, double a[][MAX])
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
            subMatriz_gaussSeidel(ordem, a, j, sub);
            s += a[0][j] * pow(-1, j) * determinante_gaussSeidel(ordem - 1, sub);
         }
   }
   return s;
}

bool checaDiagonalPrincipal_gaussSeidel(int ordem, double matriz[][MAX])
{
   for (int i = 0; i < ordem; i++)
      if (matriz[i][i] == 0)
         return false;
   return true;
}

void copiaVetor_gaussSeidel(int ordem, double a[], double copia[])
{
   for (int i = 0; i < ordem; i++)
      copia[i] = a[i];
}
double normaInf_gaussSeidel(int n, double v[])
{
   double max = 0;
   for (int i = 0; i < n; i++)
      if (abs_(v[i]) > max)
         max = abs_(v[i]);
   return max;
}
bool gaussSeidel(int n, double a[][MAX], double b[], double e, double x_ant[], int maxIte, double x[], int *ite)
{
   if (!checaDiagonalPrincipal_gaussSeidel(n, a) || determinante_gaussSeidel(n, a) == 0 || (!criterioLinhas(n, a) && !criterioSassenfeld(n, a)))
      return false;

   double *v, s, x_rec[MAX];

   copiaVetor_gaussSeidel(n, x_ant, x_rec);
   for (*ite = 1; *ite <= maxIte; (*ite)++)
   {
      for (int i = 0; i < n; i++)
      {
         s = 0;
         for (int j = 0; j < n; j++)
            if (j != i)
               s += a[i][j] * x_rec[j];

         x[i] = (b[i] - s) / a[i][i];
         x_rec[i] = x[i];
      }

      v = diferencaVet(n, x, x_ant);
      if (normaInf_gaussSeidel(n, v) / normaInf_gaussSeidel(n, x) < e)
         return true;

      copiaVetor_gaussSeidel(n, x, x_ant);
   }
   (*ite)--;
   return true;
}