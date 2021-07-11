
#include <math.h>
#include <stdio.h>
#define MAX 100


void copiaMatriz_GaussJordan(int ordem, double a[][MAX], double copia[][MAX])
{
   for (int i = 0; i < ordem; i++)
      for (int j = 0; j < ordem; j++)
         copia[i][j] = a[i][j];
}
void copiaVetor_gaussJordan(int ordem, double a[], double copia[])
{
   for (int i = 0; i < ordem; i++)
      copia[i] = a[i];
}

void subMatriz_gaussJordan(int ordem, double a[][MAX], int coluna, double sub[][MAX])
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

double determinante_gaussJordan(int ordem, double a[][MAX])
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
            subMatriz_gaussJordan(ordem, a, j, sub);
            s += a[0][j] * pow(-1, j) * determinante_gaussJordan(ordem - 1, sub);
         }
   }
   return s;
}
bool temSubMatrizesNaoSingulares_gaussJordan(int ordem, double matriz[][MAX])
{
   for (int i = 1; i <= ordem; i++)
      if (determinante_gaussJordan(i, matriz) == 0)
         return false;
   return true;
}

bool gaussJordan(int n, double a[][MAX], double b[], double x[]){
   double aL[MAX][MAX], bL[MAX], m;
   copiaMatriz_GaussJordan(n, a, aL);
   copiaVetor_gaussJordan(n, b, bL);

   if (!temSubMatrizesNaoSingulares_gaussJordan(n, a))
      return false;

   for (int k = 0; k < n; k++)
   {
      for (int i = 0; i < n; i++)
         if (i != k)
         {
            m = aL[i][k] / aL[k][k];

            for (int j = 0; j < n; j++)
               aL[i][j] -= m * aL[k][j];

            bL[i] -= m * bL[k];
         }
   }

   for (int i = 0; i < n; i++)
      x[i] = bL[i] / aL[i][i];
   return true;
}