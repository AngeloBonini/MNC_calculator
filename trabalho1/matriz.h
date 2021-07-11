
#include <math.h>
#include <stdio.h>
#define MAX 100

void subMatriz(int ordem, double a[][MAX], int coluna, double sub[][MAX])
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


void transposta(int ordem, double a[][MAX])
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


void identidade(int ordem, double a[][MAX])
{
   for (int i = 0; i < ordem; i++)
      for (int j = 0; j < ordem; j++)
         a[i][j] = (i == j) ? 1 : 0;
}

double determinante(int ordem, double a[][MAX])
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
            subMatriz(ordem, a, j, sub);
            s += a[0][j] * pow(-1, j) * determinante(ordem - 1, sub);
         }
   }
   return s;
}



