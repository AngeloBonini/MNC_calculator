#include <math.h>
#include <stdbool.h>
#define MAX 100

void subMatriz_(int ordem, double a[][MAX], int coluna, double sub[][MAX])
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

double determinante_(int ordem, double a[][MAX])
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
            subMatriz_(ordem, a, j, sub);
            s += a[0][j] * pow(-1, j) * determinante_(ordem - 1, sub);
         }
   }
   return s;
}

bool temSubMatrizesNaoSingulares(int ordem, double matriz[][MAX])
{
   for (int i = 1; i <= ordem; i++)
      if (determinante_(i, matriz) == 0)
         return false;
   return true;
}

bool checaPositivaDefinida(int ordem, double matriz[][MAX])
{
   for (int i = 1; i <= ordem; i++)
      if (determinante_(i, matriz) <= 0)
         return false;
   return true;
}

bool checaSimetria(int ordem, double matriz[][MAX])
{
   for (int i = 0; i < ordem; i++)
      for (int j = i + 1; j < ordem; j++)
         if (matriz[i][j] != matriz[j][i])
            return false;
   return true;
}

bool checaDiagonalPrincipal(int ordem, double matriz[][MAX])
{
   for (int i = 0; i < ordem; i++)
      if (matriz[i][i] == 0)
         return false;
   return true;
}
