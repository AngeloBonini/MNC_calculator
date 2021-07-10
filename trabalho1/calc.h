/**
 * @file calc.h
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2021-07-10
 * 
 * @copyright Copyright (c) 2021
 * transposed
 * @methods: transposta, identidade, submatriz, determinante  
 */
#include "read_print.h"

void geraSubMatriz(int n, double a[][MAX], int c, double sub[][MAX])
{
   for (int i = 1; i < n; i++)
      for (int j = 0, k = 0; j < n; j++, k++)
      {
         if (j != c)
            sub[i - 1][k] = a[i][j];
         else
            k--;
      }
}


void transpoe(int n, double a[][MAX])
{
   double aux;
   for (int i = 0; i < n; i++)
      for (int j = i + 1; j < n; j++)
      {
         aux = a[i][j];
         a[i][j] = a[j][i];
         a[j][i] = aux;
      }
}


void id(int n, double a[][MAX])
{
   for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
         a[i][j] = (i == j) ? 1 : 0;
}

double determinante(int n, double a[][MAX])
{
   double s = 0;
   if (n == 1)
      return a[0][0];
   else
   {
      double sub[MAX][MAX];
      for (int j = 0; j < n; j++)
         if (a[0][j] != 0)
         {
            geraSubMatriz(n, a, j, sub);
            s += a[0][j] * pow(-1, j) * determinante(n - 1, sub);
         }
   }
   return s;
}



