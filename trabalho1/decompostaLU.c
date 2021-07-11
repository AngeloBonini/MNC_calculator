#include "read_print.h"
#include <math.h>
#include <stdbool.h>
#include "matriz.h"
#include "testes_matrizes.h"
#include "norma_inferior.h"
#include "sistemaTriangularSuperior.h"
#include "sistemaTriangularInferior.h"



bool decomposicaoLU(int n, double a[][MAX], double b[], double x[])
{
   double u[MAX][MAX], l[MAX][MAX], y[MAX], s;

   if (!temSubMatrizesNaoSingulares(n, a))
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

   sistemaTriangularInferior(n, l, b, y);
   sistemaTriangularSuperior(n, u, y, x);
   return true;
}