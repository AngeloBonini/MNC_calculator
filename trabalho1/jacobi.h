bool jacobi(int n, double a[][MAX], double b[], double e, double x_ant[], int maxIte, double x[], int *ite)
{
   if (!checaDiagonalPrincipal(n, a) || determinante(n, a) == 0 || (!criterioLinhas(n, a) && !criterioColunas(n, a)))
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

      copyArray(n, x, x_ant);
   }
   (*ite)--;
   return true;
}