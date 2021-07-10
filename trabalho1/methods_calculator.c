#include "read_print.h"
#include "calc.h"
#include <math.h>
#include <stdbool.h>

// void geraSubMatriz(int n, double a[][MAX], int c, double sub[][MAX])
// {
//    for (int i = 1; i < n; i++)
//       for (int j = 0, k = 0; j < n; j++, k++)
//       {
//          if (j != c)
//             sub[i - 1][k] = a[i][j];
//          else
//             k--;
//       }
// }

// void transpoe(int n, double a[][MAX])
// {
//    double aux;
//    for (int i = 0; i < n; i++)
//       for (int j = i + 1; j < n; j++)
//       {
//          aux = a[i][j];
//          a[i][j] = a[j][i];
//          a[j][i] = aux;
//       }
// }

// void id(int n, double a[][MAX])
// {
//    for (int i = 0; i < n; i++)
//       for (int j = 0; j < n; j++)
//          a[i][j] = (i == j) ? 1 : 0;
// }

// double determinante(int n, double a[][MAX])
// {
//    double s = 0;
//    if (n == 1)
//       return a[0][0];
//    else
//    {
//       double sub[MAX][MAX];
//       for (int j = 0; j < n; j++)
//          if (a[0][j] != 0)
//          {
//             geraSubMatriz(n, a, j, sub);
//             s += a[0][j] * pow(-1, j) * determinante(n - 1, sub);
//          }
//    }
//    return s;
// }

bool temSubMatrizesNaoSingulares(int n, double a[][MAX])
{
   for (int i = 1; i <= n; i++)
      if (determinante(i, a) == 0)
         return false;
   return true;
}

bool ehDefPositiva(int n, double a[][MAX])
{
   for (int i = 1; i <= n; i++)
      if (determinante(i, a) <= 0)
         return false;
   return true;
}
bool ehSimetrica(int n, double a[][MAX])
{
   for (int i = 0; i < n; i++)
      for (int j = i + 1; j < n; j++)
         if (a[i][j] != a[j][i])
            return false;
   return true;
}
void copiaMatriz(int n, double a[][MAX], double copia[][MAX])
{
   for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
         copia[i][j] = a[i][j];
}
void copiaVet(int n, double a[], double copia[])
{
   for (int i = 0; i < n; i++)
      copia[i] = a[i];
}
bool diagPrincipalNaoNula(int n, double a[][MAX])
{
   for (int i = 0; i < n; i++)
      if (a[i][i] == 0)
         return false;
   return true;
}
bool criterioLinhas(int n, double a[][MAX])
{
   for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
         if (j != i && abs(a[i][j] / a[i][i]) >= 1)
            return false;
   return true;
}
bool criterioColunas(int n, double a[][MAX])
{
   for (int j = 0; j < n; j++)
      for (int i = 0; i < n; i++)
         if (i != j && abs(a[i][j] / a[j][j]) >= 1)
            return false;
   return true;
}
bool criterioSassenfeld(int n, double a[][MAX])
{
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
}

double normaInf(int n, double v[])
{
   double max = 0;
   for (int i = 0; i < n; i++)
      if (abs(v[i]) > max)
         max = abs(v[i]);
   return max;
}

double *diferencaVet(int n, double v1[], double v2[])
{
   double *v = (double *)malloc(sizeof(double) * n);
   for (int i = 0; i < n; i++)
      v[i] = v1[i] - v2[i];
   return v;
}

bool sistemaTriangularSuperior(int n, double a[][MAX], double b[], double x[])
{
   if (determinante(n, a) == 0)
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
bool sistemaTriangularInferior(int n, double a[][MAX], double b[], double x[])
{
   if (determinante(n, a) == 0)
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
bool gaussCompacto(int n, double a[][MAX], double b[], double x[])
{
   double u[MAX][MAX], l[MAX][MAX], bL[MAX], s;

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

      s = 0;
      for (int k = 0; k < p; k++)
         s += l[p][k] * bL[k];
      bL[p] = b[p] - s;

      for (int i = p; i < n; i++)
      {
         s = 0;
         for (int k = 0; k < p; k++)
            s += l[i][k] * u[k][p];

         l[i][p] = (a[i][p] - s) / u[p][p];
      }
   }

   sistemaTriangularSuperior(n, u, bL, x); //U.x = bL
   return true;
}
bool cholesky(int n, double a[][MAX], double b[], double x[])
{
   double l[MAX][MAX], y[MAX], s;

   if (!ehSimetrica(n, a) || !ehDefPositiva(n, a))
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

   sistemaTriangularInferior(n, l, b, y);

   transpoe(n, l);
   sistemaTriangularSuperior(n, l, y, x);
   return true;
}
bool gaussJordan(int n, double a[][MAX], double b[], double x[])
{
   double aL[MAX][MAX], bL[MAX], m;
   copiaMatriz(n, a, aL);
   copiaVet(n, b, bL);

   if (!temSubMatrizesNaoSingulares(n, a))
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
bool matrizInversa(int n, double a[][MAX], double x[][MAX])
{
   int op;
   double e[MAX][MAX];

   printf("\n[0] Decomposição LU\n[1] Gauss Compacto\n");
   do
   {
      printf("Determinar a inversa por: ");
      scanf("%d", &op);
   } while (op != 0 && op != 1);
   printf("\n");

   if (!temSubMatrizesNaoSingulares(n, a))
      return false;

   id(n, e);
   for (int i = 0; i < n; i++)
      (op == 0) ? decomposicaoLU(n, a, e[i], x[i]) : gaussCompacto(n, a, e[i], x[i]);

   transpoe(n, x);
   return true;
}

bool jacobi(int n, double a[][MAX], double b[], double e, double x_ant[], int maxIte, double x[], int *ite)
{
   if (!diagPrincipalNaoNula(n, a) || determinante(n, a) == 0 || (!criterioLinhas(n, a) && !criterioColunas(n, a)))
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

      copiaVet(n, x, x_ant);
   }
   (*ite)--;
   return true;
}
bool gaussSeidel(int n, double a[][MAX], double b[], double e, double x_ant[], int maxIte, double x[], int *ite)
{
   if (!diagPrincipalNaoNula(n, a) || determinante(n, a) == 0 || (!criterioLinhas(n, a) && !criterioSassenfeld(n, a)))
      return false;

   double *v, s, x_rec[MAX];

   copiaVet(n, x_ant, x_rec);
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
      if (normaInf(n, v) / normaInf(n, x) < e)
         return true;

      copiaVet(n, x, x_ant);
   }
   (*ite)--;
   return true;
}

int menu()
{
   int op;
   printf("   Cálculo de sistemas lineares   \n");
   printf("         e matriz inversa         \n");
   printf("..................................\n\n");
   printf("| Determinante:\n");
   printf("  01 › Laplace\n\n");
   printf("| Sistema linear:\n");
   printf("  02 › Triangular Inferior\n");
   printf("  03 › Triangular Superior\n\n");
   printf("     Métodos diretos:\n");
   printf("     04 › Decomposição LU\n");
   printf("     05 › Cholesky\n");
   printf("     06 › Gauss Compacto\n");
   printf("     07 › Gauss Jordan\n\n");
   printf("     Métodos iterativos:\n");
   printf("     08 › Jacobi\n");
   printf("     09 › Gauss Seidel\n\n");
   printf("| Matriz:\n");
   printf("  10 › Inversa\n\n");
   do
   {
      printf("Resp.: ");
      scanf("%d", &op);
   } while (op < 1 || op > 10);
   printf("\n");
   return op;
}

int main()
{
   int op, n, maxIte, ite;
   double a[MAX][MAX], i[MAX][MAX], x[MAX], x0[MAX], b[MAX], e;
   char r;
   bool ok;

   do
   {
      system("clear");
      op = menu();

      printf("Ordem da matriz: ");
      scanf("%d", &n);
      printf("\nCoeficientes da matriz %dx%d:\n", n, n);
      read2d_array(n, a);

      if (op == 1)
         printf("\nA solução é: %.4lf\n", determinante(n, a));
      else if (op >= 2 && op <= 7)
      {
         printf("\nTermos independentes:\n");
         readArray(n, b);

         switch (op)
         {
         case 2:
            ok = sistemaTriangularInferior(n, a, b, x);
            break;
         case 3:
            ok = sistemaTriangularSuperior(n, a, b, x);
            break;
         case 4:
            ok = decomposicaoLU(n, a, b, x);
            break;
         case 5:
            ok = cholesky(n, a, b, x);
            break;
         case 6:
            ok = gaussCompacto(n, a, b, x);
            break;
         case 7:
            ok = gaussJordan(n, a, b, x);
            break;
         }
         if (ok)
         {
            printf("\nO vetor solução é:\n");
            printArray(n, x);
         }
         else
            printf("\nO método não converge!\n");
      }
      else if (op == 8 || op == 9)
      {
         printf("\nTermos independentes:\n");
         readArray(n, b);
         printf("\nAproximação inicial para o vetor solução:\n");
         readArray(n, x0);
         printf("\nErro: ");
         scanf("%lf", &e);
         printf("\nNúmero máximo de iterações: ");
         scanf("%d", &maxIte);

         if (op == 8)
            ok = jacobi(n, a, b, e, x0, maxIte, x, &ite);
         else
            ok = gaussSeidel(n, a, b, e, x0, maxIte, x, &ite);

         if (ok)
         {
            printf("\nO vetor solução é:\n");
            printArray(n, x);
         }
         else
            printf("\nO método não converge!\n");
      }
      else
      {
         ok = matrizInversa(n, a, i);
         if (ok)
         {
            printf("A matriz solução é:\n");
            print2d_array(n, i);
         }
         else
            printf("\nO método não converge!\n");
      }

      do
      {
         printf("\nRetornar ao menu [s/n]? ");
         fflush(stdin);
         scanf(" %c", &r);
      } while (r != 'n' && r != 's');
   } while (r == 's');
   return 0;
}