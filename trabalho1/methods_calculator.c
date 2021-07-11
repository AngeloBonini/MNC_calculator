#include "read_print.h"
#include "cholesky.h"
#include "gauss_compacto.h"
#include "gauss_jordan.h"
#include "gauss_seidel.h"
#include "jacobi.h"
#include "matriz.h"
#include "testes_matrizes.h"
#include "sistemaTriangularSuperior.h"
#include "sistemaTriangularInferior.h"
#include "decompostaLU.h"
#include <math.h>
#include <stdbool.h>
#include <stdio.h>

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

   identidade(n, e);
   for (int i = 0; i < n; i++)
      (op == 0) ? decomposicaoLU(n, a, e[i], x[i]) : gaussCompacto(n, a, e[i], x[i]);

   transposta(n, x);
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
      system("cls");
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