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

int menu() {
   setlocale(LC_ALL, "Portuguese");
   
   moveXY (47, 0, "_______________________________________");
   moveXY (46, 1, "|                                       |");
   moveXY (44, 2, "  |   Resolução de Sistemas Lineares    |");
   moveXY (46, 3, "|     e Cálculo de matriz inversa      |");
   moveXY (46, 4, "|_______________________________________|");
   moveXY (56, 6, "Calcular Determinante");
   moveXY (50, 8, "*********************************");
   moveXY (50, 9, "* ----- Sistemas Lineares ----- *");
   moveXY (50, 10, "*                               *");
   moveXY (50, 11, "*      Triangular Inferior      *");
   moveXY (50, 12, "*      Triangular Superior      *");
   moveXY (50, 13, "*///////////////////////////////*");
   moveXY (50, 14, "*   ---- Métodos diretos ----  *");
   moveXY (50, 15, "*                               *");
   moveXY (50, 16, "*       Decomposição LU       *");
   moveXY (50, 17, "*           Cholesky            *");
   moveXY (50, 18, "*         Gauss Compacto        *");
   moveXY (50, 19, "*          Gauss Jordan         *");
   moveXY (50, 20, "*///////////////////////////////*");
   moveXY (50, 21, "* ---- Métodos iterativos ---- *");
   moveXY (50, 22, "*            Jacobi             *");
   moveXY (50, 23, "*         Gauss Seidel          *");
   moveXY (50, 24, "*********************************");
   moveXY (59, 26, "Matriz Inversa");
   moveXY (63, 28, "EXIT");
}

int main()
{
   int op, n, maxIte, ite;
   double a[MAX][MAX], i[MAX][MAX], x[MAX], x0[MAX], b[MAX], e;
   char r, xd, pronto;
   bool ok;

   moveXY(53,6, "->");
   
   do {
      menu();
      moveXY(posX, posY, "->");
      xd = toupper(getch());
      
      switch(xd){
            case 'H':
                if (posY > 6){
                    moveXY(53,posY, "  ");
                    if (posY == 28){
                       posY = 26;
                    }
                    else if (posY == 26){
                       posY = 23;
                    }
                    else if (posY == 22){
                       posY = 19;
                    }
                    else if (posY == 16){
                       posY = 12;
                    }
                    else if (posY == 11){
                       posY = 6;
                    }
                    else {
                       posY-=1;
                    }
                    moveXY(53,posY, "->");             
                    
                }
            break;
           
            case 'P':
                if (posY < 28){
                    moveXY(53,posY, "  ");
                    if (posY == 6){
                       posY = 11;
                    }
                    else if (posY == 12) {
                       posY = 16;
                    }
                    else if (posY == 19){
                       posY = 22;
                    }
                    else if (posY == 23){
                       posY = 26;
                    }
                    else if (posY == 26){
                       posY = 28;
                    }
                    else {
                       posY+=1;
                    }
                    moveXY(53,posY, "->");             
                }
            break;
           
            case 13:
                if (posY == 6){
                  system("cls");
                  moveXY(49,1, ":::::::::Fácil:::::::::: ");
                  printf("\nOrdem da matriz: ");
                  scanf("%d", &n);
                  printf("\nCoeficientes da matriz %dx%d:\n", n, n);
                  leMatriz(n, a);
                  printf("\nA solução é: %.4lf\n", determinante(n, a));
                  do{
                        pronto = getch();
                  }while (pronto != 13);
                  system ("cls");
                }
               
                if (posY == 11){
                  system("cls");
                  printf("\nOrdem da matriz: ");
                  scanf("%d", &n);
                  printf("\nCoeficientes da matriz %dx%d:\n", n, n);
                  leMatriz(n, a);
                  printf("\nTermos independentes:\n");
                  leVetor(n, b);
                  ok = sistemaTriangularInferior(n, a, b, x);
                  if (ok) {
                     printf("\nO vetor solução é:\n");
                     impVetor(n, x);
                  } else
                     printf("\nO método não converge!\n");
                  do{
                     pronto = getch();
                  }while (pronto != 13);
                  system ("cls");              
                }
                
                if (posY == 15){
                    system("cls");
                    moveXY(49,1, "TESTE");
                    
                    
                    
                }
               
                if (posY == 28){
                    system("cls");
                    exit(0);
                 }
            break;
        }
        }while(1);
   system("pause");
   return 0;
}