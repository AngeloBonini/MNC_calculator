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
#include <windows.h>
#include <locale.h>
#include <stdbool.h>
#include <stdio.h>

//Lucas Augusto Ribeiro Massarico
//RA: 191027723

//Angelo Bonini
//RA: 191026077

//Caio Regal
//RA: 181025442

short posX = 48, posY = 9;

void moveXY(int x, int y, char t[10])
{
    HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE);
    COORD position = {x,y};
    SetConsoleCursorPosition(hStdout, position);
    printf("%s", t);   
}

bool matrizInversa(int n, double a[][MAX], double x[][MAX]){
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


int menu(){
   setlocale(LC_ALL, "Portuguese");

   moveXY (47, 0, "_______________________________________");
   moveXY (46, 1, "|    Resolução de Sistemas Lineares     |");
   moveXY (46, 2, "|      e Cálculo de matriz inversa      |");
   moveXY (46, 3, "|                                       |");
   moveXY (46, 4, "|    Lucas Massarico, RA: 191027723     |");
   moveXY (46, 5, "|     Angelo Bonini, RA: 191026077      |");
   moveXY (46, 6, "|       Caio Regal, RA: 181025442       |");
   moveXY (46, 7, "|_______________________________________|");
   moveXY (46, 8, "|                                       |");
   moveXY (46, 9, "|         Calcular Determinante         |");
   moveXY (46, 10, "|                                       |");
   moveXY (46, 11, "|      Sistema Triangular Inferior      |");
   moveXY (46, 12, "|      Sistema Triangular Superior      |");
   moveXY (46, 13, "|                                       |");
   moveXY (46, 14, "|            Decomposição LU            |");
   moveXY (46, 15, "|                Cholesky               |");
   moveXY (46, 16, "|             Gauss Compacto            |");
   moveXY (46, 17, "|              Gauss Jordan             |");
   moveXY (46, 18, "|                                       |");
   moveXY (46, 19, "|                 Jacobi                |");
   moveXY (46, 20, "|              Gauss Seidel             |");
   moveXY (46, 21, "|                                       |");
   moveXY (46, 22, "|             Matriz Inversa            |");
   moveXY (46, 23, "|                                       |");
   moveXY (46, 24, "|                 SAIR                  |");
   moveXY (46, 25, "|_______________________________________|");
}

int main(){
   int op, n, maxIte, ite;
   double a[MAX][MAX], i[MAX][MAX], x[MAX], x0[MAX], b[MAX], e;
   char r, xd, pronto;
   bool ok;

   moveXY(48,9, "->");

   do{
      menu();
      moveXY(posX, posY, "->");
      xd = toupper(getch());

      switch(xd){
            case 'H':
               if (posY > 9){
                  moveXY(48,posY, "  ");
                  if (posY == 24){
                    posY = 22;
                  }
                  else if (posY == 22){
                    posY = 20;
                  }
                  else if (posY == 19){
                    posY = 17;
                  }
                  else if (posY == 14){
                     posY = 12;
                  }
                  else if (posY == 11){
                     posY = 9;
                  }
                  else {
                     posY-=1;
                  }
                  moveXY(48,posY, "->");             
                }
                break;

            case 'P':
               if (posY < 24){
                  moveXY(48,posY, "  ");
                  if (posY == 9){
                     posY = 11;
                  }
                  else if (posY == 12) {
                     posY = 14;
                  }
                  else if (posY == 17){
                     posY = 19;
                  }
                  else if (posY == 20){
                     posY = 22;
                  }
                  else if (posY == 22){
                     posY = 24;
                  }
                  else {
                     posY+=1;
                  }
                  moveXY(48,posY, "->");             
               }
               break;

               case 13:
               	if (posY == 24){
                    system("cls");
                    exit(0);
                } else {
                	system ("cls");
	               	printf("\n\tOrdem da matriz: ");
	                scanf("%d", &n);
	                printf("\n\tCoeficientes da matriz %dx%d:\n", n, n);
	                read2d_array(n, a);
                
	                if (posY == 9){
	                	system ("cls");
	                	printf ("\n\tCalculando o DETERMINANTE");
	                	printf ("\n\n\tMatriz:\n");
	                	print2d_array (n, a);
	                    printf("\tA solução do determinante é: %.4lf\n", determinante(n, a));
	                    printf ("\n\n\tPressione ENTER para voltar ao menu\n\t");
						do{
	                    	pronto = getch();
	                    }while (pronto != 13);
	                    system ("cls");
	                }
	                  
	                if (posY == 11){
	                	system ("cls");
	                	printf ("\n\tSistema Triangular Inferior");
	                	printf ("\n\n\tMatriz:\n");
	                	print2d_array (n, a);
	                    printf("\n\tTermos independentes:\n");
	                    readArray(n, b);
	                    ok = sistemaTriangularInferior(n, a, b, x);
	                    if (ok) {
	                    	printf("\n\tO vetor solução é:\n");
	                        printArray(n, x);
	                    } else
	                    	printf("\nO método não converge!\n");
	                    printf ("\n\tPressione ENTER para voltar ao menu\n\t");
	                    do{
	                    	
	                    	pronto = getch();
	                    }while (pronto != 13);
	                    system ("cls");              
	                }
	                
	                if (posY == 12){
	                    system("cls");
	                    printf ("\n\tSistema Triangular Superior");
	                	printf ("\n\n\tMatriz:\n");
	                	print2d_array (n, a);
	                    printf("\n\tTermos independentes:\n");
	                    readArray(n, b);
	                    ok = sistemaTriangularSuperior(n, a, b, x);
	                    if (ok) {
	                    	printf("\n\tO vetor solução é:\n");
	                        printArray(n, x);
	                    } else
	                    	printf("\nO método não converge!\n");
	                    printf ("\n\tPressione ENTER para voltar ao menu\n\t");
	                    do{
	                    	
	                    	pronto = getch();
	                    }while (pronto != 13);
	                    system ("cls");
	                }
	                
	                if (posY == 14){
	                	system("cls");
	                    printf ("\n\tMétodo Direto: Decomposição LU");
	                	printf ("\n\n\tMatriz:\n");
	                	print2d_array (n, a);
	                    printf("\n\tTermos independentes:\n");
	                    readArray(n, b);
	                    ok = decomposicaoLU(n, a, b, x);
	                    if (ok) {
	                    	printf("\n\tO vetor solução é:\n");
	                        printArray(n, x);
	                    } else
	                    	printf("\nO método não converge!\n");
	                    printf ("\n\tPressione ENTER para voltar ao menu\n\t");
	                    do{
	                    	
	                    	pronto = getch();
	                    }while (pronto != 13);
	                    system ("cls");
					}
					
					if (posY == 15){
	                	system("cls");
	                    printf ("\n\tMétodo direto: Cholesky");
	                	printf ("\n\n\tMatriz:\n");
	                	print2d_array (n, a);
	                    printf("\n\tTermos independentes:\n");
	                    readArray(n, b);
	                    ok = cholesky(n, a, b, x);
	                    if (ok) {
	                    	printf("\n\tO vetor solução é:\n");
	                        printArray(n, x);
	                    } else
	                    	printf("\nO método não converge!\n");
	                    printf ("\n\tPressione ENTER para voltar ao menu\n\t");
	                    do{
	                    	
	                    	pronto = getch();
	                    }while (pronto != 13);
	                    system ("cls");
					}
					
					if (posY == 16){
	                	system("cls");
	                    printf ("\n\tMétodo direto: Gauss Compacto");
	                	printf ("\n\n\tMatriz:\n");
	                	print2d_array (n, a);
	                    printf("\n\tTermos independentes:\n");
	                    readArray(n, b);
	                    ok = gaussCompacto(n, a, b, x);
	                    if (ok) {
	                    	printf("\n\tO vetor solução é:\n");
	                        printArray(n, x);
	                    } else
	                    	printf("\nO método não converge!\n");
	                    printf ("\n\tPressione ENTER para voltar ao menu\n\t");
	                    do{
	                    	
	                    	pronto = getch();
	                    }while (pronto != 13);
	                    system ("cls");
					}
					
					if (posY == 17){
	                	system("cls");
	                    printf ("\n\tMétodo direto: Gauss Jordan");
	                	printf ("\n\n\tMatriz:\n");
	                	print2d_array (n, a);
	                    printf("\n\tTermos independentes:\n");
	                    readArray(n, b);
	                    ok = gaussJordan(n, a, b, x);
	                    if (ok) {
	                    	printf("\n\tO vetor solução é:\n");
	                        printArray(n, x);
	                    } else
	                    	printf("\nO método não converge!\n");
	                    printf ("\n\tPressione ENTER para voltar ao menu\n\t");
	                    do{
	                    	
	                    	pronto = getch();
	                    }while (pronto != 13);
	                    system ("cls");
					}
					
					if (posY == 19){
						system("cls");
	                    printf ("\n\tMétodo iterativo: Jacobi");
	                	printf ("\n\n\tMatriz:\n");
	                	print2d_array (n, a);
	                    printf("\n\tTermos independentes:\n");
	                    readArray(n, b);
	                    printf("\n\tAproximação inicial para o vetor solução:\n");
	                    readArray(n, x0);
	                    printf("\n\tErro: ");
	         			scanf("%lf", &e);
	         			printf("\nNúmero máximo de iterações: ");
	         			scanf("%d", &maxIte);
	         			ok = jacobi(n, a, b, e, x0, maxIte, x, &ite);
	         			if (ok) {
	                    	printf("\n\tO vetor solução é:\n");
	                        printArray(n, x);
	                    } else
	                    	printf("\nO método não converge!\n");
	                    printf ("\n\tPressione ENTER para voltar ao menu\n\t");
	                    do{
	                    	
	                    	pronto = getch();
	                    }while (pronto != 13);
	                    system ("cls");
					}
					
					if (posY == 20){
						system("cls");
	                    printf ("\n\tMétodo iterativo: Gauss Seidel");
	                	printf ("\n\n\tMatriz:\n");
	                	print2d_array (n, a);
	                    printf("\n\tTermos independentes:\n");
	                    readArray(n, b);
	                    printf("\n\tAproximação inicial para o vetor solução:\n");
	                    readArray(n, x0);
	                    printf("\n\tErro: ");
	         			scanf("%lf", &e);
	         			printf("\nNúmero máximo de iterações: ");
	         			scanf("%d", &maxIte);
	         			ok = gaussSeidel(n, a, b, e, x0, maxIte, x, &ite);
	         			if (ok) {
	                    	printf("\n\tO vetor solução é:\n");
	                        printArray(n, x);
	                    } else
	                    	printf("\nO método não converge!\n");
	                    printf ("\n\tPressione ENTER para voltar ao menu\n\t");
	                    do{
	                    	
	                    	pronto = getch();
	                    }while (pronto != 13);
	                    system ("cls");
					}
	
					if (posY == 22){
						system("cls");
	                    printf ("\n\tMatriz Inversa");
	                	printf ("\n\n\tMatriz:\n");
	                	print2d_array (n, a);
				        ok = matrizInversa(n, a, i);
				        if (ok) {
				        	printf("\n\tA matriz solução é:\n");
				            print2d_array(n, i);
				        } else
				        	printf("\n\tO método não converge!\n");
				        printf ("\n\tPressione ENTER para voltar ao menu\n\t");
	                    do{
	                    	
	                    	pronto = getch();
	                    }while (pronto != 13);
	                    system ("cls");
					}
				}
               	
               	
                break;
   }
}while (1);
system("pause");
   return 0;
}
