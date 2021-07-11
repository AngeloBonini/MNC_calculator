#include "read_print.h"
#include <math.h>
#include <stdbool.h>


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