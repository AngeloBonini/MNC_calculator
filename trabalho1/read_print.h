#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX 100
#define abs(n) ((n) > 0 ? (n) : -(n))

#define erro(vExato, vAprox) ((abs(vExato) > 1) ? (abs(vExato - vAprox) / abs(vExato)) : (abs(vExato - vAprox)))

void printArray(int n, double b[])
{
   for (int i = 0; i < n; i++)
      printf("\t%10.4lf\n", b[i]);
   printf("\n");
}

void print2d_array(int n, double a[][MAX])
{
   for (int i = 0; i < n; i++)
   {
    for (int j = 0; j < n; j++){
        printf("\t%10.4lf ", a[i][j]);
	}
      	
      printf("\n");
   }
   printf("\n");
}

void readArray(int n, double v[])
{
   for (int i = 0; i < n; i++){
   		printf ("\t");
		scanf("%lf", &v[i]);
   }
      
}

void read2d_array(int n, double a[][MAX])
{
   for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++){
      	printf ("\tElemento[%d][%d]: ", i, j);
        scanf("%lf", &a[i][j]);
	  }
    	
}

void copy2d_Array(int ordem, double a[][MAX], double copia[][MAX])
{
   for (int i = 0; i < ordem; i++)
      for (int j = 0; j < ordem; j++)
         copia[i][j] = a[i][j];
}

void copyArray(int ordem, double a[], double copia[])
{
   for (int i = 0; i < ordem; i++)
      copia[i] = a[i];
}
