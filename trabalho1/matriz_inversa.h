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