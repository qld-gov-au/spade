int a[256], b[256], c[256];

int foo (int n, int x) {
   int i;

   /* feature: support for unknown loop bound  */
   /* feature: support for loop invariants  */
   for (i=0; i<n; i++){
      b[i] = x;
   }

   /* feature: general loop exit condition  */
   /* feature: support for bitwise operations  */
   while (n--){
      a[i] = b[i]&c[i]; i++;
   }

   return(0);
}

int main(void)
{

  foo(256,256);

  return(0);
}
      
