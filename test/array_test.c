/** Test routine for arrays
 */

#include <stdio.h>
#include "array.h"

#define ARRTEST_BLOCKSIZ 1024
#define ARRTEST_MAXNUM 34600

int main () 
{
  int *arr1, *arr2, *hp;
  int al = 0;
  int i, *ip;

  printf("1st array:\n");
  ARRCREATE(arr1, ARRTEST_BLOCKSIZ);
  al=(int) ARRLEN(arr1);
  printf("[%i]\n", al);
  i=0;
  while ((al=(int) ARRLEN(arr1)) < ARRTEST_MAXNUM) { 
    ARRCURR(arr1) = i;
    printf("[%i]\n", al);
    ARRNEXTP(hp, arr1);
    if (hp == NULL) 
      return ERRCODE_NOMEM;
    al=(int) ARRLEN(arr1);
    i++;
  }

  for (i=0; i<(int)ARRLEN(arr1); i++) 
    printf("[%i] %i\n", i, arr1[i]);

  ARRDELETE(arr1);

  printf("2nd array:\n");
  ARRCREATE(arr2, ARRTEST_BLOCKSIZ);
  al=(int) ARRLEN(arr2);
  printf("[%i]\n", al);
  i=0;
  while ((al=(int) ARRLEN(arr2)) < ARRTEST_MAXNUM) { 
    ARRNEXTP(ip, arr2);
    if (!ip)
      return ERRCODE_NOMEM;
    *ip = i++;
    al=(int) ARRLEN(arr2);
  }

  for (i=0; i<(int)ARRLEN(arr2); i++) 
    printf("[%i] %i\n", i, arr2[i]);

  ARRDELETE(arr2);
  return ERRCODE_SUCCESS;
}
